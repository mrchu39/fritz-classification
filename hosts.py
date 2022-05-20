import astropy.units as u
import imp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import requests
import sys, os
import warnings
import ztfiaenv.ztfiaenv as ztfiaenv

from astropy import coordinates
from astropy.utils.data import get_pkg_data_filename
from astropy.utils.exceptions import AstropyWarning
from astropy.io import fits
from astropy.table import Table, vstack
from astropy.visualization import PercentileInterval, AsinhStretch
from astropy.wcs import WCS
from astroquery.ned import Ned
from astroquery.sdss import SDSS
from contextlib import contextmanager
from pprint import pprint
from ztfiaenv.ztfiaenv.functions import get_DLR
from ztfiaenv.ztfiaenv import get_host_from_cat
from ztfquery import bts

from func import *

warnings.simplefilter('ignore', category=AstropyWarning)
warnings.filterwarnings('ignore')

@contextmanager
def suppress_stdout():

    """ Info : Suppresses output from host association
    """

    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:
            yield
        finally:
            sys.stdout = old_stdout

def comment_sublink(source):

    ''' Info : This is not related to hosts, but in the general algorithm we also comment a link for users to upload classifications to TNS along with the
            host association.
        Input : Source
        Returns : None
    '''

    comment_infos = get_source_api(source)['comments']

    for i in range (len(comment_infos)):

        comment_info = comment_infos[i]
        comment = comment_info['text']

        if 'Submit classification to TNS:' in comment:
            print(source + ' already has TNS link.')
            return

    resp = post_comment(source, 'Submit classification to TNS: http://gayatri.caltech.edu:88/query/tns/'+source)

    if resp['status'] == 'success':
        print(bcolors.OKGREEN + source + ' TNS link upload successful.' + bcolors.ENDC)
        return
    else:
        print(bcolors.FAIL + source + ' TNS link upload failed.' + bcolors.ENDC)
        print(bcolors.FAIL + json.dumps(resp, indent=2) + bcolors.ENDC)
        return

def cross_ref_NED(ra, dec):

    ''' Info : Indentifies potential hosts in NED within 3 arcsec of given coordinates
        Input : RA (deg), dec (deg)
        Returns : Host name, host RA, host dec, host type, host redshift
    '''

    co = coordinates.SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg), frame='icrs')

    while True:
        try:
            result_table = Ned.query_region(co, radius=3*u.arcsec)
            break
        except requests.exceptions.ConnectTimeout:
            continue

    # Gets rid of any objects if they're stars, or have 'SN' or 'AT' in they're name (i.e. they're SNe)
    non_stars = np.array([True if ('*' not in result_table['Type'][t] and result_table['Type'][t] != 'SN' and
        'SN' not in result_table['Object Name'][t] and 'AT' not in result_table['Object Name'][t]) else False for t in range(len(result_table['Type']))])

    result_table = result_table[non_stars]

    if len(result_table) == 0:
        return None, None, None, None, None

    try:
        result_table.sort('Separation')
    except ValueError:
        result_table.sort('Distance (arcmin)')

    if str(result_table['Redshift'][0]) != '--':
        # If the redshift source is SPEC or SED, it is from spectra and is of high enough quality
        if str(result_table['Redshift Flag'][0]) == 'SPEC' or str(result_table['Redshift Flag'][0]) == 'SED':
            redshift = float(result_table['Redshift'][0])
        # If the redshift source is PHOT, it is photometric and we ignore it
        elif str(result_table['Redshift Flag'][0]) == 'PHOT':
            redshift = None
        # If the source is unclear, we see if spectral data exists in NED
        else:
            try:
                r = requests.get('http://ned.ipac.caltech.edu/cgi-bin/NEDspectra?objname=' + result_table['Object Name'][0].replace('+', '%2B').replace(' ', '+')
                    + '&extend=never&detail=0&preview=1&numpp=20&refcode=ANY&bandpass=ANY&line=ANY&hconst=73.0&omegam=0.27&omegav=0.73&corr_z=1')

                if 'Spectral data in NED archive for object' in r.text:
                    redshift = float(result_table['Redshift'][0])
                else:
                    redshift = None

            except:
                redshift = None
    else:
        redshift = None

    return result_table['Object Name'][0], result_table['RA'][0], result_table['DEC'][0], result_table['Type'][0], redshift

def cross_ref_SDSS(ra, dec):

    ''' Info : Indentifies potential hosts in SDSS within 3 arcsec of given coordinates
        Input : RA (deg), dec (deg)
        Returns : Host name, host RA, host dec, host type, host redshift (always None, SDSS does not store redshifts for objects)
    '''

    co = coordinates.SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg), frame='icrs')

    while True:
        try:
            result_table = SDSS.query_region(co, radius=3*u.arcsec)
            break
        except requests.exceptions.ReadTimeout:
            continue

    if result_table == None:
        return None, None, None, None, None

    closest = np.argmin(np.sqrt((ra-np.array(result_table['ra']))**2+(dec-np.array(result_table['dec']))**2))

    c = coordinates.SkyCoord(result_table['ra'][closest]*u.degree, result_table['dec'][closest]*u.degree)

    name = 'SDSS ' + c.to_string('hmsdms', precision=2).replace('h', '').replace('m', '').replace('s', '').replace('d', '').replace(' ', '')

    return name, result_table['ra'][closest], result_table['dec'][closest], 'SDSS', None

def get_host_info(ztfname):

    ''' Info : Determines host using ztfiaenv and queries NED for info.
        Input : Name of source
        Returns : Host name, host RA, host dec, host type (e.g. gal, IRS, UVS, etc.), host redshift (if available)
    '''

    source_info = get_source_api(ztfname)

    snra = source_info['ra']
    sndec = source_info['dec']

    # Runs ztfiaenv host association
    data = ztfiaenv.GetHost(ztfname, verbose=False)

    it = 0
    while it < 3:
        try:
            data.get_info(skip_marshal_info=True)
            #with suppress_stdout():
            data.get_host_cat()

            break
        except requests.exceptions.ConnectionError:
            print(it)
            it += 1
            pass
        except (NameError, IndexError, ValueError):
            print(bcolors.FAIL + 'No host in area.' + bcolors.ENDC)
            it = 3

    if it == 3:
        return None, None, None, None, None

    try:
        conra = float(data.host_cat.best_cat.raMean)
        condec = float(data.host_cat.best_cat.decMean)
    except AttributeError:
        try:
            conra = float(data.host_cat.best_cat.ra)
            condec = float(data.host_cat.best_cat.dec)
        except AttributeError:
            return None, None, None, None, None

    # Queries NED for potential host within 3 arcsec of coordinates determined by ztfiaenv
    hostname, hostra, hostdec, hosttype, redshift = cross_ref_NED(conra, condec)

    # If nothing within 3 arcsec in NED, check SDSS
    if hostname == None:
        hostname, hostra, hostdec, hosttype, redshift = cross_ref_SDSS(conra, condec)

    if hostname == None:
        print(bcolors.FAIL + 'No host in area.' + bcolors.ENDC)
        return None, None, None, None, None

    # Run plotting
    size = int(np.round(3*3600*4*np.max([np.abs(snra-hostra), np.abs(sndec-hostdec)])))

    if size < 60:
        size = 60

    while True:
        try:
            with suppress_stdout():
                fitsurl = geturl(hostra, hostdec, size=size, filters="i", format="fits")
            break
        except FileNotFoundError:
            continue

    fh = fits.open(fitsurl[0])
    fim = fh[0].data
    wcs = WCS(fh[0].header)
    fim[np.isnan(fim)] = 0.0
    transform = AsinhStretch() + PercentileInterval(99.5)
    bfim = transform(fim)

    ax = plt.subplot(projection=wcs)

    plt.imshow(bfim,cmap="gray",origin="lower")

    plt.scatter(snra, sndec, transform=ax.get_transform('world'), c='red', marker='+', label='SNe')
    plt.scatter(hostra, hostdec, transform=ax.get_transform('world'), c='blue', label='Host')

    plt.legend(loc='best')

    plt.savefig('test_host.png')

    return hostname, hostra, hostdec, hosttype, redshift

def getimages(ra,dec,size=240,filters="grizy"):

    """ Info : Query ps1filenames.py service to get a list of images
        Input : ra (deg), dec (deg), size (in pixels, 0.25 arcsec/pixel), filters
        Returns : Table with the results
    """

    service = "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
    url = ("{service}?ra={ra}&dec={dec}&size={size}&format=fits"
           "&filters={filters}").format(**locals())
    table = Table.read(url, format='ascii')
    return table

def geturl(ra, dec, size=240, output_size=None, filters="grizy", format="jpg", color=False):

    """ Info : Get URL for images in the table
        Input : ra (deg), dec (deg), size (in pixels, 0.25 arcsec/pixel), filters, format (options are "jpg", "png" or "fits"), color
        Returns : URL as string
    """

    if color and format == "fits":
        raise ValueError("color images are available only for jpg or png formats")
    if format not in ("jpg","png","fits"):
        raise ValueError("format must be one of jpg, png, fits")
    table = getimages(ra,dec,size=size,filters=filters)
    url = ("https://ps1images.stsci.edu/cgi-bin/fitscut.cgi?"
           "ra={ra}&dec={dec}&size={size}&format={format}").format(**locals())
    if output_size:
        url = url + "&output_size={}".format(output_size)
    # sort filters from red to blue
    flist = ["yzirg".find(x) for x in table['filter']]
    table = table[np.argsort(flist)]
    if color:
        if len(table) > 3:
            # pick 3 filters
            table = table[[0,len(table)//2,len(table)-1]]
        for i, param in enumerate(["red","green","blue"]):
            url = url + "&{}={}".format(param,table['filename'][i])
    else:
        urlbase = url + "&red="
        url = []
        for filename in table['filename']:
            url.append(urlbase+filename)
    return url

def post_host(source):

    ''' Info : Identifies potential host of source, then posts information with plot on source page (or updates it if one is already there and
            different host ID'd)
        Input : Source name
        Returns : None
    '''

    flag = 0
    while flag == 0:
        try:
            hostname, hostra, hostdec, hosttype, redshift = get_host_info(source)
            flag = 1
        except requests.exceptions.ReadTimeout:
            continue

    if hostname == None:
        return

    comment_infos = get_source_api(source)['comments']

    flag = 0

    for i in range (len(comment_infos)):

        comment_info = comment_infos[i]
        comment = comment_info['text']

        if 'potential host:' in comment:
            if comment.split(':')[1].split(',')[0].strip() != hostname:
                if redshift != None:
                    resp = edit_comment(source, comment_info['id'], comment_info['author_id'], 'potential host: '+hostname+', ra = '+str(hostra)+
                        ', dec = '+str(hostdec)+', z = '+str(redshift)+', type = '+hosttype+'. host page: http://gayatri.caltech.edu:88/query/host/'+source, 'test_host.png', source+'_host.png')
                else:
                    resp = edit_comment(source, comment_info['id'], comment_info['author_id'], 'potential host: '+hostname+', ra = '+str(hostra)+
                        ', dec = '+str(hostdec)+', type = '+hosttype+'. host page: http://gayatri.caltech.edu:88/query/host/'+source, 'test_host.png', source+'_host.png')

                if resp['status'] == 'success':
                    print(bcolors.OKGREEN + source + ' host association update successful.' + bcolors.ENDC)
                    return
                else:
                    print(bcolors.FAIL + source + ' host association update failed.' + bcolors.ENDC)
                    print(bcolors.FAIL + json.dumps(resp, indent=2) + bcolors.ENDC)
                    return
            else:
                print(source + ' already has an associated host.')
                return

    if redshift != None:
        resp = post_comment(source, 'potential host: '+hostname+', ra = '+str(hostra)+', dec = '+str(hostdec)+', z = '+str(redshift)+
            ', type = '+hosttype+'. host page: http://gayatri.caltech.edu:88/query/host/'+source, 'test_host.png', source+'_host.png')
    else:
        resp = post_comment(source, 'potential host: '+hostname+', ra = '+str(hostra)+', dec = '+str(hostdec)+
            ', type = '+hosttype+'. host page: http://gayatri.caltech.edu:88/query/host/'+source, 'test_host.png', source+'_host.png')

    if resp['status'] == 'success':
        print(bcolors.OKGREEN + source + ' host association upload successful.' + bcolors.ENDC)
        return
    else:
        print(bcolors.FAIL + source + ' host association upload failed.' + bcolors.ENDC)
        print(bcolors.FAIL + json.dumps(resp, indent=2) + bcolors.ENDC)
        return

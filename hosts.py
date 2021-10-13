import numpy as np
import pandas as pd
import sys
import os
#
import requests
import warnings
import imp
import warnings
from astroquery.ned import Ned
from astroquery.sdss import SDSS
import astropy.units as u
from astropy import coordinates
from astropy.utils.exceptions import AstropyWarning
from astropy.table import Table, vstack
from astropy.io import fits
from astropy.visualization import PercentileInterval, AsinhStretch
from astropy.wcs import WCS
from astropy.utils.data import get_pkg_data_filename
warnings.simplefilter('ignore', category=AstropyWarning)
warnings.filterwarnings('ignore')
from pprint import pprint
import matplotlib.pyplot as plt
#
from ztfquery import bts
#pprint(vars(ztfiaenv.ztfiaenv))
from ztfiaenv.ztfiaenv.functions import get_DLR
from ztfiaenv.ztfiaenv import get_host_from_cat
import ztfiaenv.ztfiaenv.ztfiaenv as ztfiaenv
from contextlib import contextmanager
import sys, os
from pprint import pprint

from func import *

@contextmanager
def suppress_stdout():
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:
            yield
        finally:
            sys.stdout = old_stdout

def getimages(ra,dec,size=240,filters="grizy"):

    """Query ps1filenames.py service to get a list of images

    ra, dec = position in degrees
    size = image size in pixels (0.25 arcsec/pixel)
    filters = string with filters to include
    Returns a table with the results
    """

    service = "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
    url = ("{service}?ra={ra}&dec={dec}&size={size}&format=fits"
           "&filters={filters}").format(**locals())
    table = Table.read(url, format='ascii')
    return table

def geturl(ra, dec, size=240, output_size=None, filters="grizy", format="jpg", color=False):

    """Get URL for images in the table

    ra, dec = position in degrees
    size = extracted image size in pixels (0.25 arcsec/pixel)
    output_size = output (display) image size in pixels (default = size).
                  output_size has no effect for fits format images.
    filters = string with filters to include
    format = data format (options are "jpg", "png" or "fits")
    color = if True, creates a color image (only for jpg or png format).
            Default is return a list of URLs for single-filter grayscale images.
    Returns a string with the URL
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

def get_host_info(ztfname):

    source_info = get_source_api(ztfname)

    snra = source_info['ra']
    sndec = source_info['dec']

    data = ztfiaenv.GetHost(ztfname, verbose=False)

    it = 0
    while it < 3:
        try:
            data.get_info(skip_marshal_info=True)
            break
        except requests.exceptions.ConnectionError:
            i += 1
            pass

    if it == 3:
        return None, None, None, None, None

    with suppress_stdout():
        data.get_host_cat()

    try:
        conra = float(data.host_cat.best_cat.raMean)
        condec = float(data.host_cat.best_cat.decMean)
    except AttributeError:
        conra = float(data.host_cat.best_cat.ra)
        condec = float(data.host_cat.best_cat.dec)

    hostname, hostra, hostdec, hosttype, redshift = cross_ref_NED(conra, condec)

    if hostname == None:
        hostname, hostra, hostdec, hosttype, redshift = cross_ref_SDSS(conra, condec)

    if hostname == None:
        print(bcolors.FAIL + 'No host in area.' + bcolors.ENDC)
        return None, None, None, None, None

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

def post_host(source):

    hostname, hostra, hostdec, hosttype, redshift = get_host_info(source)

    if hostname == None:
        return

    comment_infos = get_source_api(source)['comments']

    flag = 0

    for i in range (len(get_source_api(source)['comments'])):

        comment_info = comment_infos[i]
        comment = comment_info['text']

        if 'potential host' in comment:
            #if comment.split(':')[1].split(',')[0].strip() != hostname:
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
            #else:
            #    print(source + ' already has an associated host.')
            #    break

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

def cross_ref_NED(ra, dec):
    co = coordinates.SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg), frame='icrs')
    result_table = Ned.query_region(co, radius=3*u.arcsec)

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
        if str(result_table['Redshift Flag'][0]) == 'SPEC' or str(result_table['Redshift Flag'][0]) == 'SED':
            redshift = float(result_table['Redshift'][0])
        elif str(result_table['Redshift Flag'][0]) == 'PHOT':
            redshift = None
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
    co = coordinates.SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg), frame='icrs')
    result_table = SDSS.query_region(co, radius=30*u.arcsec)

    if result_table == None:
        return None, None, None, None, None

    closest = np.argmin(np.sqrt((ra-np.array(result_table['ra']))**2+(dec-np.array(result_table['dec']))**2))

    c = coordinates.SkyCoord(result_table['ra'][closest]*u.degree, result_table['dec'][closest]*u.degree)

    name = 'SDSS ' + c.to_string('hmsdms', precision=2).replace('h', '').replace('m', '').replace('s', '').replace('d', '').replace(' ', '')

    return name, result_table['ra'][closest], result_table['dec'][closest], 'SDSS', None

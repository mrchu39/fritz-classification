import matplotlib.pyplot as plt
import numpy as np
import os
import shlex
import shutil
import sncosmo
import subprocess
import sys
import time

from astropy.io import ascii
from astropy.table import QTable
from func import *

# These numbers come from running model fits on ~500 Type Ia supernovae
z = 0.060574239946858476
z_std = 0.023994157056121096
x0 = -0.14238796934437334
x0_std = 1.4557579021314682
c = 0.08928354223298558
c_std = 0.15670291093588692
x1 = 0.0007648532623426458
x1_std = 0.0004363803462578883

with open('info.info', 'r') as f:
    SNID_loc = f.read().split('\n')[0].split(':')[1].strip()

def get_photometry(ztfname, format='flux'):

    ''' Info : Retrieves photometry data for a source from Fritz and filters out Nonetype points
        Input : Source name and brightness format ("flux" or "mag")
        Returns : Astropy QTable with data that feeds into sncosmo.fit_lc
    '''

    url = BASEURL+'api/sources/'+ztfname+'/photometry' # Access photometry

    if format == 'flux':
        data = {"format": "flux"}
    elif format == 'mag':
        data = {'format': 'mag'}

    response = api('GET', url, params=data, timeout=10)

    if format == 'flux':

        flux = []
        fluxerr = []
        band = []
        mjd = []
        zpsys = []
        zp = []

        for d in response['data']:
            if d['flux'] != None and (d['filter'] == 'ztfg' or d['filter'] == 'ztfr'):
                flux.append(d['flux'])
                fluxerr.append(d['fluxerr'])
                band.append(d['filter'])
                mjd.append(d['mjd'])
                zpsys.append(d['magsys'])
                zp.append(d['zp'])

        return QTable([mjd, band, flux, fluxerr, zp, zpsys], names=('mjd', 'filter', 'flux','fluxerr', 'zp', 'zpsys'))

    elif format == 'mag':

        mag = []
        magerr = []
        band = []
        mjd = []
        zpsys = []

        for d in response['data']:
            if d['mag'] != None and (d['filter'] == 'ztfg' or d['filter'] == 'ztfr'):
                mag.append(d['mag'])
                magerr.append(d['magerr'])
                band.append(d['filter'])
                mjd.append(d['mjd'])
                zpsys.append(d['magsys'])

        return QTable([mjd, band, mag, magerr, zpsys], names=('mjd', 'filter', 'mag','magerr', 'zpsys'))

def model_lc(source):

    ''' Info : Fits photometry data to light curve using sncosmo.
        Input : source
        Returns : photometry data, fitted parameters, plottable model
    '''

    data = get_photometry(source)

    with open('test_table.txt', 'w') as f:
        f.write(str(data))

    red, red_err = get_redshift(source, True)
    model = sncosmo.Model(source='salt2')

    if red != 'No redshift found':

        if red_err != 'No redshift error found': # If both redshift and error present, use error as bounds to fit redshift
            result, fitted_model = sncosmo.fit_lc(
                data, model,
                ['z', 't0', 'x0', 'x1', 'c'],  # parameters of model to vary
                bounds={'z':(red-red_err, red+red_err)}, minsnr=5)  # bounds on parameters (if any)
        else: # If no redshift error, don't determine redshift
            model.set(z=red)

            result, fitted_model = sncosmo.fit_lc(
                data, model,
                ['t0', 'x0', 'x1', 'c'],
                guess_z=False, minsnr=5)

    else:
        result, fitted_model = sncosmo.fit_lc(
            data, model,
            ['z', 't0', 'x0', 'x1', 'c'],
            bounds={'z':(0,0.3)}, minsnr=5)

    return data, result, fitted_model

def snid_analyze(source):

    ''' Info : Runs SNID analysis on given source and returns classification data
               Saves output files in individual folders
        Input : source
        Returns : classification, rlap score, redshift, redshift error
    '''

    fname = write_ascii_file(source, path=os.getcwd(), auto=False)[0] # Downloads spectrum data in ASCII from Fritz

    if fname == 'No Spectra Found' or fname == 'Resuming...': # Return None if no spectrum on Fritz or if user prompts to continue
        return None, None, None, None

    # Runs SNID shell command, verbose suppresses output in terminal, plot suppresses XWindow, fluxout saves template spectra with 10 highest rlap scores
    bashc = shlex.split(SNID_loc + 'snid verbose=0 plot=0 fluxout=10 ' + os.getcwd() + '/data/' + fname)
    snid = subprocess.Popen(bashc, stdout=subprocess.PIPE, stdin=subprocess.PIPE)

    # SNID takes several seconds to generate an output file, but Python moves on immediately after the subprocess
    # This waits 5 seconds, checking whether or not the output file exists in 0.01 sec intervals
    # If over 5 seconds pass, it assumes no file was generated because SNID did not converge
    count = 0
    while not os.path.exists(fname[:-6]+'_snid.output'):
        time.sleep(0.01)
        count += 1
        if count > 500:
            print('No file exists, there must not be a decent fit.')
            return None, None, None, None

    # The output file also is generated line-by-line so even if the file exists it might not have written up to where we need information
    # This waits until the file is written up to where we seek data
    while True:
        try:
            # This finds the classification with the most templates converged
            info = np.array(open(fname[:-6]+'_snid.output').read().split('\n'))
            start = int(np.argwhere(['type fraction/redshift/age' in i for i in info]))+2
            end = int(np.argwhere(['rlap-ordered template listings' in i for i in info]))-1

            # This finds the templates with the 10 highest rlap scores
            start_f = int(np.argwhere(['rlap-ordered template listings' in i for i in info]))+2
            end_f = start_f + 11

            tab = info[start:end][info[start:end] != '###']
            tab_f = info[start_f:end_f][info[start_f:end_f] != '#--- rlap cutoff']
            if len(tab_f) > 10:
                tab_f = tab_f[:10]

            temp = [s.split()[1] for s in tab_f]
            typ_f = [s.split()[2] for s in tab_f]
            rlap = [s.split()[4] for s in tab_f]
            red = [s.split()[5] for s in tab_f]
            red_err = [s.split()[6] for s in tab_f]

            break
        except IndexError:
            continue

    print(tab_f)

    # Generate outfile directory for transient
    os.mkdir('outfiles/'+fname[:-6])

    # Gather all generated files, again we need to wait for SNID to generate them all so we sleep in 0.01 sec intervals until 12 files are generated
    # It generates the 10 highest rlap templates, plus the output file, and the object spectrum (same in shape to the ASCII file)
    test = os.listdir(os.getcwd())
    fo = [f for f in test if fname[:-6] in f]
    while len(fo) < 12:
        test = os.listdir(os.getcwd())
        fo = [f for f in test if fname[:-6] in f]
        time.sleep(0.01)

    i = 0
    t_spec = ascii.read(fname[:-6]+'_snidflux.dat')
    for item in test:
        if item.startswith(fname[:-6]):
            if item.endswith('.dat') and item != fname[:-6]+'_snidflux.dat':
                fit = ascii.read(item)
                # This will plot all of the template spectra over the Fritz spectra
                plt.figure(i)
                plt.plot(t_spec['wavelength[A]'], t_spec['flux[arbitraty]'])
                plt.plot(fit['col1'], fit['col2'])
                # Title template SNe, classification, rlap score
                plt.title(temp[i] + ', ' + typ_f[i] + ', ' + rlap[i])
                # Also saves figures to same folder
                plt.savefig(os.getcwd()+'/outfiles/'+fname[:-6]+'/'+item[:-4]+'.png')
                i += 1
            # SNID generates files in the working directory, so move the files to the created directory
            shutil.move(item, os.getcwd()+'/outfiles/'+fname[:-6])

    data, result, fitted_model = model_lc(source) # Run light curve fitting on data

    if len(result.parameters) == 5:
        print('Fitted z is ' + str(np.round((result.parameters[0]-z)/z_std, 1)) + ' standard deviations from mean')
        print('Fitted x0 is ' + str(np.round((result.parameters[2]-x0)/x0_std, 1)) + ' standard deviations from mean')
        print('Fitted x1 is ' + str(np.round((result.parameters[3]-x1)/x1_std, 1)) + ' standard deviations from mean')
        print('Fitted c is ' + str(np.round((result.parameters[4]-c)/c_std, 1)) + ' standard deviations from mean')
    elif len(result.parameters) == 4:
        print('Fitted x0 is ' + str(np.round((result.parameters[1]-x0)/x0_std, 1)) + ' standard deviations from mean')
        print('Fitted x1 is ' + str(np.round((result.parameters[2]-x1)/x1_std, 1)) + ' standard deviations from mean')
        print('Fitted c is ' + str(np.round((result.parameters[3]-c)/c_std, 1)) + ' standard deviations from mean')


    plt.figure(i)
    sncosmo.plot_lc(data, model=fitted_model)

    plt.show()

    save = input('Save classification for source? [y/n] ')
    if save == 'y':
        # Templates are ordered by rlap score, so by default we would want the first one
        save2 = input(source + ' will be classified as ' + typ_f[0] + ' with rlap value of ' + rlap[0] + '. Enter nothing to proceed or enter in the number of another template. ')

        if save2 == '':
            return typ_f[0], rlap[0], red[0], red_err[0]
        else:
            try:
                save2 = int(save2)
                return typ_f[save2], rlap[save2], red[save2], red_err[save2]
            except ValueError:
                return None, None, None, None
    else:
        return None, None, None, None

    plt.close('all')
    #print(os.listdir(os.getcwd()))

    return typ[np.argmax(frac)], frac[np.argmax(frac)], red[np.argmax(frac)], red_err[np.argmax(frac)]

def submit_reds(no_reds, f):

    ''' Info : Submits redshift information to Fritz
        Input : list of sources with classifications and no redshift
        Returns : None
    '''

    transients_r, types_r, rlaps_r, reds_r, red_errs_r = run_class(no_reds)

    if len(transients_r) != 0:
        print('ZTFname\t\tRedshift\t\trlaps')
        for tr in np.arange(0,len(transients_r)):
            print(transients_r[tr] + '\t' + reds_r[tr] + ' +/- ' + red_errs_r[tr] + '\t' + rlaps_r[tr])

        up = input('Upload these redshifts? [y/n] ')

        if up == 'y':
            for tr in np.arange(0,len(transients_r)):
                fritz_redshift = submit_fritz_redshift(transients_r[tr], reds_r[tr], red_errs_r[tr])

                if fritz_redshift['status'] == 'success':
                    print(bcolors.OKGREEN + transients_r[tr] + 'classification upload successful.' + bcolors.ENDC)
                else:
                    print(bcolors.FAIL + transients_r[tr] + 'classification upload failed.' + bcolors.ENDC)
                    print(bcolors.FAIL + fritz_redshift['message'] + bcolors.ENDC)

                f['redshift'][np.argwhere(f['Source Name'] == transients_r[tr])] = reds_r[tr]

            f.write('RCF_sources.ascii', format='ascii', overwrite=True, delimiter='\t')

def submit_class(unclassifys, f):

    ''' Info : Submits classification information to Fritz
        Input : list of sources without classifications
        Returns : None
    '''

    transients, types, rlaps, reds, red_errs = run_class(unclassifys) # Runs SNID classification code

    if len(transients) != 0:
        print('ZTFname\t\tClassification\t\tRedshift\t\trlaps')
        for tr in np.arange(0,len(transients)):
            print(transients[tr] + '\t' + types[tr] + '\t\t\t' + reds[tr] + ' +/- ' + red_errs[tr] + '\t' + rlaps[tr])

        up = input('Upload these classifications? [y/n] ')

        if up == 'y':
            for tr in np.arange(0,len(transients)):
                check_r = get_redshift(transients[tr])

                if check_r == 'No redshift found':
                    fritz_redshift = submit_fritz_redshift(transients[tr], reds[tr], red_errs[tr])

                    if fritz_redshift['status'] == 'success':
                        print(bcolors.OKGREEN + transients[tr] + ' redshift upload successful.' + bcolors.ENDC)
                    else:
                        print(bcolors.FAIL + transients[tr] + ' redshift upload failed.' + bcolors.ENDC)
                        print(bcolors.FAIL + fritz_redshift['message'] + bcolors.ENDC)

                    f['redshift'][np.argwhere(f['Source Name'] == transients[tr])] = reds[tr]

                fritz_class = submit_fritz_class(transients[tr], types[tr])

                if fritz_class['status'] == 'success':
                    print(bcolors.OKGREEN + transients[tr] + ' classification upload successful.' + bcolors.ENDC)
                else:
                    print(bcolors.FAIL + transients[tr] + ' classification upload failed.' + bcolors.ENDC)
                    print(bcolors.FAIL + fritz_class['message'] + bcolors.ENDC)

                # Updates RCF source ASCII if user classified any transients
                f['Classification'][np.argwhere(f['Source Name'] == transients[tr])] = types[tr]
                f['Classification Date'][np.argwhere(f['Source Name'] == transients[tr])] = str(datetime.datetime.utcnow().date())

            f.write('RCF_sources.ascii', format='ascii', overwrite=True, delimiter='\t')

def run_class(unclassifys):

    ''' Info : Runs SNID analysis on list of sources
        Input : list of unclassified sources
        Returns : transients, SNID classifications, rlap scores, redshifts, redshift errors
    '''

    print('Clearing directories...')

    # Delete all SNID-generated files to save storage and not confuse the code
    files = os.listdir(os.getcwd())

    test = os.listdir(os.getcwd()+'/data')
    for item in test:
        if item.endswith(".ascii"):
            os.remove(os.path.join(os.getcwd()+'/data', item))

    if 'outfiles' in files:
        test = os.listdir(os.getcwd()+'/outfiles')
        for item in test:
            shutil.rmtree(os.path.join(os.getcwd()+'/outfiles', item))
    else:
        os.mkdir('outfiles')

    test = os.listdir(os.getcwd())
    for item in test:
        if item.startswith("ZTF"):
            os.remove(os.path.join(os.getcwd(), item))

    transients = []
    types = []
    rlaps = []
    reds = []
    red_errs = []

    print('There are ' + str(len(unclassifys)) + ' unclassified transients.')

    for s in np.arange(0,len(unclassifys)):
        print(bcolors.OKCYAN + str(s+1) + '/' + str(len(unclassifys)) + bcolors.ENDC + ': ' + bcolors.OKBLUE + unclassifys[s] + bcolors.ENDC)
        t, f, r, re = snid_analyze(unclassifys[s])

        if t != None:
            if t == 'II':
                t = 'Type II'
            elif t == 'Gal':
                t = 'Galactic Nuclei'
            elif t == 'Ia-csm':
                t = 'Ia-CSM'
            transients.append(unclassifys[s])
            types.append(t)
            rlaps.append(f)
            reds.append(r)
            red_errs.append(re)

    # Saves a csv of classified sources -- can be commented out if necessary
    np.savetxt('SNID_fits.csv', np.rot90(np.fliplr(np.vstack((transients, types, rlaps, reds, red_errs)))), delimiter=',', fmt='%s')

    # Clear the data directory to avoid issues with the submission code
    test = os.listdir(os.getcwd()+'/data')
    for item in test:
        if item.endswith(".ascii"):
            os.remove(os.path.join(os.getcwd()+'/data', item))

    return transients, types, rlaps, reds, red_errs

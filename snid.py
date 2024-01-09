import base64
import copy
import gc
import glob
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import os
import pandas as pd
import shlex
import shutil
import sncosmo
import subprocess
import sys
import time

from astropy.io import ascii
from astropy.table import QTable, Table
from collections import Counter
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import MultipleLocator
from mpl_toolkits.mplot3d import Axes3D
from panoptes_client import Panoptes, Project, SubjectSet, Subject

from func import *
from zooniverse import *

# These numbers come from running model fits on ~500 Type Ia supernovae
z = 0.060574239946858476
z_std = 0.023994157056121096
x1 = -0.14238796934437334
x1_std = 1.4557579021314682
c = 0.08928354223298558
c_std = 0.15670291093588692
x0 = 0.0007648532623426458
x0_std = 0.0004363803462578883

with open('info.info', 'r') as infofile:
    info = infofile.read()
    SNID_loc = info.split('\n')[0].split(':')[1].strip()
    zoo_user = info.split('\n')[5].split(':')[1].strip()
    zoo_pass = info.split('\n')[6].split(':')[1].strip()

# Grab all sources currently uploaded to Zooniverse
zoo = get_all_in_set()

def get_peak_absmag(z, x0):

    ''' Info : Calcultes peak absolute magnitude with SALT2 model parameters
        Input : SALT2-determined redshift and x0
        Returns : Peak absolute magnitude
    '''

    peak_mag = -2.5*np.log10(x0) + 10.635
    cosmo = FlatLambdaCDM(H0=70,Om0=0.3)
    mu = cosmo.distmod(z).value
    absmag = peak_mag -mu

    return absmag

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

    api_dict = api('GET', url, params=data, timeout=10) #defined in order to call values 
    status, response, vers = api_dict.values()

    if format == 'flux':

        flux = []
        fluxerr = []
        band = []
        mjd = []
        zpsys = []
        zp = []

        for d in response:
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

def listComplementElements(list1, list2):


    ''' Info : Finds common elements of two lists
        Input : Two lists
        Returns : Common elements
    '''

    storeResults = []

    for num in list1:
        if num not in list2: # this will essentially iterate your list behind the scenes
            storeResults.append(num)

    return storeResults

def model_lc(source, redshift):

    ''' Info : Fits photometry data to light curve using sncosmo.
        Input : source
        Returns : photometry data, fitted parameters, plottable model
    '''

    data = get_photometry(source)

    red = redshift
    red_err = 'No redshift error found'
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

def plot_box_spec(wave, flux):
    flux_plot = np.repeat(flux, 2)
    wv_plot = wave.copy()
    wv_plot[:-1] += np.diff(wave)/2
    wv_plot = np.append(wave[0]-(wave[1]-wave[0])/2,
                        np.append(np.repeat(wv_plot[0:-1], 2),
                                  wave[-1]+(wave[-1]-wave[-2])/2))

    return wv_plot, flux_plot

def post_lc(source, redshift):

    ''' Info : Posts LC data on Fritz as comment, along with nsigma for c and x1 and peak absolute magnitude. Plot is also attached.
        Input : ZTFname, redshift
        Returns : None
    '''

    data = get_photometry(source)
    comment_infos = get_source_api(source)['comments']

    for i in range (len(get_source_api(source)['comments'])):

        comment_info = comment_infos[i]
        comment = comment_info['text']

        # Check if LC is already posted
        if 'sncosmo light curve fit' in comment:
            if int(comment[int(comment.index('n='))+2:].split(',')[0]) != len(data) or 'gayatri' not in comment: # Check if new photometry has been uploaded

                try:
                    dfit, result, fitted_model = model_lc(source, redshift)
                except RuntimeError:
                    print(bcolors.FAIL + 'sncosmo encountered runtime error. Skipping...' + bcolors.ENDC) # Did not converge on fit
                    return
                except ValueError:
                    print(bcolors.FAIL + 'sncosmo encountered value error. Skipping...' + bcolors.ENDC) # Did not converge on fit
                    return

                x1_nstds = np.round(np.abs((result.parameters[3]-x1)/x1_std), 1)
                c_nstds = np.round(np.abs((result.parameters[4]-c))/c_std, 1)

                sncosmo.plot_lc(dfit, model=fitted_model)

                if np.max(dfit['mjd']) - np.min(dfit['mjd']) < 5: # If <5 nights of photometry, check if user wants to upload
                    plt.show(block=False)

                    res = input('There are only ' + str(np.round(np.max(dfit['mjd']) - np.min(dfit['mjd']), 1)) + ' days worth of photometry data. Do you still want to proceed? [y/n] ')

                    if res != 'y':
                        plt.close()
                        return

                    plt.close()

                plt.savefig('temp.png')

                # If comment exists but new photometry uploaded, edit comment
                resp = edit_comment(source, comment_info['id'], comment_info['author_id'], 'sncosmo light curve fit n='+str(len(data))+', M_peak = '+str(np.round(get_peak_absmag(result.parameters[0], result.parameters[2]),1))+
                    ', x1_nstds = '+str(x1_nstds)+', c_nstds = '+str(c_nstds)+'. LC page: http://gayatri.caltech.edu:88/query/lc/'+source, 'temp.png', source+'_sncosmo_lc.png')

                if resp['status'] == 'success':
                    print(bcolors.OKGREEN + source + ' LC update successful.' + bcolors.ENDC)
                else:
                    print(bcolors.FAIL + source + ' LC update failed.' + bcolors.ENDC)
                    print(bcolors.FAIL + json.dumps(resp, indent=2) + bcolors.ENDC)

                plt.close('all')

                return
            else:
                print(source + ' LC up to date.')
                return

    try:
        dfit, result, fitted_model = model_lc(source, redshift)
    except RuntimeError:
        print(bcolors.FAIL + 'sncosmo encountered runtime error. Skipping...' + bcolors.ENDC)
        return
    except ValueError:
        print(bcolors.FAIL + 'sncosmo encountered value error. Skipping...' + bcolors.ENDC) # Did not converge on fit
        return

    x1_nstds = np.round(np.abs((result.parameters[3]-x1)/x1_std), 1)
    c_nstds = np.round(np.abs((result.parameters[4]-c)/c_std), 1)

    sncosmo.plot_lc(dfit, model=fitted_model)
    plt.savefig('temp.png')

    resp = post_comment(source, 'sncosmo light curve fit n='+str(len(data))+', M_peak = '+str(np.round(get_peak_absmag(result.parameters[0], result.parameters[2]),1))+
        ', x1_nstds = '+str(x1_nstds)+', c_nstds = '+str(c_nstds)+'. LC page: http://gayatri.caltech.edu:88/query/lc/'+source, 'temp.png', source+'_sncosmo_lc.png')

    if resp['status'] == 'success':
        print(bcolors.OKGREEN + source + ' LC upload successful.' + bcolors.ENDC)
    else:
        print(bcolors.FAIL + source + ' LC upload failed.' + bcolors.ENDC)
        print(bcolors.FAIL + json.dumps(resp, indent=2) + bcolors.ENDC)

    plt.close('all')

def plot_best_5(source, output, spectra_name, z_snid, top_5, rlaps, show_redshift=False):
    source_folder = source + spectra_name

    files = np.sort(glob.glob(source_folder+"/*.dat"))
    if(len(files)==0):
        print('return, ' + spectra_name)
        return -1

    rel_files = []
    for i, f in enumerate(files):
        if 'comp' in f and int(f.split('comp')[1].split('_')[0]) in top_5:
            rel_files.append(f)
    rel_files.append(output + '/' + spectra_name + '_snidflux.dat')

    matches, spectra = read_tables(rel_files)

    for spec_num, i in enumerate(matches):
        z = i[0][3]
        snid_type = i[0][2][:-1]

        xi, yi = plot_box_spec(spectra["wavelength"], spectra["flux"])
        xi /= (1+z)
        x, y = i[1]["redshifted_wavelength"] / (1+z), i[1]["flux"]
        specplot(x,y,xi,yi,snid_type,spectra_name,output,i[0][0], z, i[0][4], z_snid, spec_num, rlaps[spec_num], show_redshift=show_redshift)
"""
        if spec_num == 1 and rlaps[spec_num] >= 9.5:
            plt.show(block=False)
        elif spec_num <= 3 and rlaps[spec_num] >= 8.5:
            plt.show(block=False)
"""
def read_tables(files):
    matches_files = files[0:len(files)-1]
    spectra = Table.read(files[-1], format = "ascii", names = ["wavelength", "flux"])
    matches = []
    for i in matches_files:
        input_data = open(i,'r').readlines()[0].split()
        row = [[int(input_data[3][:-1]), input_data[4],input_data[5][1::],float(input_data[-3].split("-")[-1]),float(input_data[-1])]]
        row.append(Table.read(i, format = "ascii", names = ["redshifted_wavelength", "flux"]))
        matches.append(row)
    return matches, spectra

def run_class(unclassifys, unclassified_reds):

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
        t, f, r, re = snid_analyze(unclassifys[s], unclassified_reds[s])

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

def snid_analyze(source, redshift):

    ''' Info : Runs SNID analysis on given source and returns classification data
               Saves output files in individual folders
        Input : source, redshift
        Returns : classification, rlap score, redshift, redshift error
    '''

    if source in np.array(zoo):
        print(source + ' already submitted to Zooniverse within the last 6 months.')
        return None, None, None, None

    fname = write_ascii_file(source, path=os.getcwd(), auto=True)[0] # Downloads spectrum data in ASCII from Fritz

    if fname == None:
        print('Unable to read spectrum.')
        return None, None, None, None

    if fname == 'No Spectra Found' or fname == 'Resuming...': # Return None if no spectrum on Fritz or if user prompts to continue
        return None, None, None, None

    # Runs SNID shell command, verbose suppresses output in terminal, plot suppresses XWindow, fluxout saves template spectra with 100 highest rlap scores
    bashc = shlex.split(SNID_loc + 'snid verbose=0 plot=0 fluxout=100 tempdir=' + SNID_loc + 'templates-2.0/ ' + os.getcwd() + '/data/' + fname)
    snid = subprocess.call(bashc)

    # SNID takes several seconds to generate an output file, but Python moves on immediately after the subprocess
    # This waits 10 seconds, checking whether or not the output file exists in 0.01 sec intervals
    # If over 10 seconds pass, it assumes no file was generated because SNID did not converge
    count = 0
    while not os.path.exists(fname[:-6]+'_snid.output'):
        time.sleep(0.01)
        count += 1
        if count > 1000:
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

    print('Determining best matches...')

    i = 0
    test = os.listdir(os.getcwd())
    t_spec = ascii.read(fname[:-6]+'_snidflux.dat')
    for item in test:
        if item.startswith(fname[:-6]):
            shutil.move(item, os.getcwd()+'/outfiles/'+fname[:-6])

    shutil.copy(os.getcwd() + '/data/' + fname, os.getcwd() + '/outfiles/' + fname[:-6])

    directory = os.getcwd() + "/outfiles/" + fname[:-6] + "/"
    file_list = glob.glob(directory + fname[:-6] + "_snid.output")

    Table_List = []
    sample_length = 1
    count = 0

    for i in file_list:
        templates_list = Table.read(i, format = "ascii", data_start = 28, data_end = 200, names = ["no", "sn", "type", "lap", "rlap", "z", "zerr", "age", "age_flag", "grade"])
        new_i = i.split("/")[-2]
        Table_List.append([new_i, templates_list])
        count += 1
    snidoutput = np.asarray(Table_List)

    ZTable_best = Table(
                    names=("Version", "ZTF_Name",
                           "z_sntemplate", "z_rlap", "z_snid", "z_snid_err", "z_level"
                           , "rank_1", "sntemplate_1", "rlap_1", "c_snid_1", "z_snid_1", "z_snid_err_1", "age_1", "age_flag_1"
                           , "rank_2", "sntemplate_2", "rlap_2", "c_snid_2", "z_snid_2", "z_snid_err_2", "age_2", "age_flag_2"
                           , "rank_3", "sntemplate_3", "rlap_3", "c_snid_3", "z_snid_3", "z_snid_err_3", "age_3", "age_flag_3"
                           , "rank_4", "sntemplate_4", "rlap_4", "c_snid_4", "z_snid_4", "z_snid_err_4", "age_4", "age_flag_4"
                           , "rank_5", "sntemplate_5", "rlap_5", "c_snid_5", "z_snid_5", "z_snid_err_5", "age_5", "age_flag_5"
                    ),
                    meta={"name": "Redscore Results"},
                    dtype=("U64", "U64",
                           "U64", "float64", "float64", "float64", "int32",
                            "int32", "U64", "float64", "U64", "float64", "float64", "float64", "int32",
                            "int32", "U64", "float64", "U64", "float64", "float64", "float64", "int32",
                            "int32", "U64", "float64", "U64", "float64", "float64", "float64", "int32",
                            "int32", "U64", "float64", "U64", "float64", "float64", "float64", "int32",
                            "int32", "U64", "float64", "U64", "float64", "float64", "float64", "int32",
                          )
                    )

    count = 0

    for j in snidoutput:
        j[1] = copy.deepcopy(j[1][0:100])
        nocopys = []
        for linenum in range(len(j[1])):
            name = j[1][linenum]["sn"].split('sn')[-1].split('_b')[0]
            for i in range(linenum, len(j[1])):
                other_name = j[1][i]["sn"].split('sn')[-1].split('_b')[0]
                if(((linenum != i) and ((name in other_name) or (other_name in name)) and (np.abs(j[1][linenum]["age"] - j[1][i]["age"]) < 3))):
                    nocopys.append(i)
        mocopys_unique = listComplementElements(list(range(len(j[1]))),np.unique(nocopys))
        j[1] = j[1][mocopys_unique]

        Top50 = copy.deepcopy(j[1][0:50])
        Top5 = copy.deepcopy(Top50[0:5])
        types = np.unique(Top50["type"])
        top5Types = np.unique(Top5["type"])

        #print(Top5)

        row = []
        row.append(j[0] + ".ascii")
        row.append(j[0].split("_")[0])

        good = j[1][np.where(j[1]["grade"] == "good")]
        if("SLSN" in str(j[1][0]["type"])):
            good = good[np.where(good["z"] <= .5)]
        else:
            good = good[np.where(good["z"] <= .2)]
        if(len(good) != 0):
            row.append(good[0]["sn"])
            row.append(good[0]["rlap"])
            row.append(float(good[0]["z"]))
            row.append(float(good[0]["zerr"]))
            row.append(1)
        else:
            row.append(j[1][0]["sn"])
            row.append(j[1][0]["rlap"])
            row.append(float(j[1][0]["z"]))
            row.append(float(j[1][0]["zerr"]))
            row.append(0)

        if(len(top5Types) >= 3):
            for i in Top5:
                row.append(i["no"])
                row.append(i["sn"])
                row.append(i["rlap"])
                row.append(i["type"])
                row.append(i["z"])
                row.append(i["zerr"])
                row.append(i["age"])
                row.append(i["age_flag"])
        elif(len(top5Types) == 2):
            if(len(types) == 2):
                for i in Top5:
                    row.append(i["no"])
                    row.append(i["sn"])
                    row.append(i["rlap"])
                    row.append(i["type"])
                    row.append(i["z"])
                    row.append(i["zerr"])
                    row.append(i["age"])
                    row.append(i["age_flag"])
            else:
                newType = -1
                for i in range(len(Top50)):
                    line = Top50[i]

                    #print(str(line['type']) in top5Types)
                    #print(len(top5Types))
                    #print(np.array(top5Types))
                    #print(top5Types)

                    if(str(line["type"]) not in top5Types):
                        newType = i
                        break
                for i in range(len(Top5)-1,-1,-1):
                    line = Top5[i]
                    unique = list(Top5["type"]).count(line["type"])
                    if(unique > 1):
                        Top5[i] = Top50[newType]
                        break
                Top5.sort("no")
                for i in Top5:
                    row.append(i["no"])
                    row.append(i["sn"])
                    row.append(i["rlap"])
                    row.append(i["type"])
                    row.append(i["z"])
                    row.append(i["zerr"])
                    row.append(i["age"])
                    row.append(i["age_flag"])
        elif(len(top5Types) == 1):
            if(len(types) == 1):
                for i in Top5:
                    row.append(i["no"])
                    row.append(i["sn"])
                    row.append(i["rlap"])
                    row.append(i["type"])
                    row.append(i["z"])
                    row.append(i["zerr"])
                    row.append(i["age"])
                    row.append(i["age_flag"])
            else:
                for k in range(4,2,-1):
                    top5Types = np.unique(Top5["type"])
                    for i in Top50:
                        if(str(i["type"]) not in top5Types):
                            Top5[k] = i
                            break
                Top5.sort("no")
                for i in Top5:
                    row.append(i["no"])
                    row.append(i["sn"])
                    row.append(i["rlap"])
                    row.append(i["type"])
                    row.append(i["z"])
                    row.append(i["zerr"])
                    row.append(i["age"])
                    row.append(i["age_flag"])
        try:
            ZTable_best.add_row(row)
        except ValueError:
            print(j[0])
        count += 1

    datasource = os.getcwd() + "/outfiles/"
    output = directory

    sample_remaining = ZTable_best

    for i in np.arange(1,6):
        print(str(sample_remaining[0]['rank_' + str(i)]) + '\t' + str(sample_remaining[0]['sntemplate_' + str(i)]) + ' '*(14-len(str(sample_remaining[0]['sntemplate_' + str(i)]))) + '\t' + str(sample_remaining[0]['c_snid_' + str(i)]) + ' '*(10-len(str(sample_remaining[0]['c_snid_' + str(i)]))) + '\t' + str(sample_remaining[0]['rlap_' + str(i)]))

    top_5 = [int(sample_remaining['rank_1']), int(sample_remaining['rank_2']), int(sample_remaining['rank_3']), int(sample_remaining['rank_4']),
        int(sample_remaining['rank_5'])]

    sample_remaining.to_pandas().to_csv(directory + fname[:-6] + '_samp.csv', index = False)

    for i in sample_remaining:
        spectra_name = i["Version"].split(".")[0]
        z_snid = i["z_snid"]
        plot_best_5(datasource,output,spectra_name,z_snid, top_5, [sample_remaining['rlap_1'][0], sample_remaining['rlap_2'][0], sample_remaining['rlap_3'][0],
            sample_remaining['rlap_4'][0], sample_remaining['rlap_5'][0]], show_redshift = False)
        gc.collect()

    #print(sample_remaining)

    try:
        data, result, fitted_model = model_lc(source, redshift) # Run light curve fitting on data

        print('Fitted z is ' + str(np.round((result.parameters[0]-z)/z_std, 1)) + ' standard deviations from mean')
        print('Fitted x0 is ' + str(np.round((result.parameters[2]-x0)/x0_std, 1)) + ' standard deviations from mean')
        print('Fitted x1 is ' + str(np.round((result.parameters[3]-x1)/x1_std, 1)) + ' standard deviations from mean')
        print('Fitted c is ' + str(np.round((result.parameters[4]-c)/c_std, 1)) + ' standard deviations from mean')
    except RuntimeError:
        pass

    save = input('Save classification for source? [y/n] ')
    #save == 'n'
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

        zoo_sub = input('Submit to Zooniverse? [y/n] ')

        if zoo_sub == 'y':

            RedshiftClass = sample_remaining

            indicies = []
            taken = []
            counter = 0
            for i in RedshiftClass:
                if not("Ia" in i["c_snid_1"] and i["rlap_1"] > 10):
                    indicies.append(counter)
                    taken.append(i["ZTF_Name"])
                counter+=1

            manifest = Table(
                names=("!ZTF_Version", "!ZTF_Name", "#Image0", "#Image1", "#Image2", "#Image3", "#Image4", "z_snid", "z_snid_err"
                ),
                meta={"name": "Spectrum Results after SNID"},
                dtype=("U64", "U64", "U64", "U64", "U64", "U64", "U64", "float64", "float64"
                      )
                )

            for i in RedshiftClass:
                row = []
                version = i["Version"].split('.')[0]
                row.append(version)
                row.append(i["ZTF_Name"])

                rank = i["rank_1", "rank_2", "rank_3", "rank_4", "rank_5"]
                list_index = [rank["rank_1"], rank["rank_2"], rank["rank_3"], rank["rank_4"], rank["rank_5"]]

                for j in list_index:
                    row.append("snidfits_emclip_" + version + "_" + str(j)+".png")
                row.append(round(i["z_snid"],3))
                row.append(round(i["z_snid_err"],3))
                manifest.add_row(row)

            manifest.to_pandas().to_csv(directory + "manifest.csv", index = False)

            Panoptes.connect(username=zoo_user, password=zoo_pass)

            project = Project.find(12959)

            subject_set = SubjectSet()

            subject_set.links.project = project
            subject_set.display_name = 'Newly Unclassified'

            subject_set = SubjectSet.find(99282)

            ImageLoc = directory

            subject_set.save()

            count = 1

            new_subjects = []
            f = open(directory + "noFile2.txt", "w")
            for i in manifest:
                try:
                    subject = Subject()

                    subject.links.project = project

                    subject.add_location(ImageLoc + i["#Image0"])
                    subject.add_location(ImageLoc + i["#Image1"])
                    subject.add_location(ImageLoc + i["#Image2"])
                    subject.add_location(ImageLoc + i["#Image3"])
                    subject.add_location(ImageLoc + i["#Image4"])

                    subject.metadata.update({"!ZTF_Version": i["!ZTF_Version"], "!ZTF_Name": i["!ZTF_Name"], "z_snid": i["z_snid"], "z_snid_err": i["z_snid_err"],
                        'rlaps': [sample_remaining['rlap_1'][0], sample_remaining['rlap_2'][0], sample_remaining['rlap_3'][0], sample_remaining['rlap_4'][0], sample_remaining['rlap_5'][0]]})

                    print(subject.metadata)

                    subject.save()
                    new_subjects.append(subject)
                except FileNotFoundError:
                    f.write(i["!ZTF_Version"] + "\n")

            f.close()

            subject_set.add(new_subjects)

            #plt.close('all')
            #print(os.listdir(os.getcwd()))


        return None, None, None, None

def specplot(x, y, xi, yi, snid_type, fname, output, best_num, z_template, z_template_unc, z_snid, spec_num, rlap, show_redshift=False):
    fig, ax = plt.subplots(figsize=(8,4.5))
    ax.plot(xi,yi,color='#32384D',alpha=0.5,
             label='New SN')
    ax.plot(x,y,color='#217CA3',
             label='SNID template', lw=3)
    if show_redshift:
        ax.plot(x[-3],y[-3],color='white',lw=0,
                 label=r'$z_\mathrm{} = $ {:.3f}$\,\pm\,${:.3f}'.format("{SNID}", z_template, z_template_unc))
        ax.text(0.78, 0.955, r'$z_\mathrm{} = ${:.4f}'.format("{SN}", z_snid),
                va='center',
                fontsize=15, transform=plt.gcf().transFigure)
    else:
        ax.text(0.78, 0.955, 'Match #' + str(spec_num+1),
                va='center',
                fontsize=15, transform=plt.gcf().transFigure)

    ax.plot(x[-3],y[-3],color='#217CA3', lw=3)
    ax.set_xlabel(r'Rest Frame Wavelength ($\mathrm{\AA}$)', fontsize=17)
    ax.set_ylabel('Relative Flux', fontsize=17)
    ax.tick_params(which='both',labelsize=15)

    ax.grid(axis='x', color='0.7', ls=':')
    ax.xaxis.set_minor_locator(MultipleLocator(250))
    ax.set_yticklabels([])


    ax.text(0.105, 0.955, 'SNID type: ',
            va='center',
            fontsize=15, transform=plt.gcf().transFigure)
    ax.text(0.245, 0.955, snid_type,
            color='#217CA3', weight='bold', va='center',
            fontsize=23, transform=plt.gcf().transFigure)

    ax.legend(fancybox=True)
    fig.subplots_adjust(left=0.055,right=0.99,top=0.925,bottom=0.145)
    fig.savefig(output + 'snidfits_emclip_' + fname + "_" + str(best_num) + '.png', dpi = 600)
    #print(output + 'snidfits_emclip_' + fname + "_" + str(best_num) + '.png')

    #plt.close(fig)
    #plt.close()

def submit_class(unclassifys, unclassified_reds, f):

    ''' Info : Submits classification information to Fritz
        Input : list of sources without classifications
        Returns : None
    '''

    transients, types, rlaps, reds, red_errs = run_class(unclassifys, unclassified_reds) # Runs SNID classification code

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

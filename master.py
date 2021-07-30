from __future__ import print_function
import numpy as np
import pandas as pd
import os, glob2
import requests, json
from lxml import html
from astropy.time import Time
import pickle
from astropy.table import Table
from astropy.io import ascii
import os
import requests
import json
import datetime
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from astropy import constants as const
import sys, getopt, argparse
import re
from time import sleep
from astropy.io import fits
from subprocess import call
from lxml import html
import webbrowser as wb
from urllib.error import HTTPError
import xlsxwriter
from tqdm import tqdm
from func import *
from snid import *

# Create data directory if one does not exist
test = os.listdir(os.getcwd())
if 'data' not in test:
    os.mkdir('data')

# Opens file with last date and time the script was run
if 'last_date.txt' in test:
    with open('last_date.txt', 'r+') as file:
        last_date = datetime.datetime.strptime(file.read(), '%Y-%m-%d')
        file.seek(0)
        file.truncate(0)
        file.write(str(datetime.date.today()))
else:
    with open('last_date.txt', 'w') as file:
        file.write(str(datetime.date.today() - datetime.timedelta(days=1)))

print('Welcome to the master script!')
if 'RCF_sources.ascii' not in test:
    since = input('Enter earliest date to download sources from (YYYY-MM-DD) or enter nothing to set it to 6 months ago: ')
    get_source_file('RCF_sources', since)
else:
    dl = input('Download new list of RCF sources? ([y]/n) ')

    if dl == 'y':
        since = input('Enter earliest date to download sources from (YYYY-MM-DD) or enter nothing to set it to 6 months ago: ')
        get_source_file('RCF_sources', since)

f = ascii.read("RCF_sources.ascii") #ascii file containing the names of sources and their saved dates

entered = input('Enter in the earliest date you want to check classifications or saves (YYYY-MM-DD), \'y\' for yesterday at midnight, or enter nothing to select last time program was run: ')

if entered == '':
    startd = last_date
elif entered == 'Y' or entered == 'y':
    startd = datetime.datetime.combine((datetime.date.today() - datetime.timedelta(days=1)), datetime.datetime.min.time())
else:
    startd = datetime.datetime.strptime(entered, '%Y-%m-%d')

sources, tns_names, savedates, classifys, class_dates, reds, unclassifys = read_ascii(f, startd) # Parses data from ASCII file according to previous input

print(bcolors.OKGREEN + 'Checking for missing redshifts...' + bcolors.ENDC)

if len(reds[reds=='No redshift found']) == 0:
    print('All classified sources have listed redshifts.')
else:
    red_resp = input(str(len(reds[reds=='No redshift found'])) + ' classified sources have no redshift listed. Use SNID to determine? [y/n] ')

    if red_resp == 'y':
        transients_r, types_r, rlaps_r, reds_r, red_errs_r = run_class(sources[reds=='No redshift found'])

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

                    f['redshift'][np.argwhere(f['Source Name'] == transients_r[tr])] = reds_r[tr]

                f.write('RCF_sources.ascii', format='ascii', overwrite=True, delimiter='\t')

                sources, tns_names, savedates, classifys, class_dates, reds, unclassifys = read_ascii(f, startd) # Reload RCF source file with newly classified transients

print(bcolors.OKGREEN + 'Checking for unclassified transients...' + bcolors.ENDC)

if len(unclassifys) == 0:
    print('There are no unclassified sources since the date entered.')
else:
    c = input(str(len(unclassifys)) + ' unclassified transients on Fritz have been saved. Would you like to classify them? [y/n] ')

    if c == 'y':

        transients, types, rlaps, reds, red_errs = run_class(unclassifys) # Runs SNID classification code

        if len(transients) != 0:
            print('ZTFname\t\tClassification\t\tRedshift\t\trlaps')
            for tr in np.arange(0,len(transients)):
                print(transients[tr] + '\t' + types[tr] + '\t\t\t' + reds[tr] + ' +/- ' + red_errs[tr] + '\t' + rlaps[tr])

            up = input('Upload these classifications? [y/n] ')

            if up == 'y':
                for tr in np.arange(0,len(transients)):
                    check_r = get_redshift(transients[tr])

                    if check_r != 'No redshift found':
                        fritz_redshift = submit_fritz_redshift(transients[tr], reds[tr], red_errs[tr])

                        if fritz_redshift['status'] == 'success':
                            print(bcolors.OKGREEN + transients[tr] + ' redshift upload successful.' + bcolors.ENDC)
                        else:
                            print(bcolors.FAIL + transients[tr] + ' redshift upload failed.' + bcolors.ENDC)

                        f['redshift'][np.argwhere(f['Source Name'] == transients[tr])] = reds[tr]

                    fritz_class = submit_fritz_class(transients[tr], types[tr])

                    if fritz_class['status'] == 'success':
                        print(bcolors.OKGREEN + transients[tr] + ' classification upload successful.' + bcolors.ENDC)
                    else:
                        print(bcolors.FAIL + transients[tr] + ' classification upload failed.' + bcolors.ENDC)

                    # Updates RCF source ASCII if user classified any transients
                    f['Classification'][np.argwhere(f['Source Name'] == transients[tr])] = types[tr]
                    f['Classification Date'][np.argwhere(f['Source Name'] == transients[tr])] = str(datetime.date.today())

                f.write('RCF_sources.ascii', format='ascii', overwrite=True, delimiter='\t')

                sources, tns_names, savedates, classifys, class_dates, reds, unclassifys = read_ascii(f, startd) # Reload RCF source file with newly classified transients

print(bcolors.OKGREEN + 'Beginning TNS submissions...' + bcolors.ENDC)

print('There are ' + str(len(sources)) + ' objects saved or classified later than ' + str(startd) + ' with classifications.')

class_submission(sources, tns_names, classifys, class_dates) # Runs TNS submission script

print('Submission complete. Goodbye!')

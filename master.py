from __future__ import print_function

import astropy.units as u
import datetime
import json
import os, glob2
import numpy as np
import pandas as pd
import pickle
import pytz
import re
import requests, json, simplejson
import sys, getopt, argparse
import webbrowser as wb
import xlsxwriter

from astropy import constants as const
from astropy.cosmology import FlatLambdaCDM
from astropy.io import ascii, fits
from astropy.table import Table
from astropy.time import Time
from lxml import html
from subprocess import call
from time import sleep
from tqdm import tqdm
from urllib.error import HTTPError

if 'info.info' not in os.listdir(os.getcwd()): # Retrieves API key info and location of SNID
    print('No info file in directory! "info.info" has been generated, enter in the location of SNID, and API information.')
    with open('info.info', 'w') as f:
        f.write('SNID loc: \nFritz API key: \nTNS API Key: \nTNS Bot ID: \nZooniverse username: \nZoonviverse Password: ')
    exit()

from func import *
from snid import *
from hosts import *
from zooniverse import *

# Create data directory if one does not exist
test = os.listdir(os.getcwd())
if 'data' not in test:
    os.mkdir('data')

print('Welcome to the master script!')

entered = input('Enter in the earliest date you want to check classifications or saves (YYYY-MM-DD) or \'y\' for yesterday at midnight: ')

if entered == 'Y' or entered == 'y':
    startd = datetime.datetime.combine((datetime.datetime.utcnow().date() - datetime.timedelta(days=1)), datetime.datetime.min.time()).replace(tzinfo=datetime.timezone.utc)
else:
    startd = datetime.datetime.strptime(entered, '%Y-%m-%d').replace(tzinfo=datetime.timezone.utc)

if input('Check for new Zooniverse classifications? [y/n] ') == 'y':
    pull_class(startd)

print(bcolors.OKGREEN + 'Continuing...' + bcolors.ENDC)

if 'RCF_sources.ascii' not in test:
    since = input('Enter earliest date to download sources from (YYYY-MM-DD) or enter nothing to set it to 6 months ago: ')
    get_source_file('RCF_sources', since)
else:
    dl = input('Download new list of RCF sources? ([y]/n) ')

    if dl == 'y':
        since = input('Enter earliest date to download sources from (YYYY-MM-DD) or enter nothing to set it to 6 months ago: ')
        get_source_file('RCF_sources', since)

f = ascii.read("RCF_sources.ascii") #ascii file containing the names of sources and their saved dates

sources, tns_names, savedates, classifys, class_dates, reds, unclassifys = read_ascii(f, startd) # Parses data from ASCII file according to previous input

if 'zooniverse.ascii' not in os.listdir(os.getcwd()):
    with open('zooniverse.ascii', 'w') as file:
        file.write('ztfname\tdate\n')

converters = {'ztfname': [ascii.convert_numpy(np.str)], 'date': [ascii.convert_numpy(np.str)]}
zoo = ascii.read('zooniverse.ascii', converters=converters)

remove_indices = []
for z in range(len(zoo)):
    if datetime.datetime.strptime(zoo['date'][z], '%Y-%m-%d').replace(tzinfo=datetime.timezone.utc).date() < datetime.datetime.utcnow().date() - datetime.timedelta(days=180):
        remove_indices.append(z)

zoo.remove_rows(remove_indices)
ascii.write(zoo, 'zooniverse.ascii', overwrite=True)

option = ''
while option != 0: # Select options
    print('1: Check for missing redshifts\n2: Classify unclassified sources\n3: Check and upload light curve data\n4: Associate hosts with saved sources\n5: Submit Fritz classifications to TNS')
    option = int(input('Enter in what you want to do, 0 to exit, or "all" to do all: '))

    if option == 1 or option == 'all':

        print(bcolors.OKGREEN + 'Checking for missing redshifts...' + bcolors.ENDC)

        if len(reds[reds=='No redshift found']) == 0:
            print('All classified sources have listed redshifts.')
        else:
            print(str(len(reds[reds=='No redshift found'])) + ' classified sources have no redshift listed.')

            submit_reds(sources[reds=='No redshift found'], f)

            sources, tns_names, savedates, classifys, class_dates, reds, unclassifys = read_ascii(f, startd) # Reload RCF source file with transients with newly determined redshifts

    if option == 2 or option == 'all':

        print(bcolors.OKGREEN + 'Checking for unclassified transients...' + bcolors.ENDC)

        if len(unclassifys) == 0:
            print('There are no unclassified sources since the date entered.')
        else:
            print(str(len(unclassifys)) + ' unclassified transients on Fritz have been saved. Would you like to classify them? [y/n] ')

            submit_class(unclassifys, f)

            sources, tns_names, savedates, classifys, class_dates, reds, unclassifys = read_ascii(f, startd) # Reload RCF source file with newly classified transients

    if option == 3 or option == 'all':

        print(bcolors.OKGREEN + 'Checking for updated photometry data...' + bcolors.ENDC)

        phot_sources = np.concatenate(([sources[s] for s in np.arange(0,len(sources)) if 'Ia' in classifys[s]], unclassifys))

        for p in np.arange(0,len(phot_sources)):
            print(bcolors.OKCYAN + str(p+1) + '/' + str(len(phot_sources)) + bcolors.ENDC + ': ' + bcolors.OKBLUE + phot_sources[p] + bcolors.ENDC)
            post_lc(phot_sources[p])

    if option == 4 or option == 'all':

        print(bcolors.OKGREEN + 'Checking for hosts...' + bcolors.ENDC)

        saved_sources = np.append(sources, unclassifys)

        for i in np.arange(0,len(saved_sources)):
            print(bcolors.OKCYAN + str(i+1) + '/' + str(len(saved_sources)) + bcolors.ENDC + ': ' + bcolors.OKBLUE + saved_sources[i] + bcolors.ENDC)
            post_host(saved_sources[i])

    if option == 5 or option == 'all':
        print(bcolors.OKGREEN + 'Beginning TNS submissions...' + bcolors.ENDC)

        print('There are ' + str(len(sources)) + ' objects saved or classified later than ' + str(startd) + ' with classifications.')

        class_submission(sources, tns_names, classifys, class_dates) # Runs TNS submission script

print('Submission complete. Goodbye!')

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

if 'info.info' not in os.listdir(os.getcwd()): # Retrieves API key info and location of SNID
    print('No info file in directory! "info.info" has been generated, enter in the location of SNID, and API information.')
    with open('info.info', 'w') as f:
        f.write('SNID loc: \nFritz API key: \nTNS API Key: \nTNS Bot ID: ')
    exit()

from func import *
from snid import *

# Create data directory if one does not exist
test = os.listdir(os.getcwd())
if 'data' not in test:
    os.mkdir('data')

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

entered = input('Enter in the earliest date you want to check classifications or saves (YYYY-MM-DD) or \'y\' for yesterday at midnight: ')

if entered == 'Y' or entered == 'y':
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
        submit_reds(sources[reds=='No redshift found'])

        sources, tns_names, savedates, classifys, class_dates, reds, unclassifys = read_ascii(f, startd) # Reload RCF source file with transients with newly determined redshifts

print(bcolors.OKGREEN + 'Checking for unclassified transients...' + bcolors.ENDC)

if len(unclassifys) == 0:
    print('There are no unclassified sources since the date entered.')
else:
    c = input(str(len(unclassifys)) + ' unclassified transients on Fritz have been saved. Would you like to classify them? [y/n] ')

    if c == 'y':
        submit_class(sources)

        sources, tns_names, savedates, classifys, class_dates, reds, unclassifys = read_ascii(f, startd) # Reload RCF source file with newly classified transients

print(bcolors.OKGREEN + 'Beginning TNS submissions...' + bcolors.ENDC)

print('There are ' + str(len(sources)) + ' objects saved or classified later than ' + str(startd) + ' with classifications.')

class_submission(sources, tns_names, classifys, class_dates) # Runs TNS submission script

print('Submission complete. Goodbye!')

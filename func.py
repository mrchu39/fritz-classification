from __future__ import print_function

import astropy.units as u
import base64
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

with open('info.info', 'r') as f:
    info = f.read().split('\n')
    ft = info[1].split(':')[1].strip()
    tns_apikey = info[2].split(':')[1].strip()
    tns_botid = info[3].split(':')[1].strip()

global TOKEN, BASEURL
GETTOKEN = ft      # Fritz API Key, retrieves from info file
BASEURL = 'https://fritz.science/'                     # Fritz base url

API_KEY = tns_apikey     # TNS API Key from info file
YOUR_BOT_ID = tns_botid
YOUR_BOT_NAME="ZTF_Bot1"

# TNS URLs for real uploads
TNS_BASE_URL = "https://www.wis-tns.org/api/"
upload_url = "https://www.wis-tns.org/api/file-upload"
report_url = "https://www.wis-tns.org/api/bulk-report"
reply_url = "https://www.wis-tns.org/api/bulk-report-reply"

# SANDBOX URLs for TNS upload trials
SAND_TNS_BASE_URL = "https://sandbox-tns.org/api/"
SAND_upload_url = "https://sandbox-tns.org/api/"
SAND_report_url = "https://sandbox-tns.org/api/bulk-report"
SAND_reply_url = "https://sandbox-tns.org/api/bulk-report-reply"

class bcolors:

    ''' Info : Colors for console output.
        Attributes: UTF-8 codes for colors as strings
    '''

    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

class TNSClassificationReport:

    ''' Info : TNS Classification Report object, contains information to be submitted to TNS
        Attributes: self-explanatory
    '''

    def __init__(self):
        self.name = ''
        self.fitsName = ''
        self.asciiName = ''
        self.classifierName = ''
        self.classificationID = ''
        self.redshift = ''
        self.classificationComments = ''
        self.obsDate = ''
        self.instrumentID = ''
        self.expTime = ''
        self.observers = ''
        self.reducers = ''
        self.specTypeID = ''
        self.spectrumComments = ''
        self.groupID = ''
        self.spec_proprietary_period_value = ''
        self.spec_proprietary_period_units = ''


    def fill(self):

        ''' Info : Generate TNS submission dict
            Input : self
            Returns : Dictionary to submit to TNS
        '''

        spectrumdict = {
            'obsdate': self.obsDate,
            'instrumentid': self.instrumentID,
            'exptime': self.expTime,
            'observer': self.observers,
            'reducer': self.reducers,
            'spectypeid': self.specTypeID,
            'ascii_file': self.asciiName,
            'fits_file': self.fitsName,
            'remarks': self.spectrumComments,
            'spec_proprietary_period' : self.spec_proprietary_period_value}

        classificationdict =  {
            'classification_report': {
                '0': {
                    'name': self.name,
                    'classifier': self.classifierName,
                    'objtypeid': self.classificationID,
                    'redshift': self.redshift,
                    'groupid': self.groupID,
                    'remarks': self.classificationComments,
                    'spectra': {
                        'spectra-group': {
                            '0': spectrumdict
                        }
                    }
                }
            }
        }

        return classificationdict

    def as_json(self):

        ''' Info : Returns the dictionary as a formatted string
            Input : self
            Returns : classification dictionary as formatted string
        '''

        return json.dumps(self.fill())

def api(method, endpoint, data=None, params=None, timeout=10):
    ''' Info : Basic API query, takes input the method (eg. GET, POST, etc.), the endpoint (i.e. API url)
               and additional data for filtering
        Returns : response in json format
        CAUTION! : If the query doesn't go through, try putting the 'data' input in 'data' or 'params'
                    argument in requests.request call
    '''

    headers = {'Authorization': f'token {GETTOKEN}'}

    while True:
        try:
            response = requests.request(method, endpoint, json=data, headers=headers, params=params, timeout=timeout)
            re_dict = response.json()

            if '429 Too Many Requests' not in response.text:
                return re_dict
            else:
                #print('429')
                continue
        except requests.exceptions.Timeout:
            print('Timeout Exception, restarting...')
            continue
        except (json.decoder.JSONDecodeError, simplejson.errors.JSONDecodeError):
            #print('JSON Decode Error, restarting...')
            continue

def APO(specid):

    ''' Info : Retrieves spectrum's observing run info from Fritz API
        Input : specid
        Returns : observation date, exposure time, observers, reducers
    '''

    a = get_spectrum_api(specid)
    inst = (a['data']['instrument_name'])

    obsdate = a['data']['original_file_string'].split('#')[6]
    a,b = obsdate.split(' ', 1)
    c,obsdate = b.split(' ', 1)

    OBSDATE = obsdate

    a = get_spectrum_api(specid)

    exptime = a['data']['original_file_string'].split('#')[9]
    a,b = exptime.split(' ', 1)
    c,EXPTIME = b.split(' ', 1)

    a = get_spectrum_api(specid)

    observers = (a['data']['original_file_string'].split('#')[10])
    a,b = observers.split(' ', 1)
    c,OBSERVERS = b.split(' ', 1)

    a = get_spectrum_api(specid)

    reducers = a['data']['original_file_string'].split('#')[11]
    a,b = reducers.split(' ', 1)
    c,d = b.split(' ', 1)
    REDUCERS,e = d.split('\n', 1)

    return OBSDATE.split(' \n')[0], EXPTIME.split(' \n')[0], OBSERVERS.split(' \n')[0], REDUCERS

def check_TNS_class(ztfname):

    ''' Info : Checks TNS page for other classification reports
        Input : ZTFname
        Returns : Group that reported classification
    '''

    if len(get_IAUname(ztfname)) == 0:
        return
    tns_name = get_IAUname(ztfname)[0]["name"][3:]
    data = {'api_key' : API_KEY}
    headers={'User-Agent':'tns_marker{"tns_id":'+str(YOUR_BOT_ID)+', "type":"bot", "name":"'+YOUR_BOT_NAME+'"}'}
    response = requests.get('https://www.wis-tns.org/object/'+tns_name, headers=headers, data=data)
    class_data = response.text.split('Classification Reports', 2)[1]

    if 'no-data' in class_data.split('class="clear"')[0]:
        return

    class_info = class_data.split('class="odd"')[1]

    print(ztfname + ' classified as ' + class_info.split('"cell-type">')[1].split('<')[0] + ' by ' + class_info.split('"cell-user_name">')[1].split('<')[0] + ' at ' +
          class_info.split('"cell-time_received">')[1].split('<')[0] + '.')

    return class_info.split('"cell-source_group_name">')[1].split('<')[0]

def class_submission(sources, tns_names, classifys, class_dates):

    ''' Info : Takes source list and submits those with classifications that have not been submitted prior by ZTF
        Input : source list, TNS names of sources, classifications, classifiction dates
        Returns : None
    '''

    sc = -1
    for source in sources:
        sc += 1

        print(bcolors.OKCYAN + str(sc+1) + '/' + str(len(sources)) + bcolors.ENDC + ': ' + bcolors.OKBLUE + source + bcolors.ENDC)

        flag = 0
        ztfname = source

        for i in range (len(get_source_api(source)['comments'])):

            comment = get_source_api(source)['comments'][i]['text']

            if comment == 'Uploaded to TNS':
                print(ztfname + ' already uploaded to TNS.')
                flag = 1
                continue

            if comment == 'Classification from TNS':
                print(source + ' classified from TNS.')
                flag = 1
                continue

        if flag == 0:

            if tns_names[sc] == 'Not reported to TNS':
                print(ztfname + ' not reported to TNS yet.')
                continue

            if check_TNS_class(ztfname) == 'ZTF':
                if input(ztfname + ' already uploaded to TNS by ZTF, submit another classification? [y/n] ') != 'y':
                    continue

            class_date = class_dates[sc]
            classify = classifys[sc]

            #print (classify, class_date)

            #print(savedates)

            spectrum_info = write_ascii_file(ztfname) #returns "spectrum_name"
            spectrum_name = spectrum_info[0]

            if spectrum_name != 'No Spectra Found' and spectrum_name != 'Resuming...':

                if flag == 0:

                    print(ztfname + '\t' + tns_names[sc] + '\t' + classify + '\t' + class_date)

                    path = os.getcwd()

                    specfile = (path+'/data/'+spectrum_name)

                    files = specfile

                    specid = spectrum_info[1]

                    a = get_spectrum_api(specid)

                    inst = (a['data']['instrument_name'])


                    if inst == 'SEDM':

                        classifiers = 'M. Chu, A. Dahiwale, C. Fremling(Caltech) on behalf of the Zwicky Transient Facility (ZTF)'### Change accordingly
                        source_group = 48 ### Require source group id from drop down list, 0 is for None
                        spectypes = np.array(['object','host','sky','arcs','synthetic'])

                        #proprietary_period = int(input("Proprietary period in years:", x)
                        proprietary_period = '0'
                        proprietary_units = "years"
                        spec_comments =''
                        classification_comments = ''
                        spectype='object'
                        spectype_id = ['object', 'host', 'sky', 'arcs', 'synthetic'].index(spectype) + 1

                        header = (a['data']['altdata'])
                        obsdate = str((header['UTC']).split('T')[0])+' '+str((header['UTC']).split('T')[1])

                        classificationReport = TNSClassificationReport()
                        classificationReport.name = get_IAUname(ztfname)[0]['name'][3:]
                        classificationReport.fitsName = ''
                        classificationReport.asciiName = spectrum_name
                        classificationReport.classifierName = classifiers
                        classificationReport.classificationID = get_TNS_classification_ID(classify)
                        classificationReport.redshift = get_redshift(ztfname)
                        classificationReport.classificationComments = classification_comments
                        classificationReport.obsDate = obsdate
                        classificationReport.instrumentID = get_TNS_instrument_ID(inst)
                        classificationReport.expTime = (header['EXPTIME'])
                        classificationReport.observers = 'SEDmRobot'
                        classificationReport.reducers = (header['REDUCER'])
                        classificationReport.specTypeID = spectype_id
                        classificationReport.spectrumComments = spec_comments
                        classificationReport.groupID = source_group
                        classificationReport.spec_proprietary_period_value = proprietary_period
                        classificationReport.spec_proprietary_period_units = proprietary_units

                        pprint(classificationReport.fill(), tab='  ')
                        proceed = input("\nProceed with classification and upload? ([y]/n) : ")
                        if proceed == 'y' and not proceed.strip() == '':

                            # ASCII FILE UPLOAD
                            print ("\n")
                            response = upload_to_TNS(files)
                            print (response)

                            if not response:
                                print("File upload didn't work")
                                print(response)
                                #return False

                            print(response['id_code'], response['id_message'],
                                  "\nSuccessfully uploaded ascii spectrum")
                            #classificationReport.asciiName = response['data'][-1]

                            report_id = tns_classify(classificationReport)
                            if tns_feedback(report_id) == True:
                                post_comment(ztfname, 'Uploaded to TNS')


                    if inst == 'SPRAT':

                        classifiers = 'D. A. Perley(LJMU), M. Chu, A. Dahiwale, C. Fremling (Caltech) on behalf of the Zwicky Transient Facility (ZTF)'### Change accordingly
                        source_group = 48 ### Require source group id from drop down list, 0 is for None
                        spectypes = np.array(['object','host','sky','arcs','synthetic'])

                        #proprietary_period = int(input("Proprietary period in years:", x)
                        proprietary_period = '0'
                        proprietary_units = "years"
                        spec_comments =''
                        classification_comments = ''
                        spectype='object'
                        spectype_id = ['object', 'host', 'sky', 'arcs', 'synthetic'].index(spectype) + 1

                        header = (a['data']['altdata'])
                        obsdate = str(header['OBSDATE'].split('T')[0])+' '+str(header['OBSDATE'].split('T')[1])

                        classificationReport = TNSClassificationReport()
                        classificationReport.name = get_IAUname(ztfname)[0]['name'][3:]
                        classificationReport.fitsName = ''
                        classificationReport.asciiName = spectrum_name
                        classificationReport.classifierName = classifiers
                        classificationReport.classificationID = get_TNS_classification_ID(classify)
                        classificationReport.redshift = get_redshift(ztfname)
                        classificationReport.classificationComments = classification_comments
                        classificationReport.obsDate = obsdate
                        classificationReport.instrumentID = get_TNS_instrument_ID(inst)
                        classificationReport.expTime = (header['EXPTIME'])
                        classificationReport.observers = 'LTRobot'
                        classificationReport.reducers = 'D. Perley'
                        classificationReport.specTypeID = spectype_id
                        classificationReport.spectrumComments = spec_comments
                        classificationReport.groupID = source_group
                        classificationReport.spec_proprietary_period_value = proprietary_period
                        classificationReport.spec_proprietary_period_units = proprietary_units

                        pprint(classificationReport.fill(), tab='  ')
                        proceed = input("\nProceed with classification and upload? ([y]/n) : ")
                        if proceed == 'y' and not proceed.strip() == '':

                            # ASCII FILE UPLOAD
                            print ("\n")
                            response = upload_to_TNS(files)
                            print (response)

                            if not response:
                                print("File upload didn't work")
                                print(response)
                                #return False

                            print(response['id_code'], response['id_message'],
                                  "\nSuccessfully uploaded ascii spectrum")
                            #classificationReport.asciiName = response['data'][-1]

                            report_id = tns_classify(classificationReport)
                            if tns_feedback(report_id) == True:
                                post_comment(ztfname, 'Uploaded to TNS')



                    if inst == 'ALFOSC':

                        classifiers = 'M. Chu, A. Dahiwale, C. Fremling(Caltech) on behalf of the Zwicky Transient Facility (ZTF)'### Change accordingly
                        source_group = 48 ### Require source group id from drop down list, 0 is for None
                        spectypes = np.array(['object','host','sky','arcs','synthetic'])

                        #proprietary_period = int(input("Proprietary period in years:", x)
                        proprietary_period = '0'
                        proprietary_units = "years"
                        spec_comments =''
                        classification_comments = ''
                        spectype='object'
                        spectype_id = ['object', 'host', 'sky', 'arcs', 'synthetic'].index(spectype) + 1

                        header = (a['data']['altdata'])
                        obsdate = str(a['data']['observed_at'].split('T')[0])+' '+str(a['data']['observed_at'].split('T')[1])

                        classificationReport = TNSClassificationReport()
                        classificationReport.name = get_IAUname(ztfname)[0]['name'][3:]
                        classificationReport.fitsName = ''
                        classificationReport.asciiName = spectrum_name
                        classificationReport.classifierName = classifiers
                        classificationReport.classificationID = get_TNS_classification_ID(classify)
                        classificationReport.redshift = get_redshift(ztfname)
                        classificationReport.classificationComments = classification_comments
                        classificationReport.obsDate = obsdate
                        classificationReport.instrumentID = get_TNS_instrument_ID(inst)
                        classificationReport.expTime = (header['EXPTIME'])
                        classificationReport.observers = (str(a['data']['observers'][0]['first_name'])+' '+str(a['data']['observers'][0]['last_name']))
                        classificationReport.reducers = (str(a['data']['reducers'][0]['first_name'])+' '+str(a['data']['reducers'][0]['last_name']))
                        classificationReport.specTypeID = spectype_id
                        classificationReport.spectrumComments = spec_comments
                        classificationReport.groupID = source_group
                        classificationReport.spec_proprietary_period_value = proprietary_period
                        classificationReport.spec_proprietary_period_units = proprietary_units

                        pprint(classificationReport.fill(), tab='  ')
                        proceed = input("\nProceed with classification and upload? ([y]/n) : ")
                        if proceed == 'y' and not proceed.strip() == '':

                            # ASCII FILE UPLOAD
                            print ("\n")
                            response = upload_to_TNS(files)
                            print (response)

                            if not response:
                                print("File upload didn't work")
                                print(response)
                                #return False

                            print(response['id_code'], response['id_message'],
                                  "\nSuccessfully uploaded ascii spectrum")
                            #classificationReport.asciiName = response['data'][-1]

                            report_id = tns_classify(classificationReport)
                            if tns_feedback(report_id) == True:
                                post_comment(ztfname, 'Uploaded to TNS')



                    if inst == 'DBSP':

                        classifiers = 'M. Chu, A. Dahiwale, C. Fremling(Caltech) on behalf of the Zwicky Transient Facility (ZTF)'### Change accordingly
                        source_group = 48 ### Require source group id from drop down list, 0 is for None
                        spectypes = np.array(['object','host','sky','arcs','synthetic'])

                        #proprietary_period = int(input("Proprietary period in years:", x)
                        proprietary_period = '0'
                        proprietary_units = "years"
                        spec_comments =''
                        classification_comments = ''
                        spectype='object'
                        spectype_id = ['object', 'host', 'sky', 'arcs', 'synthetic'].index(spectype) + 1

                        OBSDATE = str(a['data']['observed_at'].split('T')[0])+' '+str(a['data']['observed_at'].split('T')[1])

                        classificationReport = TNSClassificationReport()
                        classificationReport.name = get_IAUname(ztfname)[0]['name'][3:]
                        classificationReport.fitsName = ''
                        classificationReport.asciiName = spectrum_name
                        classificationReport.classifierName = classifiers
                        classificationReport.classificationID = get_TNS_classification_ID(classify)
                        classificationReport.redshift = get_redshift(ztfname)
                        classificationReport.classificationComments = classification_comments
                        classificationReport.obsDate = OBSDATE
                        classificationReport.instrumentID = get_TNS_instrument_ID(inst)
                        #classificationReport.expTime = '900'
                        classificationReport.observers = (str(a['data']['observers'][0]['first_name'])+' '+str(a['data']['observers'][0]['last_name']))
                        classificationReport.reducers = (str(a['data']['reducers'][0]['first_name'])+' '+str(a['data']['reducers'][0]['last_name']))
                        classificationReport.specTypeID = spectype_id
                        classificationReport.spectrumComments = spec_comments
                        classificationReport.groupID = source_group
                        classificationReport.spec_proprietary_period_value = proprietary_period
                        classificationReport.spec_proprietary_period_units = proprietary_units

                        pprint(classificationReport.fill(), tab='  ')
                        proceed = input("\nProceed with classification and upload? ([y]/n) : ")
                        if proceed == 'y' and not proceed.strip() == '':

                            # ASCII FILE UPLOAD
                            print ("\n")
                            response = upload_to_TNS(files)
                            print (response)

                            if not response:
                                print("File upload didn't work")
                                print(response)
                                #return False

                            print(response['id_code'], response['id_message'],
                                  "\nSuccessfully uploaded ascii spectrum")
                            #classificationReport.asciiName = response['data'][-1]

                            report_id = tns_classify(classificationReport)
                            if tns_feedback(report_id) == True:
                                post_comment(ztfname, 'Uploaded to TNS')


                    if inst == 'LRIS':

                        classifiers = 'M. Chu, A. Dahiwale, C. Fremling(Caltech) on behalf of the Zwicky Transient Facility (ZTF)'### Change accordingly
                        source_group = 48 ### Require source group id from drop down list, 0 is for None
                        spectypes = np.array(['object','host','sky','arcs','synthetic'])

                        header = (a['data']['altdata'])

                        #proprietary_period = int(input("Proprietary period in years:", x)
                        proprietary_period = '0'
                        proprietary_units = "years"
                        spec_comments =''
                        classification_comments = ''
                        spectype='object'
                        spectype_id = ['object', 'host', 'sky', 'arcs', 'synthetic'].index(spectype) + 1

                        OBSDATE = str(a['data']['observed_at'].split('T')[0])+' '+str(a['data']['observed_at'].split('T')[1])

                        classificationReport = TNSClassificationReport()
                        classificationReport.name = get_IAUname(ztfname)[0]['name'][3:]
                        classificationReport.fitsName = ''
                        classificationReport.asciiName = spectrum_name
                        classificationReport.classifierName = classifiers
                        classificationReport.classificationID = get_TNS_classification_ID(classify)
                        classificationReport.redshift = get_redshift(ztfname)
                        classificationReport.classificationComments = classification_comments
                        classificationReport.obsDate = OBSDATE
                        classificationReport.instrumentID = get_TNS_instrument_ID(inst)
                        classificationReport.expTime = '300'
                        classificationReport.observers = (str(a['data']['observers'][0]['first_name'])+' '+str(a['data']['observers'][0]['last_name']))
                        classificationReport.reducers = (str(a['data']['reducers'][0]['first_name'])+' '+str(a['data']['reducers'][0]['last_name']))
                        classificationReport.specTypeID = spectype_id
                        classificationReport.spectrumComments = spec_comments
                        classificationReport.groupID = source_group
                        classificationReport.spec_proprietary_period_value = proprietary_period
                        classificationReport.spec_proprietary_period_units = proprietary_units

                        pprint(classificationReport.fill(), tab='  ')
                        proceed = input("\nProceed with classification and upload? ([y]/n) : ")
                        if proceed == 'y' and not proceed.strip() == '':

                            # ASCII FILE UPLOAD
                            print ("\n")
                            response = upload_to_TNS(files)
                            print (response)

                            if not response:
                                print("File upload didn't work")
                                print(response)
                                #return False

                            print(response['id_code'], response['id_message'],
                                  "\nSuccessfully uploaded ascii spectrum")
                            #classificationReport.asciiName = response['data'][-1]

                            report_id = tns_classify(classificationReport)
                            if tns_feedback(report_id) == True:
                                post_comment(ztfname, 'Uploaded to TNS')


                    if inst == 'DIS':

                        classifiers = 'Melissa L. Graham (UW), M. Chu, A. Dahiwale, C. Fremling(Caltech) on behalf of the Zwicky Transient Facility (ZTF)'### Change accordingly
                        source_group = 48 ### Require source group id from drop down list, 0 is for None
                        spectypes = np.array(['object','host','sky','arcs','synthetic'])
                        #proprietary_period = int(input("Proprietary period in years:", x)
                        proprietary_period = '0'
                        proprietary_units = "years"
                        spec_comments =''
                        classification_comments = ''
                        spectype='object'
                        spectype_id = ['object', 'host', 'sky', 'arcs', 'synthetic'].index(spectype) + 1

                        #obsdate = APO(specid)[0]
                        #exptime = APO(specid)[1]
                        #observers = APO(specid)[2]
                        #reducers = APO(specid)[3]

                        obsdate = str(a['data']['observed_at'].split('T')[0])+' '+str(a['data']['observed_at'].split('T')[1])

                        classificationReport = TNSClassificationReport()
                        classificationReport.name = get_IAUname(ztfname)[0]['name'][3:]
                        classificationReport.fitsName = ''
                        classificationReport.asciiName = spectrum_name
                        classificationReport.classifierName = classifiers
                        classificationReport.classificationID = get_TNS_classification_ID(classify)
                        classificationReport.redshift = get_redshift(ztfname)
                        classificationReport.classificationComments = classification_comments
                        classificationReport.obsDate = obsdate
                        classificationReport.instrumentID = get_TNS_instrument_ID(inst)
                        #classificationReport.expTime = exptime
                        #classificationReport.observers = observers
                        #classificationReport.reducers = reducers
                        classificationReport.observers = (str(a['data']['observers'][0]['first_name'])+' '+str(a['data']['observers'][0]['last_name']))
                        classificationReport.reducers = (str(a['data']['reducers'][0]['first_name'])+' '+str(a['data']['reducers'][0]['last_name']))
                        classificationReport.specTypeID = spectype_id
                        classificationReport.spectrumComments = spec_comments
                        classificationReport.groupID = source_group
                        classificationReport.spec_proprietary_period_value = proprietary_period
                        classificationReport.spec_proprietary_period_units = proprietary_units

                        pprint(classificationReport.fill(), tab='  ')

                        proceed = input("\nProceed with classification and upload? ([y]/n) : ")
                        if proceed == 'y' and not proceed.strip() == '':

                            #ASCII FILE UPLOAD
                            print ("\n")
                            response = upload_to_TNS(files)
                            print (response)

                            if not response:
                                print("File upload didn't work")
                                print(response)
                                #return False

                            print(response['id_code'], response['id_message'],
                                  "\nSuccessfully uploaded ascii spectrum")
                            #classificationReport.asciiName = response['data'][-1]

                            report_id = tns_classify(classificationReport)
                            if tns_feedback(report_id) == True:
                                post_comment(ztfname, 'Uploaded to TNS')

def convert_to_jd(date):

    d = Time(date, format='fits')
    dat = d.jd
    return dat

def get_all_group_sources(group_id):
    ''' Info : Query all sources saved in a group
        Input : group id
        Returns : List of jsons of all sources in group(s)
        Comment : Takes a long time
    '''

    sources = []

    for i in range (get_number_of_sources(group_id)):

        url = BASEURL+'api/sources?saveSummary=true&group_ids='+group_id
        response = api('GET',url)
        ztfname = response['data']['sources'][i]['obj_id']
        sources.append(ztfname)

    return sources

def get_all_spectra_id(ztfname):
    ''' Info : Query all spectra corresponding to a source, takes input ZTF name
        Returns : list of spectrum jsons
    '''

    spec_id = []

    for i in range (get_all_spectra_len(ztfname)):

        url = BASEURL+'api/sources/'+ztfname+'/spectra'
        response = api('GET',url)

        specid = response['data']['spectra'][i]['id']
        spec_id.append(specid)

    return spec_id

def get_all_spectra_len(ztfname):

    url = BASEURL+'api/sources/'+ztfname+'/spectra'
    response = api('GET',url)
    return len(response['data']['spectra'])

def get_classi(ztfname):

    ''' Info : Query the classification and classification date for any source
        Input : ZTFname
        Returns : Classification and Classification date
        Comment : You need to choose the classification if there are multiple classifications
    '''

    try:

        url = BASEURL+'api/sources/'+ztfname+'/classifications'
        response = api('GET',url)
        output = response['data']

        if (len(output)< 1):
            classification = "Not Classified"
            classification_date = "None"

        if (len(output) >= 1):

            classification = response['data'][(len(output)-1)]['classification']
            classification_date = response['data'][(len(output)-1)]['created_at']

    except KeyError as e:
        classification = "Error"
        classification_date = "Error"

    return classification, classification_date

def get_classification(ztfname):

    ''' Info : Query the classification and classification date for any source
        Input : ZTFname
        Returns : Classification and Classification date
        Comment : You need to choose the classification if there are multiple classifications
    '''

    url = BASEURL+'api/sources/'+ztfname+'/classifications'
    response = api('GET',url)
    output = response['data']

    if (len(output)< 1):
        classification = "No Classification found"
        probability = "None"
        classification_date = "None"

    if (len(output)==1):

        classification = response['data'][0]['classification']
        probability = response['data'][0]['probability']
        classification_date = response['data'][0]['created_at'].split('T')[0]

        if (probability <= 0.6):

            print ("Low probability, do you want to proceed with the classification?\a")

            user_input = input("y/n: ")

            if user_input == 'y':

                classification = classification
                probability =  probability
                classification_date = classification_date

    if (len(output) > 1):

        classification = []
        classification_date = []
        probability = []

        for i in range (len(output)):

            classify = response['data'][i]['classification']
            classify_date = response['data'][i]['created_at']
            prob = response['data'][i]['probability']

            classification.append(classify)
            probability=np.append(probability, prob)
            classification_date.append(classify_date)

        for i in range (len(classification)):

            print ((i+1),")", "Classification: ", classification[i],  "\t Probability:", str(probability[i]), "\t Classification date:", classification_date[i].split('T')[0])

        user_input = input("Choose classification: \a")

        if (probability[int(user_input)-1] <= 0.6):

            print ("Low probability, do you want to proceed with the classification?\a")

            user_choice = input("y/n: ")

            if user_choice == 'y':

                classification = classification[int(user_input)-1]
                probability = probability[int(user_input)-1]
                classification_date = classification_date[int(user_input)-1].split('T')[0]
        else:

            classification = classification[int(user_input)-1]
            probability = probability[int(user_input)-1]
            classification_date = classification_date[int(user_input)-1].split('T')[0]

    return classification, probability, classification_date

def get_group_ids(groupnames=['Redshift Completeness Factor', 'Census of the Local Universe Caltech']):
    ''' Info : Query group ids of groups specified
        Input : Name or names of groups in an array []
        Returns : List of group  names and their group ids
    '''

    url = BASEURL+'api/groups'
    headers = {'Authorization': f'token {GETTOKEN}'}
    groupnames = np.atleast_1d(groupnames)
    grpids = []
    for grpname in groupnames:
        response = requests.request('GET',url,params={'name':grpname}, headers=headers).json()
        answer = str(grpname)+' = '+str(response['data'][0]['id'])
        grpids.append(answer)

    return grpids

def get_group_sources(group_id, date):
    ''' Info : Query all sources saved in a group after a certain date
        Input : group id, date [yyyy-mm-dd]
        Returns : List of jsons of all sources in group(s)
        Comment : Takes a little time based on the date
    '''

    sources = []

    for i in range (get_number_of_sources(group_id, date)):

        url = BASEURL+'api/sources?saveSummary=true&group_ids='+group_id+'&savedAfter='+date+'T00:00:00.000001'
        response = api('GET',url)
        ztfname = response['data']['sources'][i]['obj_id']
        sources.append(ztfname)

    return sources

def get_IAUname(ztfname):

    ''' Info : Query the TNS name for any source
        Input : ZTFname
        Returns : ATname
    '''

    url = BASEURL+'api/alerts_aux/'+ztfname
    response = api('GET',url)
    return response["data"]["cross_matches"]["TNS"]

def get_number(group_id, date):
    ''' Info : Query number of sources saved in a group after a certain date
        Input : group id, date [yyyy-mm-dd]
        Returns : Number of sources saved after a given date to the specified group
    '''

    url = BASEURL+'api/sources?saveSummary=true&group_ids='+group_id+'&savedAfter='+date+'T00:00:00.000001'
    response = api('GET',url)

    return len(response['data']['sources'])

def get_number_of_sources(group_id, date):
    ''' Info : Query number of sources saved in a group after a certain date
        Input : group id, date [yyyy-mm-dd]
        Returns : Number of sources saved after a given date to the specified group
    '''

    url = BASEURL+'api/sources?saveSummary=true&group_ids='+group_id+'&savedAfter='+date+'T00:00:00.000001'
    response = api('GET',url)
    return len(response['data']['sources'])

def get_pprint(item, indent=0, tab=' '*4, maxwidth=float('inf')):
    """
    it's just like 'from pprint import pprint', except instead of
    having dictionaries use hanging indents dependent on the length
    of their key, if the value is a list or dict it prints it indented
    by the current indent plus tab

    params:
        item <di or li> (the thing to be printed)
        indent <int> (the number of times it's been indented so far)
        tab <str> (how an indent is represented)
        maxwidth <int> (maximum characters per line in ouptut)

    returns:
        result <str>
    """
    def get_pprint_di(di, indent, tab=' '*4):
        """
        pprints a dictionary

        params:
            di <dict>
            indent <int> (the number of indents so far)

        returns:
            di_str <str>
        """
        di_str = ''
        for i, (key, item) in enumerate(di.items()):
            di_str += tab*indent
            di_str += repr(key) + ': ' + get_pprint(item, indent, tab)
            if i+1 < len(di):
                # everything until the last item has a trailing comma
                di_str += ',\n'
            else:
                di_str += '\n'
        return di_str

    def get_pprint_li(li, indent, tab=' '*4):
        """
        pprints a list

        params:
            li <list>
            indent <int> (the number of indents so far)

        returns:
            current_result <str>
        """
        li_str = ''
        for i, item in enumerate(li):
            li_str += tab*indent
            pprint(item, indent, tab)
            if i+1 < len(li):
                li_str += ',\n'
            else:
                li_str += '\n'
        return li_str

    result = ''
    if isinstance(item, dict):
        result += '{\n'
        result += get_pprint_di(item, indent+1, tab)
        result += tab*indent + '}'
    elif isinstance(item, list):
        result += '[\n'
        result += get_pprint_li(item, indent+1, tab)
        result += tab*indent + ']'
    else:
        result += repr(item)


    # this gets rid of too-long lines, but only supports space tabs
    lines = result.split('\n')
    for i, line in enumerate(lines):
        while max([len(li) for li in line.split('\n')]) > maxwidth:
            tabs = line[:-len(line.lstrip())]
            if len(tabs) > maxwidth - 8:
                break # giving up
            line = line[:78] + '\\\n' + tabs + 2*tab + line[78:]
            lines[i] = line
    result = '\n'.join(lines)

    return result

def get_redshift(ztfname, return_err=False):

    ''' Info : Query the redshift for any source
        Input : ZTFname
        Returns : redshift
    '''

    url = BASEURL+'api/sources/'+ztfname
    response = api('GET',url)

    redshift = response['data']['redshift']
    redshift_err = response['data']['redshift_error']

    if (redshift == None):
        redshift = "No redshift found"
    if redshift_err == None:
        redshift_err = 'No redshift error found'

    if return_err == False:
        return redshift
    else:
        return redshift, redshift_err

def get_required_spectrum_id(ztfname, auto=False):

    ''' Info : Requests spectrum from Fritz based on user selection
        Input : ZTFname, auto (if True and only one spectrum on Fritz, returns it without prompt)
        Returns : spectrum ID on Fritz
    '''

    flag = 0

    spec = (get_all_spectra_len(ztfname))

    name = []
    instids = []
    date = []

    if spec == 0:

        specid = "No Spectra Found"
        flag = 1

    if flag == 0:

        spec_id = get_all_spectra_id(ztfname)

        if auto==True and len(spec_id) == 1:
            print(ztfname + ' automatically selected spectrum with ID ' + str(spec_id[0]))
            return spec_id[0]

        for s in range (spec):

            url = BASEURL+'api/sources/'+ztfname+'/spectra'
            response = api('GET',url)

            spec_name = response['data']['spectra'][s]['original_file_filename']
            spec_date = response['data']['spectra'][s]['observed_at']
            instid = response['data']['spectra'][s]['id']

            name.append(spec_name)
            date.append(spec_date.split('T')[0])
            instids.append(instid)

        print ("Please choose from the following spectra (enter 0 to resume): \a\n")

        for i in range (len(name)):
            print ((i+1),")", "instrument id: ", instids[i], "spectrum name: ", name[i], "spectrum date:", date[i])

        wb.open(BASEURL+'source/'+ztfname, new=2)

        while True:
            user_input = input("Choose spectrum to upload: ")

            if user_input == '0':
                return 0
            else:
                try:
                    specid = spec_id[int(user_input)-1]

                    return specid
                except (IndexError, ValueError):
                    continue

    return specid

def get_source_api(ztfname):
    ''' Info : Query a single source, takes input ZTF name
        Returns : all basic data of that source (excludes photometry and spectra,
                  includes redshift, classification, comments, etc.)
    '''
    url = BASEURL+'api/sources/'+ztfname+'?includeComments=true'
    response = api('GET',url, timeout=None)
    return response['data']

def get_source_file(outfile, since):

    ''' Info : Runs sourceclassification based on 'since' parameter
        Input : outfile, date since
        Returns : None
    '''

    if since == '':
        sourceclassification(outfile) #download the updated list of sources saved to RCF in descending order
    else:
        sourceclassification(outfile, since)

def get_sources(group_id, date):
    ''' Info : Query all sources saved in a group after a certain date
        Input : group id, date [yyyy-mm-dd]
        Returns : List of jsons of all sources in group(s)
        Comment : Takes a little time based on the date
    '''

    sources = []

    for i in range (get_number(group_id, date)):

        url = BASEURL+'api/sources?saveSummary=true&group_ids='+group_id+'&savedAfter='+date+'T00:00:00.000001'
        response = api('GET',url)
        ztfname = response['data']['sources'][i]['obj_id']
        sources.append(ztfname)

    return sources

def get_spectrum_api(spectrum_id):
    ''' Info : Query all spectra corresponding to a source, takes input ZTF name
        Returns : list of spectrum jsons
    '''
    url = BASEURL+'api/spectrum/'+str(spectrum_id)
    response = api('GET',url)
    return response

def get_TNS_classification_ID(classification):

    ''' Info : Retrieves TNS classification ID based on Fritz classification
        Input : Fritz classification
        Returns : TNS ID
    '''

    class_ids = {"Other": 0, "Supernova": 1, "Type I": 2, "Ia": 3, "Ia-norm": 3, "Ib": 4, "Ib-norm": 4, "Ic": 5, "Ic-norm": 5, "Ib/c": 6, "Ic-BL": 7, "Ib-Ca-rich": 8, "Ibn": 9, "Type II": 10, "II-norm": 10, "IIP": 11, "IIL": 12, "IIn": 13, "IIb": 14,
        "I-faint": 15, "I-rapid": 16, "SLSN-I": 18, 'Ic-SLSN': 18, "SLSN-II": 19, "SLSN-R": 20, "Afterglow": 23, "LBV": 24, "ILRT": 25, "Novae": 26, "Cataclysmic": 27, "Stellar variable": 28, "AGN": 29, "Galactic Nuclei": 30, "QSO": 31, "Light-Echo": 40,
        "Std-spec": 50, "Gap": 60, "Gap I": 61, "Gap II": 62, "LRN": 65, "FBOT": 66, "kilonova": 70, "Impostor-SN": 99, "Ia-pec": 100, "Ia-SC": 102, "Ia-03fg": 102, "Ia-91bg": 103, "Ia-91T": 104, "Ia-02cx": 105,
        "Ia-CSM": 106, "Ib-pec": 107, "Ic-pec": 108, "Icn": 109, "Ibn/Icn": 110, "II-pec": 111, "IIn-pec": 112, "Tidal Disruption Event": 120, "FRB": 130, "Wolf-Rayet": 200, "WR-WN": 201, "WR-WC": 202, "WR-WO": 203, "M dwarf": 210,
        "Computed-Ia": 1003, "Computed-IIP": 1011, "Computed-IIb": 1014, "Computed-PISN": 1020, "Computed-IIn": 1021}

    #keys = np.array(class_ids.keys())
    for keys in class_ids:
        if (keys == classification):
            classkey = class_ids[keys]
            return classkey

def get_TNS_information(ztfname):

    ''' Info : Retrieves relevant info for submission to TNS from Fritz (IAU name, classification, redshift)
        Input : ZTFname
        Returns : ZTFname (redundant), IAU name, classification, redshift
    '''

    url = BASEURL+'api/sources/'+ztfname
    response = api('GET',url)

    IAU = get_IAUname(ztfname)

    if not IAU:
        IAU = "Not reported to TNS"

    else:
        IAU = IAU[0]['name']

    clas = get_classification(ztfname)

    if clas[1] == 'None':
        clas = "Not classified yet"

    else:
        clas = ('Classification: '+str(clas[0])+','+ ' Probability: '+str(clas[1])+','+' Classification date: '+str(clas[2]))

    redshift = get_redshift(ztfname)

    if redshift == None:
        redshift = "No redshift found"

    else:
        redshift = ('redshift:'+str(redshift))

    return ztfname, IAU, clas, redshift

def get_TNS_instrument_ID(inst):

    ''' Info : Retrieves the TNS instrument ID based on Fritz instrument
        Input : instrument
        Returns : TNS instrument ID
    '''

    inst_ids = {'DBSP':1, 'ALFOSC': 41, 'LRIS': 3, 'DIS': 70, 'SEDM': 149, 'SPRAT': 156, 'GMOS': 6, 'Lick-3m': 10, 'LFC': 2, 'TSPEC': 109}

    #keys = np.array(class_ids.keys())
    for keys in inst_ids:
        if (keys == inst):
            instkey = inst_ids[keys]
            return instkey

def get_TNSname(ztfname):

    ''' Info : Query the TNS name for any source
        Input : ZTFname
        Returns : ATname
    '''

    try:

        url = BASEURL+'api/alerts_aux/'+ztfname
        response = api('GET',url)
        IAU = response['data']['cross_matches']['TNS']


        if not IAU:
            IAU = "Not reported to TNS"

        else:
            IAU = IAU[0]['name']

    except KeyError as e:
        IAU = "Error"

    return IAU

def get_total_number_of_sources(group_id):
    ''' Info : Query total number of sources saved in a group
        Input : group id
        Returns : Total number of sources saved in a group
    '''

    url = BASEURL+'api/sources?saveSummary=true&group_ids='+group_id
    response = api('GET',url)
    return len(response['data']['sources'])

def get_user(fritz_id):
    resp = api('GET', BASEURL+'api/user/'+str(fritz_id))

    return resp['data']['username']

def write_ascii_file(ztfname, path=os.getcwd(), auto=False):

    ''' Info : Generates ASCII file with data from selected Fritz spectrum
        Input : ZTFname, path, auto
        Returns : spectrum_name (name of file), specid
    '''

    if auto == True:
        specid = get_required_spectrum_id(ztfname, auto=True)
    else:
        specid = get_required_spectrum_id(ztfname)

    if (specid == 'No Spectra Found'):
        spectrum_name = 'No Spectra Found'
        print (ztfname + ' has no spectrum on Fritz.')

        return spectrum_name, specid

    if specid == 0:
        spectrum_name = 'Resuming...'

        return spectrum_name, specid

    a = get_spectrum_api(specid)

    inst = (a['data']['instrument_name'])
    #print (inst)

    if inst == 'SEDM':

        header = (a['data']['altdata'])


        s = (ztfname+'_'+str(header['OBSDATE'])+'_'+str(inst)+'.ascii')

        with open(path+'/data/'+s,'w') as f:
            f.write(a['data']['original_file_string'])
        f.close()

        #print (s,'\n')
        spectrum_name = s


    elif inst == 'SPRAT':


        header = (a['data']['altdata'])



        s = (ztfname+'_'+str(header['OBSDATE'].split('T')[0])+'_'+str(inst)+'.ascii')

        with open(path+'/data/'+s,'w') as f:
            f.write(a['data']['original_file_string'])
        f.close()

        #print (s,'\n')
        spectrum_name = s


    elif inst == 'ALFOSC':

        OBSDATE = a['data']['observed_at'].split('T')[0]


        s = (ztfname+'_'+str(OBSDATE)+'_'+str(inst)+'.ascii')

        with open(path+'/data/'+s,'w') as f:
            f.write(a['data']['original_file_string'])
        f.close()

        #print (s,'\n')
        spectrum_name = s


    elif inst == 'DBSP':

        wav = (a['data']['wavelengths'])
        flux = (a['data']['fluxes'])
        err = (a['data']['errors'])

        OBSDATE = a['data']['observed_at'].split('T')[0]



        s = (ztfname+'_'+str(OBSDATE)+'_'+str(inst)+'.ascii')

        if err == None:

            with open(path+'/data/'+s,'w') as f:

                for i in range(len(wav)):
                    f.write(str(wav[i])+'\t'+str(flux[i])+'\n')
            f.close()

            #print (s,'\n')
            spectrum_name = s

        else:

            with open(path+'/data/'+s,'w') as f:

                for i in range(len(wav)):
                    f.write(str(wav[i])+'\t'+str(flux[i])+'\t'+str(err[i])+'\n')
            f.close()

            #print (s,'\n')
            spectrum_name = s


    elif inst == 'DIS':

        #obsdate = a['data']['original_file_string'].split('#')[6]
        #a,b = obsdate.split(' ', 1)
        #c,OBSDATE = b.split(' ', 1)
        #OBSDATE = OBSDATE.split('T')[0]

        obsdate = a['data']['observed_at']
        OBSDATE = obsdate.split('T')[0]


        path = path+'/data/'

        s = (ztfname+'_'+str(OBSDATE)+'_'+str(inst)+'.ascii')

        a = get_spectrum_api(specid)

        with open(path+s,'w') as f:
            f.write(a['data']['original_file_string'])
        f.close()

        #print (s,'\n')
        spectrum_name = s


    elif inst == 'LRIS':

        wav = (a['data']['wavelengths'])
        flux = (a['data']['fluxes'])
        err = (a['data']['errors'])

        OBSDATE = a['data']['observed_at'].split('T')[0]



        s = (ztfname+'_'+str(OBSDATE)+'_'+str(inst)+'.ascii')

        if err == None:

            with open(path+'/data/'+s,'w') as f:

                for i in range(len(wav)):
                    f.write(str(wav[i])+'\t'+str(flux[i])+'\n')
            f.close()

            #print (s,'\n')
            spectrum_name = s

        else:

            with open(path+'/data/'+s,'w') as f:

                for i in range(len(wav)):
                    f.write(str(wav[i])+'\t'+str(flux[i])+'\t'+str(err[i])+'\n')
            f.close()

            #print (s,'\n')
            spectrum_name = s


    elif inst == 'EFOSC2':

        spectrum_name = 'TNS_spectrum'

    else:
        spectrum_name = None

    return spectrum_name, specid

def post_comment(ztfname, text, attach=None, attach_name=None):

    ''' Info : Posts a comment on transient's Fritz page
        Input : ZTFname, text
        Returns : API response
    '''

    data = {
            "text": text,
           }

    if attach != None:
        with open(attach, "rb") as img_file:
            at_str = base64.b64encode(img_file.read()).decode('utf-8')

        data['attachment'] = {'body': at_str, 'name': attach_name}

    url = BASEURL+'api/sources/'+ztfname+'/comments'

    response = api('POST', url, data=data)

    return response

def edit_comment(ztfname, comment_id, author_id, text, attach=None, attach_name=None):

    ''' Info : Posts a comment on transient's Fritz page
        Input : ZTFname, text
        Returns : API response
    '''

    data = {  "obj_id": ztfname,
              "text": text,
              'author_id': author_id
           }

    if attach != None:
        with open(attach, "rb") as img_file:
            at_str = base64.b64encode(img_file.read()).decode('utf-8')

        data['attachment_name'] = attach_name
        data['attachment_bytes'] = at_str

    url = BASEURL+'api/sources/'+ztfname+'/comments/'+str(comment_id)

    response = api('PUT', url, data=data)

    return response

def pprint(*args, **kwargs):
    """
    slightly more convenient function instead of print(get_pprint)

    params:
        *args (arguments to pass to get_pprint)
        **kwargs (keyword arguments to pass to get_pprint)
    """
    print(get_pprint(*args, **kwargs))

def read_ascii(f, startd):

    ''' Info : Reads ASCII table for classified or saved transients passed specified date
        Input : ASCII table, earliest date to filter
        Returns : sources (ZTF names), dates saved, classifications, classification dates, unclassified transients, redshifts
    '''

    sources_r = np.asarray(f['Source Name'])
    tns_names_r = np.asarray(f['TNS Name'])
    savedates_r = np.asarray(f['Saved Date'])
    classifys_r = np.asarray(f['Classification'])
    class_dates_r = np.asarray(f['Classification Date'])
    reds_r = np.asarray(f['redshift'])

    sources = np.array([])
    tns_names = np.array([])
    savedates = np.array([])
    classifys = np.array([])
    class_dates = np.array([])
    reds = np.array([])
    unclassifys = np.array([])

    for i in np.arange(0,len(sources_r)):
        if classifys_r[i] != 'Not Classified' and (datetime.datetime.strptime(class_dates_r[i], '%Y-%m-%d').replace(tzinfo=datetime.timezone.utc) >= startd or datetime.datetime.strptime(savedates_r[i], '%Y-%m-%d').replace(tzinfo=datetime.timezone.utc) >= startd):
            sources = np.append(sources, sources_r[i])
            tns_names = np.append(tns_names, tns_names_r[i])
            savedates = np.append(savedates, savedates_r[i])
            classifys = np.append(classifys, classifys_r[i])
            class_dates = np.append(class_dates, class_dates_r[i])
            reds = np.append(reds, str(reds_r[i]))
        if classifys_r[i] == 'Not Classified' and datetime.datetime.strptime(savedates_r[i], '%Y-%m-%d').replace(tzinfo=datetime.timezone.utc) >= startd:
            unclassifys = np.append(unclassifys, sources_r[i])

    return sources, tns_names, savedates, classifys, class_dates, reds, unclassifys

def sourceclassification(outfile, dat=str(datetime.datetime.utcnow().date() - datetime.timedelta(days=180))):

    ''' Info : Downloads list of transients on Fritz saved after specified date (or since 180 days prior if no input)
               Saves ZTF names, TNS names, dates saved, classifications, classifications, redshifts as ASCII file
        Input : outfile name, date to check after
        Returns : None
    '''

    path = 'https://fritz.science/api/sources?group_ids=41&saveSummary=true&savedAfter='+str(dat)+'T00:00:00.000001'

    response = api('GET',path)

    #print('data received')

    srcs = []
    dates = []
    classify = []
    class_date = []
    TNS = []
    reds = []

    listdir = os.getcwd()
    f = open (listdir+'/'+outfile+'.ascii','w')
    f.write('Source Name'+'\t'+'TNS Name'+'\t'+'Saved Date'+'\t'+'Classification'+'\t'+'Classification Date'+'\t'+'redshift'+'\n')

    num_rcf = range (get_number('41', dat))

    for i in tqdm(num_rcf, desc='Downloading sources...'):

        source_name = response['data']['sources'][i]['obj_id']
        saved_date = response['data']['sources'][i]['saved_at']
        classification = get_classi(source_name)[0]
        date = get_classi(source_name)[1]
        IAU = get_TNSname(source_name)
        red = str(get_redshift(source_name))

        #print (i, source_name)

        srcs.append(source_name)
        TNS.append(IAU)
        dates.append(saved_date.split('T')[0])
        classify.append(classification)
        class_date.append(date.split('T')[0])
        reds.append(red)

    output = sorted(zip(class_date, srcs, TNS, dates, classify, reds), reverse=True)

    for i in num_rcf:

        f.write(output[i][1]+'\t'+output[i][2]+'\t'+output[i][3]+'\t'+output[i][4]+'\t'+output[i][0]+'\t'+output[i][5]+'\n')

    f.close()

def submit_fritz_class(ztfname, clas):

    ''' Info : Uploads classification to Fritz
        Input : ZTFname, classification
        Returns : API response
    '''

    data = {    "obj_id": ztfname,
                "classification": clas,
                "taxonomy_id": 3,
                "probability": 1}

    url = BASEURL+'api/classification'

    response = api('POST', url, data=data)

    return response

def submit_fritz_redshift(ztfname, redshift,redshift_err):

    ''' Info : Uploads redshift to Fritz
        Input : ZTFname, redshift
        Returns : API response
    '''

    data = {    "redshift": redshift,
                "redshift_error": redshift_err}

    url = BASEURL+'api/sources/'+ztfname

    response = api('PATCH', url, data=data)

    return response

def tns_classify(classificationReport, base_url= report_url, api_key=API_KEY):

    ''' Info : Uploads classification report to TNS
        Input : classification report, transient URL on TNS, API KEY
        Returns : API response
    '''

    url = base_url

    headers={'User-Agent':'tns_marker{"tns_id":'+str(YOUR_BOT_ID)+', "type":"bot", "name":"'+YOUR_BOT_NAME+'"}'}

    data = {'api_key' : api_key, 'data' : classificationReport.as_json()}
    response = requests.post(url, headers=headers, data=data).json()
    if not response:
        return False

    res_code = response['id_code']
    report_id = response['data']['report_id']
    print("ID:", report_id)
    print(str(res_code) + ' ' + response['id_message'] + ' ' + "reporting finished")
    if res_code == 200:
        return report_id
    else:
        print("Result reporting didn't work")
        pprint(response)
        print("re-submit classification, but don't re-upload files")
        return False

def tns_feedback(report_id):

    ''' Info : Verifies that report has been uploaded
        Input : ID of report
        Returns : API response
    '''

    data = {'api_key': API_KEY, 'report_id': report_id}

    headers={'User-Agent':'tns_marker{"tns_id":'+str(YOUR_BOT_ID)+', "type":"bot", "name":"'+YOUR_BOT_NAME+'"}'}

    response = requests.post(reply_url, headers=headers, data=data).json()

    feedback_code = response['id_code']
    print(feedback_code, response['id_message'], "feedback finished")
    if feedback_code == 200:
        print(bcolors.OKGREEN + 'Feedback successful. Continuing...' + bcolors.ENDC)
        return True
    elif feedback_code == 404:
        print(bcolors.WARNING + "Waiting and retrying..." + bcolors.ENDC)
        sleep(2)
        try:
            return tns_feedback(report_id)
        except KeyboardInterrupt:
            return False
    elif feedback_code == 400:
        print(bcolors.FAIL + json.dumps(response, indent=2) + bcolors.ENDC)
        return False
    else:
        # error receiving the feedback from TNS about the upload
        print("Something went wrong with the feedback, but the report may",
              "still have been fine?")
        return False

def upload_to_TNS(filename, base_url = upload_url, api_key = API_KEY, filetype='ascii'):

    ''' Info : Uploads spectrum to TNS
        Input : spectrum file name, transient URL on TNS, API KEY, spectrum file type
        Returns : API response
    '''

    url = base_url
    data = {'api_key' : api_key}

    headers={'User-Agent':'tns_marker{"tns_id":'+str(YOUR_BOT_ID)+', "type":"bot", "name":"'+YOUR_BOT_NAME+'"}'}

    if filetype == 'ascii':
        files = [('files[]', (filename, open(filename), 'text/plain'))]

    elif filetype == 'fits':
        files = [('files[0]', (filename, open(filename, 'rb'),
                               'application/fits'))]

    if filename:
        response = requests.post(url, headers=headers, data=data, files=files)
        try:
            return response.json()
        except:
            print(url, data, files, response.content, sep='\n')
            return False
    else:
        return {}

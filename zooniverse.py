from panoptes_client import Panoptes, Project, SubjectSet, Subject, Classification, Workflow
import csv
import json
import datetime
import sys
import os
import numpy as np
from scipy import stats
import requests
from PIL import Image
from func import *
from pprint import pprint
import time
from superfit_func import *

with open('info.info', 'r') as infofile:
    info = infofile.read()
    zoo_user = info.split('\n')[4].split(':')[1].strip()
    zoo_pass = info.split('\n')[5].split(':')[1].strip()

def pull_class(startd):

    old_ims = os.listdir(os.getcwd() + '/zooniverse')
    for o in old_ims:
        os.remove(os.getcwd() + '/zooniverse/' + o)

    Panoptes.connect(username=zoo_user, password=zoo_pass)

    project = Project.find(12959)

    subject_set = SubjectSet()

    subject_set.links.project = project
    subject_set.display_name = 'Newly Unclassified'

    subject_set = SubjectSet.find(99282)

    news = []
    image_urls = []
    sub_ids = []

    for subject in subject_set.subjects:
        if subject.subject_workflow_status(16969).retired_at == None:
            continue

        if datetime.datetime.strptime(subject.subject_workflow_status(16969).retired_at, '%Y-%m-%dT%H:%M:%S.%fZ').replace(tzinfo=datetime.timezone.utc) > startd:
            #print(subject.id)
            sub_ids.append(subject.id)
            news.append(subject.metadata['!ZTF_Name'])
            #pprint(subject.locations)
            image_urls.append(subject.locations)

    subject_set.save()

    cresult = Classification.where(scope='project', project_id=12959, workflow_id=16969, last_id=373078596)
    df = pd.DataFrame([c.raw for c in list(cresult)])

    classifications = []
    class_sources = []
    class_users = []

    startd = datetime.datetime.strptime('2021-11-10', '%Y-%m-%d').replace(tzinfo=datetime.timezone.utc)

    for d in range(len(df)):
        if df['links'][d]['subjects'][0] in np.array(sub_ids):
            class_sources.append(news[int(np.argwhere(np.array(sub_ids) == df['links'][d]['subjects'][0]))])
            classifications.append(df['annotations'][d][0]['value'])
            class_users.append(df['links'][d]['user'])

    classifications = [int(0 if value is None else value) for value in classifications]

    sources = np.unique(class_sources)

    for n, new in enumerate(sources):
        print(bcolors.OKCYAN + str(n+1) + '/' + str(len(sources)) + ': ' + bcolors.ENDC + bcolors.OKBLUE + new + bcolors.ENDC)

        un_class = np.array(classifications)[np.array(class_sources) == new]
        un_class_users = np.array(class_users)[np.array(class_sources) == new]
        #print(un_class)
        #print(un_class_users)
        #print(stats.mode(un_class).mode)

        images_available = False

        for m in stats.mode(un_class).mode:
            if m > 0:
                images_available = True
                image_url = image_urls[int(np.argwhere(np.array(news) == new))][m-1]['image/png']
                img_data = requests.get(image_url).content
                with open('zooniverse/' + new + '.png', 'wb') as handler:
                    handler.write(img_data)
                image = Image.open('zooniverse/' + new + '.png')
            elif m == 0:
                print('No good match')
            else:
                print('Issues')

        if images_available:
            comment_infos = get_source_api(new)['comments']

            uploaded = False
            zoo_class = False
            for i in range (len(get_source_api(new)['comments'])):

                comment_info = comment_infos[i]
                comment = comment_info['text']

                if 'Uploaded to TNS' in comment:
                    uploaded = True

                if 'zooniverse classification' in comment:
                    zoo_class = True

            if uploaded == True:
                print(new + ' has already been uploaded to TNS.')
                continue

            if zoo_class == True:
                print(new + ' has already been classified from Zooniverse.')
                continue

            width, height = image.size
            image = image.resize((width//5, height//5))
            image.show()

            upload = input('Enter in the name of the best classification: ')

            if upload == 'II':
                upload = 'Type II'
            elif upload == 'Gal':
                upload = 'Galactic Nuclei'
            elif upload == 'Ia-csm':
                upload = 'Ia-CSM'

            if 'II' not in upload:
                resp = post_comment(new, 'zooniverse classification: ' + upload + ', ' + str(stats.mode(un_class).count) + '/' + str(len(un_class)) +
                    ' classifications', 'zooniverse/'+new+'.png', new+'_zooniverse.png')

                continue

            print('Running superfit...')

            run_superfit(new)

            match = input('Does the superfit classification match SNID? [y/n] ')

            if match != 'y':
                continue

            pre_class = get_classi(new)[0]
            if pre_class == upload:
                print(new + ' already classified with the same classification on Fritz.')
                #subject_set.remove(news_ids[n])
                continue
            elif pre_class != 'Not Classified':
                if input(new + ' already classified with classification ' + pre_class + '. Submit another? [y/n] ') == 'y':
                    pass
                else:
                    continue

            fritz_class = submit_fritz_class(new, upload)

            if fritz_class['status'] == 'success':
                print(bcolors.OKGREEN + new + ' classification upload successful.' + bcolors.ENDC)
            else:
                print(bcolors.FAIL + new + ' classification upload failed.' + bcolors.ENDC)
                print(bcolors.FAIL + fritz_class['message'] + bcolors.ENDC)

            resp = post_comment(new, 'zooniverse classification: ' + upload + ', ' + str(stats.mode(un_class).count) + '/' + str(len(un_class)) +
                ' classifications', 'zooniverse/'+new+'.png', new+'_zooniverse.png')

            if resp['status'] == 'success':
                print(bcolors.OKGREEN + new + ' comment upload successful.' + bcolors.ENDC)
            else:
                print(bcolors.FAIL + new + ' comment failed.' + bcolors.ENDC)
                print(bcolors.FAIL + json.dumps(resp, indent=2) + bcolors.ENDC)

            #subject_set.remove(news_ids[n])

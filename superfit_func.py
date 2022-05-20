from func import *
import json
import sys
import os
import time

home = os.getcwd()

with open('info.info', 'r') as infofile:
    info = infofile.read()
    superfit_loc = info.split('\n')[1].split(':')[1].strip()

sys.path.insert(1, superfit_loc)

os.chdir(superfit_loc)

import run

os.chdir(home)

def run_superfit(source):

    ''' Info : Runs Superfit on given source
        Input : Source name
        Returns : None
    '''

    os.chdir(superfit_loc)

    if 'data' not in os.listdir(superfit_loc):
        os.mkdir('data')

    fname = write_ascii_file(source, path=superfit_loc, auto=False)[0] # Downloads spectrum data in ASCII from Fritz
    redshift = get_redshift(source)

    if fname == 'No Spectra Found' or fname == 'Resuming...': # Return None if no spectrum on Fritz or if user prompts to continue
        os.chdir(home)
        return

    with open('parameters.json') as f: # Load parameters for superfit into dictionary, then modify for given source
        parameters = json.load(f)

    parameters['object_to_fit'] = 'data/' + fname

    if redshift != 'No redshift found':
        parameters['use_exact_z'] = 1
        parameters['z_exact'] = redshift
    else:
        parameters['use_exact_z'] = 0

    with open('parameters.json', 'w') as f: # Send updated dict back into json file
        json.dump(parameters, f, indent=4)

    time.sleep(1) # Not sure if this is necessary, but wait for file to be written. Probably can get rid of this

    run.run()

    os.chdir(home)

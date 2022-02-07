from func import *
import json
import sys
import os

home = os.getcwd()

with open('info.info', 'r') as infofile:
    info = infofile.read()
    superfit_loc = info.split('\n')[1].split(':')[1].strip()

sys.path.insert(1, superfit_loc)

os.chdir(superfit_loc)

import run

os.chdir(home)

def run_superfit(source):

    os.chdir(superfit_loc)

    fname = write_ascii_file(source, path='superfit', auto=False)[0] # Downloads spectrum data in ASCII from Fritz
    redshift = get_redshift(source)

    if fname == 'No Spectra Found' or fname == 'Resuming...': # Return None if no spectrum on Fritz or if user prompts to continue
        return

    with open('parameters.json') as f:
        parameters = json.load(f)

    parameters['object_to_fit'] = 'superfit/data/' + fname

    if redshift != 'No redshift found':
        parameters['use_exact_z'] = 1
        parameters['z_exact'] = redshift
    else:
        parameters['use_exact_z'] = 0

    with open('parameters.json', 'w') as f:
        json.dump(parameters, f, indent=4)

    run.run()

    os.chdir(home)

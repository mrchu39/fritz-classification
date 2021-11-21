from func import *
import run
import json

def run_superfit(source):
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

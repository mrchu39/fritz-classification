# fritz-classification
This tool pulls data from Fritz, uses SNID to classify transients, and uploads classification data to TNS.

## Requirements

* [Anaconda for Python 3.7+](https://www.anaconda.com/products/individual)
* [astropy](https://www.astropy.org/)
* [pandas](https://pandas.pydata.org/)
* [SNID](https://people.lam.fr/blondin.stephane/software/snid/index.html)

## Installation

Installation of SNID is detailed on the above link. Ensure that its dependencies (PGPLOT) are installed. Enter in the location the SNID executable is located in `snid.py`.

Generate your unique Fritz API token by going to your [profile](https://fritz.science/profile), finding "Generate New Token for Command-Line Authentication", entering in a name (this is not important), checking all the boxes, and clicking "generate token". Copy the token and enter it in `func.py`.

You will also need the ID and API key for the bot that submits information to TNS. Email mrchu@caltech.edu for them and enter them in `func.py`

## Usage

This script should be run about daily. Run `python master.py` in terminal.

### Data Download

The script will prompt the user to download an ASCII file containing a list of sources saved since the inputted date:

```
Download new list of RCF sources? ([y]/n)
```

Each source is saved with its TNS name, save date, classification (if available), classification date (if available), and redshift (if available). By default, the last 180 days of saved sources will be saved. This is necessary because some sources are saved much earlier and potentially classified recently. This process can take a while, so be patient. This data will be saved as `RCF_sources.ASCII`. If the code fails for whatever reason, you can run from the beginning but skip downloading the ASCII file by entering `n` to the prompt.

The code will then prompt to enter in a date:

```
Enter in the earliest date you want to check classifications or saves (YYYY-MM-DD),
'y' for yesterday at midnight, or enter nothing to select last time program was run:
```

In most cases, you will want to look for newly classified or saved transients from the previous day, so entering `y` is generally acceptable (if more than one day since last running, such as Monday, enter in the date of last use). However, since the code is often run multiple times in one day, the last time the program was run is often too recent. The code also saves the last time it was run, so if you want to choose this enter in nothing. If for whatever reason you want to look at objects saved earlier, you can enter in this date in the appropriate format.

### Redshift Determination for Classified Sources

The script will first search for sources saved since the inputted date that are classified but do not have redshifts in Fritz. If any exist, the user can proceed to run SNID on these sources to determine redshift. If the source has a spectrum on Fritz, the user will be prompted to enter in a spectrum to use in SNID. The Fritz page will open in your browser. Some helpful things to check are:

* **The spectra themselves**: check the dates of the uploaded spectra, typically a classification will occur a day or so after the spectrum was uploaded. If there are duplicate uploaded spectra, it does not matter. Generally, the "constep" spectra have lower SNR and were used for classification. If there is none in this time frame, check below.
* **Classifications**: Since these sources have classifications, check the date they were classified, and choose the spectrum that led to classification.
* **Comments**: Check comments for information on the source, often there will be information about what led to classification. Check these to ensure the correct spectrum is selected, sometimes they will explicitly mention when they classified. The SEDMbot will also upload reports from observing runs. If the classification came from TNS, this will likely be conveyed in the comments, and if so, none of the spectra led to classification. If this is the case, select the spectrum based on the SNID criteria in the next set of bullet points.

Once the correct spectrum is uploaded, SNID will run on the uploaded spectrum, at which point it will either converge on certain SN types or not. If it does not converge, the code will move on to the next source. If it does, it will output the ten templates with the highest `rlap` score. Of these, check the top few and identify if they match with the Fritz classification, and if they do, enter `y` when prompted by the code to save. After running through all sources that are classified without redshifts, the code will ask whether you want to upload the ones you saved to Fritz. If you upload, the code will output the API response to verify whether or not they have successfully been uploaded to Fritz.

### SNID Analysis on Unclassified Sources

After completing this, the script will select all unclassified sources that have been saved since the inputted date. The user can proceed to run SNID on unclassified transients, the process of which is the same as with the redshifts. Select the spectrum to use carefully using the following critera:

* **The spectra themselves**: check the plots on Fritz. Generally, SEDM spectra which were taken most recently will have lower SNR, and "constep" spectra will generally also have lower SNR than "robot" spectra. If there are spectra from more powerful telescopes than SEDM (such as LRIS or NOT), use these.
* **Comments**: Check comments for information on the source, sometimes there will be discussion on several potential classifications. If this is the case, SNID will likely not converge with a high `rlap` on any specific template.

The spectral flux data will be downloaded from Fritz as an ASCII fileand saved in `/data` (this directory will be generated if it does not exist). SNID will then run on the downloaded spectrum and pull this spectrum from the folder. If SNID converges, the script will open up ten plots, those with the highest `rlap` score. Check the rlap score of the first few, if they are greater than 10 and all converge on the same classification, then consider submitting a classification to Fritz. Check the plots to ensure that the spectrum closely matches that of the templates. If ever in doubt, **do not submit a classification**, further observing runs can provide better spectra.

The results of SNID will be saved in `/outfiles/<ZTFname>` (this directory also will be generated if it does not exist). These include the `.output` file, which includes the converged templates and their `rlap` scores, plus the redshift and redshift error from the fit of that spectrum. The spectral flux data of the ten templates will also be saved, along with plots of their spectra plotted over the source spectrum.

Again, at the end of the source list the user will be prompted to upload the classifications to Fritz. It will also upload redshifts for the transients if they are not already on Fritz, but will not overwrite if it does already have them. The API responses for each upload will again be returned. In the case that a SNID classification is not an option in Fritz, it will return a failure code, but in this situation just ignore the source and move on. The successfully uploaded classifications and redshifts will be updated in the source ASCII file, so it does not need to be downloaded again to include the updated classifications.

### TNS Submission

The code will then move on to submit to TNS any classified transients that have not been previously submitted by ZTF.

The code will run through the sources labeled as "classified" in the ASCII file. Some may not have spectra on Fritz, if this is the case then they were classified externally from publications or by other groups. If this is the case, there is no need to submit. However, if there is a spectrum, it will first check the comments to ensure a report has not already been submitted by ZTF. If it was, typically one of the comments will say "Uploaded to TNS" verbatim. If not, it will then check the transient's TNS page and see if there is a report from ZTF. If either of the two are found, it will move on to the next source.

If not, it will again prompt the user to enter in a spectrum. Like with the redshift analysis, we want to select the spectrum that led to classification. Follow the same criteria to ensure you selected the correct spectrum. Similarly, we want to check the comments to see if the classification came from TNS. If it did, enter 0 to proceed onto the next source. **We do not want to send classification reports if our classification came from TNS**. We send reports if we independently classify a transient.

Once you select the correct spectrum, the code will generate a TNS report, which will be submitted to TNS through the API. The script might indicate that the source has already been classified and submitted to TNS by a different group. You can check whether the classification is the same, but regardless we want to send our own reports with our own classification analyses to TNS. The code will print out the report as a dictionary, check to make sure that things look correct. If the Fritz classification has no corresponding classification on TNS or there is no redshift available in Fritz, it will return an error in the API response. Do not worry, and move on to the next source. If the submission was successful, it will print a 200 success message in green. The "Uploaded to TNS" comment will then be posted on the source's Fritz page, and the code will move on to the next source.

After completing all in the list, the script will indicate that the submission process is complete.

## Acknowledgements

A majority of the TNS submission code has been provided for use by Aishwarya Dahiwale from her own [repository](https://github.com/Aishwarya113/TNS-Classification).

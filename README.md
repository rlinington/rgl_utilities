# rgl_utilities

Basic scripts for lab data management functions

# func001_extractor.py

This tool extracts func001.csv files from the MSeXpress output directory
This is useful because data from each sample in MSeXpress is exported as a set of files in a directory that matches the sample name. Every file within that directory has standard nomenclauture (e.g. func001.csv)
The output directory contains renamed func001.csv files that are named using the directory name

To use it type:
'''
python func001_extractor.py
'''

The script will request source and destination directories.  
The source directory is the parent MSeXpress directory containing the sample directories for all the processed samples
The destination directory is the place you want the func001 files to end up.

## On Mac OSX:
Select the directory in the finder window, right click, press and hold the 'option' key and select 'Copy "DIR" as Pathname'
paste in to terminal window and hit return

## On Windows:
Select the directory in the Explorer window, right click, and select 'copy as path'
paste in to the terminal window and hit return


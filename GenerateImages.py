#!/usr/bin/env python

##########################################################################################
#
# File: GenerateImages.py
# Author: Matthew Leeds <mwleeds@crimson.ua.edu>
# Last Edit: 2015-03-29
# Purpose: Generate color RGB images in JPEG format from Sunrise's broadband
# output using Patrik's make_color module. Please run this:
# hg clone https://bitbucket.org/lutorm/python
# rename that folder 'patrikpython', add an empty '__init__.py' file, and put
# it in your PYTHONPATH. Also have numpy and astropy installed.
# Config parameters should be set in 'config.ini' under 'Generate Images'.
# Broadband output files are expected to be in the format:
# SUNRISE_DIR/GALAXY_NAME-TIME_STEP-sunrise/
# GALAXY_NAME.TIME_STEP.DIAMETERkpc.phys|sim.broadband.fits
# Image files will be written to 
# SUNRISE_DIR/GALAXY_NAME-TIME_STEP-sunrise/ (or IMAGE_DIR/)
# <galaxy name>.<time step>.<cut diameter>kpc.phys|sim-filter<i>-camera<i>[-redshift].jpg
#
##########################################################################################

from os import chdir, listdir, curdir, pardir, rename, sep
from sys import stdout, stderr, exit
from re import match
from astropy.io import fits
from patrikpython import make_color
from numpy import transpose
from ConfigParser import SafeConfigParser

CONFIG_FILE = "config.ini"
SECTION_NAME = "Generate Images"

def main():
    config = SafeConfigParser()
    stdout.write("Reading " + CONFIG_FILE + "\n")
    config.read(CONFIG_FILE)
    if not config.has_section(SECTION_NAME):
        stderr.write("Error: No '" + SECTION_NAME + "' section found in " + CONFIG_FILE + ".\n")
        exit(1)
    GALAXY_NAME = config.get(SECTION_NAME, "GALAXY_NAME") 
    SUNRISE_DIR = config.get(SECTION_NAME, "SUNRISE_DIR")
    if len(SUNRISE_DIR) == 0: SUNRISE_DIR = curdir + sep
    chdir(SUNRISE_DIR)
    listOfTimesteps = []
    TIME_STEPS = config.get(SECTION_NAME, "TIME_STEPS")
    if len(TIME_STEPS) > 0: 
        listOfTimesteps = TIME_STEPS.split(",")
    # If a time step wasn't specified, look for them.
    if len(listOfTimesteps) == 0:
        for f in listdir(curdir):
            if match(r"^" + GALAXY_NAME + "-\d+-sunrise$", f) != None:
                listOfTimesteps.append(f.split("-")[1])
    IMAGE_DIR = config.get(SECTION_NAME, "IMAGE_DIR")
    if len(IMAGE_DIR) == 0: IMAGE_DIR = curdir + sep
    FILTER_NAMES = config.get(SECTION_NAME, "FILTER_NAMES")
    AUTOPERCENTILE = config.getfloat(SECTION_NAME, "AUTOPERCENTILE")
    # Make images for every perspective and filter set combination for every time step.
    for timeStep in listOfTimesteps:
        folderName = GALAXY_NAME + "-" + timeStep + "-sunrise"
        chdir(folderName)
        # Find the broadband FITS file for this time step, and the redshifted one if it's there.
        fitsFiles = []
        for f in listdir(curdir):
            if match(r"^" + GALAXY_NAME + "\." + timeStep + "\.\d+kpc\.(phys|sim)\.broadband(-redshift)?\.fits$", f) != None:
                fitsFiles.append(f)
        if len(fitsFiles) == 0:
            stderr.write("Error: unable to find broadband FITS file for " + GALAXY_NAME + "-" + timeStep + "\n")
            continue
        for fitsFilename in fitsFiles:
            redShift = ("redshift" in fitsFilename)
            stdout.write("Reading " + fitsFilename + "\n")
            f = fits.open(fitsFilename)
            numCameras = f["MCRX"].header["N_CAMERA"]
            filters = [f["FILTERS"].data[i][0] for i in range(len(f["FILTERS"].data))]
            f.close()
            filterSets = []
            filterNameSets = FILTER_NAMES.split(";")
            # Convert the filter names to indices.
            for filterNames in filterNameSets:
                filterIndices = [str(filters.index(filterName + ".res")) for filterName in filterNames.split(",")]
                filterSets.append(",".join(filterIndices))
            stdout.write("Generating RGB images for " + str(numCameras) + " perspectives with " + str(len(filters)) + " filters.\n")
            # For each specified set of filters, make images from every perspective.
            for i, filterSet in enumerate(filterSets):
                outFile = (fitsFilename[:-24] if redShift else fitsFilename[:-15])
                outFile += "-filter" + str(i) + "-camera%d"
                outFile += ("-redshift.jpg" if redShift else ".jpg")
                make_color.write_all(fitsFilename, "-BROADBAND", IMAGE_DIR + outFile, band=filterSet, 
                                     scale="auto", autopercentile=AUTOPERCENTILE, overwrite=True)
        chdir(pardir)
    exit(0)

if __name__=="__main__":
    main()

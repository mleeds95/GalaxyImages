#!/usr/bin/env python

###############################################################################
#
# File: GenerateImages.py
# Author: Matthew Leeds <mwleeds@crimson.ua.edu>
# License: GNU GPL v3 <gnu.org/licenses>
# Last Edit: 2015-04-14
# Purpose: Generate color RGB images in JPEG format from Sunrise's broadband
# output using Patrik's make_color module. Please run this:
# hg clone https://bitbucket.org/lutorm/python
# rename that folder 'patrikpython', add an empty '__init__.py' file, and put
# it in your PYTHONPATH. Also have numpy, astropy, and PIL installed.
# Config parameters should be set in 'config.ini' under 'Generate Images'.
# Broadband output files are expected to be in the format:
# SUNRISE_DIR/GALAXY_NAME-TIME_STEP-sunrise/
# <galaxy name>.<time step>.h<halo id>.<cut diameter>kpc.phys|sim.broadband[-redshift].fits
# Image files will be written to 
# SUNRISE_DIR/GALAXY_NAME-TIME_STEP-sunrise/ (or IMAGE_DIR/)
# <galaxy name>.<time step>.h<halo id>.<cut diameter>kpc-filter<i>-camera<i>-ap<autopercentile>[-redshift].jpg
#
###############################################################################

from os import chdir, listdir, curdir, pardir, rename, sep
from sys import stdout, stderr, exit
from re import match
from astropy.io import fits
from patrikpython import make_color
from numpy import transpose
from ConfigParser import SafeConfigParser
from math import degrees

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
    AUTOPERCENTILES = config.get(SECTION_NAME, "AUTOPERCENTILES").split(",")
    failedTimesteps = []
    # Make images for every perspective, filter set, and autopercentile combination for every time step.
    for timeStep in listOfTimesteps:
        folderName = GALAXY_NAME + "-" + timeStep + "-sunrise"
        chdir(folderName)
        # Go ahead and make sure we have write access to the output directory.
        if not access(IMAGE_DIR, os.W_OK):
            stderr.write("Error no write access to " + os.path.abspath(IMAGE_DIR) + ". Check the permissions.\n")
            exit(1)
        # Find the broadband FITS file for this time step, and the redshifted one if it's there.
        fitsFiles = []
        matchRE = r"^" + GALAXY_NAME + "\." + timeStep + "\.h\d+\.\d+kpc\.(phys|sim)\.broadband(-redshift)?\.fits$"
        for f in listdir(curdir):
            if match(matchRE, f) != None:
                fitsFiles.append(f)
        if len(fitsFiles) == 0:
            stderr.write("Error: unable to find broadband FITS file for " + GALAXY_NAME + "-" + timeStep + "\n")
            failedTimesteps.append(timeStep)
            chdir(pardir)
            continue
        for fitsFilename in fitsFiles:
            # Read filter and campos data from the FITS file.
            redShift = ("redshift" in fitsFilename)
            stdout.write("Reading " + fitsFilename + "\n")
            f = fits.open(fitsFilename)
            numCameras = f["MCRX"].header["N_CAMERA"]
            filters = [f["FILTERS"].data[i][0] for i in range(len(f["FILTERS"].data))]
            camposList = []
            # Assume campos were specified as integers in degrees.
            for i in range(numCameras):
                camposList.append(str(int(round(degrees(f["CAMERA" + str(i) + "-PARAMETERS"].header["THETA"])))) + " " + 
                                  str(int(round(degrees(f["CAMERA" + str(i) + "-PARAMETERS"].header["PHI"])))))
            f.close()
            filterSets = []
            filterNameSets = FILTER_NAMES.split(";")
            # Convert the filter names to indices.
            for filterNames in filterNameSets:
                filterIndices = [str(filters.index(filterName + ".res")) for filterName in filterNames.split(",")]
                filterSets.append(",".join(filterIndices))
            # Write filters and campos configurations to files.
            filtersFile = folderName[:-7] + "filters"
            camposFile = folderName[:-7] + "campos"
            stdout.write("Writing " + filtersFile + " and " + camposFile + "\n")
            with open(IMAGE_DIR + filtersFile, "w") as f:
                for filterNames in filterNameSets:
                    f.write(filterNames + "\n")
            with open(IMAGE_DIR + camposFile, "w") as f:
                for campos in camposList:
                    f.write(campos + "\n")
            stdout.write("Generating RGB images for " + str(numCameras) + " perspectives with " + str(len(filters)) + " filters.\n")
            # For each specified filter set, perspective, and autopercentile, make images.
            for autoPercentile in AUTOPERCENTILES:
                for i, filterSet in enumerate(filterSets):
                    outFile = (fitsFilename[:-24] if redShift else fitsFilename[:-15])
                    outFile = (outFile[:-5] if "phys" in outFile else outFile[:-4])
                    outFile += "-filter" + str(i) + "-camera%d"
                    outFile += "-ap" + autoPercentile
                    outFile += ("-redshift.jpg" if redShift else ".jpg")
                    make_color.write_all(fitsFilename, "-BROADBAND", IMAGE_DIR + outFile, band=filterSet, 
                                         scale="autolum", autopercentile=float(autoPercentile), overwrite=True)
        chdir(pardir)
    stdout.write("\n\nFinished. Of " + str(len(listOfTimesteps)) + " time steps, " + str(len(failedTimesteps)) + " failed.\n")
    if len(failedTimesteps) > 0: stdout.write("Failed: " + str(failedTimesteps) + "\n")
    exit(0)

if __name__=="__main__":
    main()

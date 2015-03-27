#!/usr/bin/env python

################################################################################
#
# File: GenerateImages.py
# Author: Matthew Leeds <mwleeds@crimson.ua.edu>
# Last Edit: 2015-03-26
# Purpose: Generate color RGB images in JPEG format from Sunrise's broadband
# output using Patrik's make_color module. Please run this:
# hg clone https://bitbucket.org/lutorm/python
# rename that folder 'patrikpython', add an empty '__init__.py' file, and put
# it in your PYTHONPATH. Also have numpy and astropy installed.
# Config parameters should be set in 'config.ini' under 'Generate Images'.
# Broadband output files are expected to be in the format:
# SUNRISE_DIR/GALAXY_NAME-TIME_STEP-sunrise/
# GALAXY_NAME.TIME_STEP.DIAMETERkpc.phys|sim.broadband.fits
#
################################################################################

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
    galaxyName = config.get(SECTION_NAME, "GALAXY_NAME") 
    sunriseDir = config.get(SECTION_NAME, "SUNRISE_DIR")
    if len(sunriseDir) == 0: sunriseDir = curdir + sep
    listOfTimesteps = []
    timeSteps = config.get(SECTION_NAME, "TIME_STEP")
    if len(timeSteps) > 0: 
        listOfTimesteps = timeSteps.split(",")
    # If a time step wasn't specified, look for them.
    if len(listOfTimesteps) == 0:
        for f in listdir(curdir):
            if match(r"^" + galaxyName + "-\d+-sunrise$", f) != None:
                listOfTimesteps.append(f.split("-")[1])
    outFolder = config.get(SECTION_NAME, "IMAGE_DIR")
    if len(outFolder) == 0: outFolder = curdir + sep
    chdir(outFolder)
    # Make images for every perspective and filter set combination for every time step.
    for timeStep in listOfTimesteps:
        folderName = galaxyName + "-" + timeStep + "-sunrise"
        chdir(folderName)
        fitsFilename = ""
        # Find the broadband FITS file for this time step.
        for f in listdir(curdir):
            if match(r"^" + galaxyName + "\." + timeStep + "\.\d+kpc\.(phys|sim)\.broadband(-redshift)?\.fits$", f) != None:
                fitsFilename = f
                break
        if len(fitsFilename) == 0:
            stderr.write("Error: unable to find broadband FITS file for " + galaxyName + "-" + timeStep + "\n")
            continue
        stdout.write("Reading " + fitsFilename + "\n")
        f = fits.open(fitsFilename)
        numCameras = f["MCRX"].header["N_CAMERA"]
        filters = [f["FILTERS"].data[i][0] for i in range(len(f["FILTERS"].data))]
        f.close()
        filterSets = []
        filterNameSets = config.get(SECTION_NAME, "FILTER_NAMES").split(";")
        # Convert the filter names to indices.
        for filterNames in filterNameSets:
            filterIndices = [str(filters.index(filterName + ".res")) for filterName in filterNames.split(",")]
            filterSets.append(",".join(filterIndices))
        stdout.write("Generating RGB images for " + str(numCameras) + " perspectives with " + str(len(filters)) + " filters.\n")
        autoPercentile = float(config.get(SECTION_NAME, "AUTOPERCENTILE"))
        # For each specified set of filters, make images from every perspective.
        for i, filterSet in enumerate(filterSets):
            outFile = fitsFilename[:-15] + "-filter" + str(i)
            make_color.write_all(fitsFilename, "-BROADBAND", outFolder + outFile + "-camera%d.jpg", 
                                 band=filterSet, scale="auto", autopercentile=autoPercentile, overwrite=True)
        chdir(pardir)
    exit(0)

if __name__=="__main__":
    main()

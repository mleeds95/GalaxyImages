#!/usr/bin/env python

################################################################################
#
# File: GenerateImages.py
# Author: Matthew Leeds <mwleeds@crimson.ua.edu>
# Last Edit: 2015-03-25
# Purpose: Generate color RGB images in JPEG format from Sunrise's broadband
# output using Patrik's make_color module. Please run this:
# hg clone https://bitbucket.org/lutorm/python
# rename that folder 'patrikpython', add an empty '__init__.py' file, and put
# it in your PYTHONPATH. Also have numpy and astropy installed.
# Config parameters should be set in 'config.ini' under 'Generate Images'.
#
################################################################################

import os
import sys
import re
from astropy.io import fits
from patrikpython import make_color
from numpy import transpose
from ConfigParser import SafeConfigParser, NoOptionError

CONFIG_FILE = "config.ini"
SECTION_NAME = "Generate Images"

def main():
    config = SafeConfigParser()
    sys.stdout.write("Reading " + CONFIG_FILE + "\n")
    config.read(CONFIG_FILE)
    if not config.has_section(SECTION_NAME):
        sys.stderr.write("Error: No '" + SECTION_NAME + "' section found in " + CONFIG_FILE + ".\n")
        sys.exit(1)
    galaxyName = config.get(SECTION_NAME, "GALAXY_NAME") 
    os.chdir(config.get(SECTION_NAME, "SUNRISE_DIR"))
    listOfTimesteps = []
    try:
        timeStep = config.get(SECTION_NAME, "TIME_STEP")
        if len(timeStep) > 0: listOfTimesteps.append(timeStep)
    except NoOptionError:
        pass
    # If a time step wasn't specified, look for them.
    if len(listOfTimesteps) == 0:
        for f in os.listdir("."):
            if re.match(r"^" + galaxyName + "-\d+-sunrise$", f) != None:
                listOfTimesteps.append(f.split("-")[1])
    # Make images for every perspective and filter set combination for every time step.
    outFolder = config.get(SECTION_NAME, "IMAGE_DIR")
    for timeStep in listOfTimesteps:
        folderName = galaxyName + "-" + timeStep + "-sunrise"
        fitsFilename = ""
        # Find the broadband FITS file for this time step.
        for f in os.listdir(folderName):
            #TODO improve regex
            if re.match(r"^" + galaxyName + "\." + timeStep + ".*\.broadband(-redshift)?\.fits$", f) != None:
                fitsFilename = f 
                break
        if len(fitsFilename) == 0:
            sys.stderr.write("Error: unable to find broadband FITS file for " + galaxyName + "-" + timeStep + "\n")
            continue
        sys.stdout.write("Reading " + fitsFilename + "\n")
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
        sys.stdout.write("Generating RGB images for " + str(numCameras) + " perspectives with " + str(len(filters)) + " filters.\n")
        autoPercentile = float(config.get(SECTION_NAME, "AUTOPERCENTILE"))
        # For each specified set of filters, make images from every perspective.
        for i, filterSet in enumerate(filterSets):
            outFile = fitsFilename[:-15] + "-filter" + str(i)
            make_color.write_all(fitsFilename, "-BROADBAND", outFolder + outFile + "-camera%d.jpg", 
                                 band=filterSet, scale="auto", autopercentile=autoPercentile, overwrite=True)
    sys.exit(0)

if __name__=="__main__":
    main()

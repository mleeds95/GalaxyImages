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
# Image files will be written to the same directory as the input fits file.
# Usage: 'python GenerateImages.py <broadband-outfile>.fits'
#
################################################################################

import sys
from astropy.io import fits
from patrikpython import make_color
from numpy import transpose
import ConfigParser

if len(sys.argv) != 2:
    print "Usage: 'python GenerateImages.py <broadband-outfile>.fits'"
    sys.exit(1)

CONFIG_FILE = "config.ini"

def main(args):
    fitsFilename = args[0]
    config = ConfigParser.SafeConfigParser()
    sys.stdout.write("Reading " + CONFIG_FILE + "\n")
    config.read(CONFIG_FILE)
    if not config.has_section("Generate Images"):
        sys.stderr.write("Error: No \"Generate Images\" section found in " + CONFIG_FILE + ".\n")
        sys.exit(1)
    sys.stdout.write("Reading " + fitsFilename + "\n")
    f = fits.open(fitsFilename)
    numCameras = f['MCRX'].header['N_CAMERA']
    filters = [f['FILTERS'].data[i][0] for i in range(len(f['FILTERS'].data))]
    sys.stdout.write("Generating RGB images for " + str(numCameras) + " perspectives with " + str(len(filters)) + " filters.\n")
    f.close()
    autoPercentile = float(config.get('Generate Images', 'AUTOPERCENTILE'))
    filterSets = config.get('Generate Images', 'FILTER_SETS').split(';')
    for filterSet in filterSets:
        make_color.write_all(fitsFilename, '-BROADBAND', fitsFilename[:-5] + '-%d.jpg', 
                             band=filterSet, scale='auto', autopercentile=autoPercentile, overwrite=True)
    sys.exit(0)

if __name__=='__main__':
    main(sys.argv[1:])

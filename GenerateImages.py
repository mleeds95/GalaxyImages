#!/usr/bin/env python

###############################################################################
#
# File: GenerateImages.py
# Author: Matthew Leeds <mwleeds@crimson.ua.edu>
# Last Edit: 2015-03-25
# Purpose: Generate color RGB images in JPEG format from Sunrise's broadband
# output using Patrik's make_color module. Please make sure this code:
# https://bitbucket.org/lutorm/python
# is in your PATH and you have astropy and numpy installed.
# Usage: 'python GenerateImages.py <broadband-outfile>.fits'
#
###############################################################################

import sys
from astropy.io import fits
from make_color import make_image
from numpy import transpose
import ConfigParser

if len(sys.argv) != 2:
    print "Usage: 'python GenerateImages.py <broadband-outfile>.fits'"
    sys.exit(1)

CONFIG_FILE = "config.ini"
BANDS=(3,2,1) # temporarily hardcoded

def main(args):
    config = ConfigParser.SafeConfigParser()
    sys.stdout.write("Reading " + CONFIG_FILE + "\n")
    config.read(CONFIG_FILE)
    if not config.has_section("Generate Images"):
        sys.stderr.write("Error: No \"Generate Images\" section found in " + CONFIG_FILE + ".\n")
        sys.exit(1)
#TODO TODO TODO TODO remove dependence on camposFile, etc.
    broadbandFile = fits.open(args[0])
    for i, line in enumerate(camposFile):
        im = transpose(broadbandFile['CAMERA' + str(i) + '-BROADBAND'].data, axes=(1,2,0))
        img = make_image(im, band=BANDS, scale='auto', autopercentile=AUTOPERCENTILE, return_jpeg=True) 
        filename = args[0] + '_' + str(i) + '.jpg'
        with open(filename, 'w') as f:
            f.write(img)
        print(filename + ': ' + line)
    broadbandFile.close()
    sys.exit(0)

if __name__=='__main__':
    main(sys.argv[1:])

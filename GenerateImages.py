#!/usr/bin/python
# File: GenerateImages.py
# Author: Matthew Leeds
# Last Edit: 2015-03-04
# Purpose: Generate color images from sunrise's broadband output
# using Patrik's make_color module.
# Usage: ./GenerateImages broadband_out.fits campos

import sys
sys.path.append('/home/matthew/Desktop/GalaxyImages/patrik-python')
import pyfits
from make_color import make_image
from numpy import transpose

BANDS=(3,2,1) # temporarily hardcoded
AUTOPERCENTILE=0.1

def main(args):
    broadbandFile = pyfits.open(args[0])
    camposFile = open(args[1])
    for i, line in enumerate(camposFile):
        im = transpose(broadbandFile['CAMERA' + str(i) + '-BROADBAND'].data, axes=(1,2,0))
        img = make_image(im, band=BANDS, scale='auto', autopercentile=AUTOPERCENTILE, return_jpeg=True) 
        filename = args[0] + '_' + str(i) + '.jpg'
        with open(filename, 'w') as f:
            f.write(img)
        print(filename + ': ' + line)
    broadbandFile.close()
    camposFile.close()
    sys.exit(0)

if __name__=='__main__':
    main(sys.argv[1:])

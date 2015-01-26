#!/usr/bin/python

###################################################################
#
# File: ExtractGalaxy.py
# Author: Matthew Leeds
# Last Edit: 2015-01-20
# Purpose: Use PyNBody to extract the galaxy and its surroundings
# for a certain radius.
#
###################################################################

import pynbody

def main():
    pynbody.config['verbose'] = True
    print "Loading g15784 at t=1024"
    galaxy = pynbody.load("./MUGS-g15784/g15784.01024")
    print "ngas = %e, ndark = %e, nstar = %e"%(len(galaxy.g), len(galaxy.d), len(galaxy.s))
    print "Converting to physical units"
    galaxy.physical_units()
    print "Finding the largest halo"
    h1 = galaxy.halos()[1]    
    print "ngas = %e, ndark = %e, nstar = %e"%(len(h1.g), len(h1.d), len(h1.s))
    print "Centering the simulation around the main halo"
    hcenter = pynbody.analysis.halo.center(h1, mode='hyb')
    print(hcenter)
    print "Filtering out a sphere"
    sphere = h1[pynbody.filt.Sphere('50 kpc')]
    #print(sphere)
    print "Writing the sphere to the disk in Tipsy format"
    sphere.write(filename=sphere.filename, fmt=pynbody.tipsy.TipsySnap)
    
main()

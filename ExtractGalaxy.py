#!/usr/bin/python

###################################################################
#
# File: ExtractGalaxy.py
# Author: Matthew Leeds
# Last Edit: 2015-01-26
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
    pynbody.analysis.halo.center(h1, mode='hyb')
    print "Filtering out a sphere"
    sphere = galaxy[pynbody.filt.Sphere('50 kpc')]
    print "ngas = %e, ndark = %e, nstar = %e"%(len(sphere.g), len(sphere.d), len(sphere.s))

    print "Writing the star and gas particles from the sphere to the disk in Tipsy format"
    sphere.st.union(sphere.gas).write(filename=sphere.filename, fmt=pynbody.tipsy.TipsySnap)
    
if __name__=="__main__":
    main()

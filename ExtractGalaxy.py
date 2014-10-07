#!/usr/bin/python

import pynbody

def main():
    print "Loading g15784 at t=1024"
    galaxy = pynbody.load("./MUGS-g15784/g15784.01024")
    print "ngas = %e, ndark = %e, nstar = %e"%(len(galaxy.g), len(galaxy.d), len(galaxy.s))
    print "Converting to physical units"
    galaxy.physical_units()
    print "Finding the largest halo"
    h1 = galaxy.halos()[1]    
    print "ngas = %e, ndark = %e, nstar = %e"%(len(h1.g), len(h1.d), len(h1.s))
    print "Centering the simulation around the main halo"
    hcenter = pynbody.analysis.halo.center(h1, mode='com')
    print(hcenter)
    sphere = h1[pynbody.filt.Sphere('100 kpc')]
    print(sphere)

main()

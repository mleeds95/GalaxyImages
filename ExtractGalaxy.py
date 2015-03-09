#!/usr/bin/env python

########################################################################
#
# File: ExtractGalaxy.py
# Author: Matthew Leeds
# Last Edit: 2015-03-08
# Purpose: Use PyNBody to extract the galaxy and its surroundings
# for a certain radius from the larger simulation file, write it out
# in standard tipsy and ascii formats, and run SMOOTH on it. Please 
# ensure pynbody is installed, std2ascii and smooth are in your PATH, 
# and you're a directory above the simulation files.
# Usage: 
# $ python ExtractGalaxy.py [-v] <sim_name> sim|phys cube|sphere "x kpc"
#
########################################################################

import pynbody
import subprocess
import sys
import os
import math

if len(sys.argv) != 5 and len(sys.argv) != 6::
    print "Usage: python ExtractGalaxy.py [-v] <sim_name> sim|phys cube|sphere \"x kpc\""
    sys.exit(1)

# determine config parameters from command line input
VERBOSE = False
if sys.argv[1] == "-v":
    VERBOSE = True
    sys.argv.pop(1)
SNAP_NAME = sys.argv[1] # g15784.01024 for example
PHYS = (sys.argv[2] == "phys") # else simulation units
CUBE = (sys.argv[3] == "cube") # else sphere
CUT_RADIUS = sys.argv[4] # '25 kpc' for example

# determine output file names
G_NAME = SNAP_NAME.split(".")[0]
OUTFILE = SNAP_NAME + "." + str(2 * int(CUT_RADIUS.split(" ")[0]))
OUTFILE += (".phys" if PHYS else ".sim")
OUTFILESTD = OUTFILE + ".stdtipsy"
OUTFILEASC = OUTFILE + ".ascii"

def main():
    if VERBOSE:
        pynbody.config["verbose"] = True
    print "Loading galaxy simulation file " + SNAP_NAME
    sim = pynbody.load("./MUGS-" + G_NAME + "/" + SNAP_NAME)
    if VERBOSE: print "ngas = %e, ndark = %e, nstar = %e"%(len(sim.g), len(sim.d), len(sim.s))
    if PHYS:
        if VERBOSE: print "Converting to physical units"
        sim.physical_units()
    if VERBOSE: print "Finding the latest time value"
    snaptime = str(float(sim.s['tform'].max()))
    if VERBOSE: print "Accounting for redshift"
    a = sim.properties['a']
    if VERBOSE: print "Finding Unit conversion values and m_star_creation"
    if PHYS:
        UnitMass = 1.9889e33
        UnitLength = 3.08567758e21
        UnitVelocity = 1e5
        MStarCreation = 63231.4
    else:
        with open("./MUGS-" + G_NAME + "/" + G_NAME + ".param") as paramfile:
            for line in paramfile:
                if line.startswith("dMsolUnit"):
                    UnitMass = 1.98892e33 * float(line.split("=")[1].strip())
                if line.startswith("dKpcUnit"):
                    UnitLength = (3.086e21 * float(line.split("=")[1].strip())) / a
                if line.startswith("dInitStarMass"):
                    MStarCreation = float(line.split("=")[1].strip()) 
        UnitVelocity = 1.72756e8 * math.sqrt(a)
    # TODO generate config file
    if VERBOSE: print "Finding the largest halo"
    h1 = sim.halos()[1]    
    if VERBOSE: print "ngas = %e, ndark = %e, nstar = %e"%(len(h1.g), len(h1.d), len(h1.s))
    if VERBOSE: print "Centering the simulation around the main halo"
    pynbody.analysis.halo.center(h1, mode="hyb")
    if VERBOSE: print "Rotating to a face-on view"
    pynbody.analysis.angmom.faceon(h1, cen=(0,0,0))
    print "Filtering out a " + CUT_RADIUS + " radius"
    if CUBE:
        cut = h1[pynbody.filt.Cuboid("-" + CUT_RADIUS)]
    else:
        cut = h1[pynbody.filt.Sphere(CUT_RADIUS)]
    if VERBOSE: print "ngas = %e, ndark = %e, nstar = %e"%(len(cut.g), len(cut.d), len(cut.s))
    FNULL = open(os.devnull, 'w')
    # write out the snapshot in standard tipsy format
    cut.write(filename=OUTFILESTD, fmt=pynbody.tipsy.TipsySnap)
    print "Writing " + OUTFILEASC
    cmd = "cat " + OUTFILESTD + " | std2ascii > " + OUTFILEASC
    p1 = subprocess.Popen(cmd, shell=True, stdout=FNULL, stderr=FNULL)
    p1.wait()
    if p1.returncode != 0:
        print "Error running std2ascii. Perhaps it's not in your PATH?"
        FNULL.close()
        sys.exit(1)
    print "Writing smoothing lengths to smooth.hsm"
    cmd = "smooth hsmooth < " + OUTFILESTD
    p2 = subprocess.Popen(cmd, shell=True, stdout=FNULL, stderr=FNULL)
    p2.wait()
    if p2.returncode != 0:
        print "Error running smooth. Perhaps it's not in your PATH?"
        FNULL.close()
        sys.exit(1)
    # Sunrise just needs the ASCII version
    cmd = "rm " + OUTFILESTD
    subprocess.Popen(cmd, shell=True, stdout=FNULL, stderr=FNULL).wait()
    FNULL.close()
    sys.exit(0)
    
if __name__=="__main__":
    main()

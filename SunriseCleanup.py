#!/usr/bin/env python

###############################################################################
#
# File: SunriseCleanup.py
# Author: Matthew Leeds <mwleeds@crimson.ua.edu>
# License: GNU GPL v3 <gnu.org/licenses>
# Last Edit: 2015-04-05
# Purpose: Assuming Sunrise has finished running for the specified time steps,
# delete the intermediate files and tar up the useful ones (broadband output
# and sfrhist input). Without doing this, each galaxy would use roughly 
# a Terabyte of disk space.
#
###############################################################################

import os
from sys import stdout, stderr, exit
from re import match
from ConfigParser import SafeConfigParser
from subprocess import Popen

CONFIG_FILE = "config.ini"
SECTION_NAME = "Sunrise Cleanup"

def main():
    config = SafeConfigParser()
    stdout.write("Reading " + CONFIG_FILE + "\n")
    config.read(CONFIG_FILE)
    if not config.has_section(SECTION_NAME):
        stderr.write("Error: No '" + SECTION_NAME + "' section found in " + CONFIG_FILE + ".\n")
        exit(1)
    GALAXY_NAME = config.get(SECTION_NAME, "GALAXY_NAME") 
    RUN_DIR = config.get(SECTION_NAME, "RUN_DIR")
    if len(RUN_DIR) == 0: RUN_DIR = os.curdir + os.sep
    os.chdir(RUN_DIR)
    listOfTimesteps = []
    TIME_STEPS = config.get(SECTION_NAME, "TIME_STEPS")
    if len(TIME_STEPS) > 0: 
        listOfTimesteps = TIME_STEPS.split(",")
    # If a time step wasn't specified, look for them.
    if len(listOfTimesteps) == 0:
        for f in os.listdir(os.curdir):
            if match(r"^" + GALAXY_NAME + "-\d+-sunrise$", f) != None:
                listOfTimesteps.append(f.split("-")[1])
    OUT_DIR = config.get(SECTION_NAME, "OUT_DIR")
    if len(OUT_DIR) == 0: OUT_DIR = os.curdir + os.sep
    listOfFolderNames = [] # use later for making the tarball
    # For each time step, delete unnecessary files (unless the run failed).
    for timeStep in listOfTimesteps:
        folderName = GALAXY_NAME + "-" + timeStep + "-sunrise"
        listOfFolderNames.append(folderName)
        os.chdir(folderName)
        # Try to find the broadband FITS file for this time step.
        foundBroadband = False
        broadbandFitsRE = r"^" + GALAXY_NAME + "\." + timeStep + "\.\d+kpc\.(phys|sim)\.broadband(-redshift)?\.fits$"
        for f in os.listdir(os.curdir):
            if match(broadbandFitsRE, f) != None:
                foundBroadband = True
        # Look at the broadband output files to see if it succeeded.
        broadbandSuccess = False
        with open(os.devnull, "w") as FNULL:
            for fname in ("broadband.out", "broadband-redshift.out"):
                if os.isfile(fname):
                    cmd = "head -n30 " + fname + " | grep \"Successfully completed.\""
                    if (Popen(cmd, shell=True, stdout=FNULL, stderr=FNULL).wait()) == 0:
                        broadbandSuccess = True
        if not foundBroadband or not broadbandSuccess:
            print GALAXY_NAME + " appears not to have finished successfully. Leaving its files alone." 
            continue
        # If Broadband succeeded delete everything but broadband output and sfrhist input 
        sfrhistSnapRE = r"^" + GALAXY_NAME + "\." + timeStep + "\.\d+kpc\.(phys|sim)\.ascii"
        for f in os.listdir(os.curdir):
            if match(broadbandFitsRE, f) == None and match(sfrhistSnapRE, f) == None:
                os.remove(f)
        os.chdir(os.pardir) # go back up to RUN_DIR
    # End for loop over time steps.
    # Make a gzipped tarball with all the slimmed down output folders.
    if config.getboolean(SECTION_NAME, "TARBALL"):
        listOfTimestepInts = [int(x) for x in listOfTimesteps]
        minTimestep = min(listOfTimestepInts)
        maxTimestep = max(listOfTimestepInts)
        tarball = "sunrise-output-" + GALAXY_NAME + "-t" + format(minTimestep, "05")
        if len(listOfTimesteps) > 1: tarball += "-" + format(maxTimestep, "05")
        tarball += ".tgz"
        cmd = "tar czf " + OUT_DIR + tarball + " " + " ".join(listOfFolderNames) + " core.*"
        p = Popen(cmd, shell=True).wait()
        if p != 0:
            print "Erroring creating tarball " + tarball + ". (" + str(p) + ")"
            exit(1)
    print "Successfully checked " + len(listOfTimesteps) + " directories."
    exit(0)

if __name__=="__main__":
    main()

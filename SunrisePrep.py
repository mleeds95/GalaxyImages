#!/usr/bin/env python

###############################################################################
#
# File: SunrisePrep.py
# Author: Matthew Leeds <mwleeds@crimson.ua.edu>
# Contributor: Owain Snaith
# License: GNU GPL v3 <gnu.org/licenses>
# Last Edit: 2015-04-07
# Purpose: Determine the trunk of the merger tree, extract the galaxy 
# for a certain radius from the larger simulation file, run SMOOTH on it,
# generate the Sunrise config files, and write job submission scripts. Please 
# ensure pynbody is installed, std2ascii and smooth are in your PATH, and all 
# necessary parameters are set. Campos, filters, sfrhist/mcrx stub files, and 
# 'config.ini' should be in the current directory. Note that the format 
# <galaxy name>.<time step> for simulation files is hard-coded, as is the 
# format for param files. You probably want to use the "screen" utility to
# run this so it will continue if your session disconnects, and it's probably
# also a good idea to write the output to a log file.
# Usage: $ python SunrisePrep.py >> sunrise-prep.log 2>&1
#
###############################################################################

import pynbody
import subprocess
import sys
import os
import shutil
import math
import ConfigParser
import re
import matplotlib
import scipy
import pickle
import numpy
import datetime

CONFIG_FILE = "config.ini"
SECTION_NAME = "Sunrise Prep"
WORKING_DIR = os.getcwd() + os.sep

def main():
    #
    #Step 0: Read the config file and initialize listOfTimesteps
    #
    sys.stdout.write("======== Running Sunrise Prep at " + str(datetime.datetime.now()) + " ========\n")
    config = ConfigParser.SafeConfigParser()
    sys.stdout.write("Reading " + CONFIG_FILE + "\n")
    config.read(CONFIG_FILE)
    if not config.has_section(SECTION_NAME):
        sys.stderr.write("Error: No \"" + SECTION_NAME + "\" section found in " + CONFIG_FILE + ".\n")
        sys.exit(1)
    VERBOSE = config.getboolean(SECTION_NAME, "VERBOSE")
    if VERBOSE: pynbody.config["verbose"] = True
    SIM_DIR = config.get(SECTION_NAME, "SIM_DIR")
    RUN_DIR = config.get(SECTION_NAME, "RUN_DIR")
    OUT_DIR = config.get(SECTION_NAME, "OUT_DIR")
    SHARE_DIR = config.get(SECTION_NAME, "SHARE_DIR")
    SRC_DIR = config.get(SECTION_NAME , "SRC_DIR")
    BIN_DIR = config.get(SECTION_NAME, "BIN_DIR")
    # Enforce the use of absulute pathing (assuming a Unix-like file system).
    for directory in (SIM_DIR, RUN_DIR, OUT_DIR, SHARE_DIR, SRC_DIR, BIN_DIR):
        if directory[0] != os.sep or directory[-1] != os.sep:
            sys.stderr.write("Error: incorrectly formatted directory config parameter. It should start and end with a '" + os.sep + "'\n")
            sys.exit(1)
    # Make sure the local directories exist.
    if not os.path.isdir(SIM_DIR):
        sys.stderr.write("Error: simulation directory '" + SIM_DIR + "' does not exist.\n")
        sys.exit(1)
    if not os.path.isdir(OUT_DIR):
        sys.stdout.write("Attempting to create directory '" + OUT_DIR + "'.\n")
        os.makedirs(OUT_DIR)
    GALAXY_NAME = config.get(SECTION_NAME, "GALAXY_NAME")
    # Find every time step for the galaxy (for halo tracking).
    allTimesteps = []
    for f in os.listdir(SIM_DIR):
        if re.match(r"^" + GALAXY_NAME + "\.\d+$", f) != None:
            allTimesteps.append(f.split(".")[1])
    if len(allTimesteps) == 0:
        sys.stderr.write("Error: No time steps specified or found in " + SIM_DIR + "\n")
        sys.exit(1)
    # Check if any time steps were specified.
    listOfTimesteps = []
    try:
        TIME_STEPS = config.get("Sunrise Prep", "TIME_STEPS")
        if len(TIME_STEPS) > 0:
            listOfTimesteps = TIME_STEPS.split(",")
    except ConfigParser.NoOptionError:
        pass
    # If no time steps were specified, process them all.
    if len(listOfTimesteps) == 0:
        listOfTimesteps = allTimesteps[:]
    # Read the rest of the parameters.
    GALAXY_SET = config.get(SECTION_NAME, "GALAXY_SET")
    PARAM_FILE = config.get(SECTION_NAME, "PARAM_FILE")
    PHYS = config.getboolean(SECTION_NAME, "PHYS")
    CUBE = config.getboolean(SECTION_NAME, "CUBE")
    CUT_RADIUS = config.get(SECTION_NAME, "CUT_RADIUS")
    cutDiameter = str(int(CUT_RADIUS.split(" ")[0]) * 2)
    CAMPOS_FILE = config.get(SECTION_NAME, "CAMPOS_FILE")
    FILTERS_FILE = config.get(SECTION_NAME, "FILTERS_FILE")
    QUEUE_NAME = config.get(SECTION_NAME, "QUEUE_NAME")
    MAX_LEVEL = config.get(SECTION_NAME, "MAX_LEVEL")
    N_THREADS = config.get(SECTION_NAME, "N_THREADS")
    SFRHIST_STUB = config.get(SECTION_NAME, "SFRHIST_STUB")
    MCRX_STUB = config.get(SECTION_NAME, "MCRX_STUB")
    AUTO_RUN = config.getboolean(SECTION_NAME, "AUTO_RUN")
    TARBALL = config.getboolean(SECTION_NAME, "TARBALL")
    #
    #Step 1: Traverse the time steps to find the trunk of the merger tree.
    #
    sys.stdout.write("Traversing all time steps for this galaxy to find the trunk of the merger tree.\n")
    numTimesteps = len(allTimesteps)
    # Ensure that we traverse the time steps in reverse chronological order.
    allTimestepInts = [int(x) for x in allTimesteps]
    allTimestepInts.sort(reverse=True)
    allTimesteps = [format(x, "05") for x in allTimestepInts]
    # Initialize arrays to hold the star halo IDs.
    starIDs = numpy.full(numTimesteps, -999, dtype=numpy.int32)
    starHalo = 1
    # Traverse the time steps, comparing each to its predecessor.
    os.chdir(SIM_DIR)
    starIDsDict = {}
    for i in range(numTimesteps - 1):
        t1 = startTimer()
        if VERBOSE: sys.stdout.write("Comparing " + allTimesteps[i] + " and " + allTimesteps[i+1] + "\n")
        # Load the ith snapshot.
        snap1 = pynbody.load(GALAXY_NAME + "." + allTimesteps[i])
        snap1.physical_units()
        idx = numpy.zeros(shape=(2, len(snap1.d)))
        halo1 = snap1.halos()
        halo1.make_grp()
        if len(halo1) == 0:
            if VERBOSE: sys.stdout.write(GALAXY_NAME + "." + allTimesteps[i] + " has no halos.\n")
            continue
        # Use a modified version of pynbody's center, which doesn't always do what you expect.
        mycenter(halo1[1], 1, mode="myhyb")
        pynbody.analysis.angmom.faceon(halo1[1], cen=[0,0,0])
        # Load the (i+1)th snapshot. (really the previous one)
        snap2 = pynbody.load(GALAXY_NAME + "." + allTimesteps[i+1])
        snap2.physical_units()
        halo2 = snap2.halos()
        halo2.make_grp()
        if len(halo2) == 0:
            if VERBOSE: sys.stdout.write(GALAXY_NAME + "." + allTimesteps[i+1] + " has no halos.\n")
            continue
        mycenter(halo2[1], 1, mode='myhyb')
        pynbody.analysis.angmom.faceon(halo2[1], cen=[0,0,0])
        # Use the stars to walk across the snapshots.
        ind = snap1.s['grp'] == starHalo 
        nstar = len(snap2.s)
        ind2 = ind[:nstar]
        ids = snap2.s['grp'][ind2]
        uids = numpy.unique(ids)
        if(numpy.size(uids) > 0):
            muids = numpy.concatenate((uids, [uids.max()+1]))        
            hh, hids = numpy.histogram(ids, bins=muids)
            ind = hids[:-1] != 0
            if(hh[ind].sum() > 0):
                imax = hh[ind].argmax()
                starHalo = hids[ind][imax]
                starIDs[i] = starHalo
                starIDsDict[allTimesteps[i]] = starHalo
        if VERBOSE: stopTimer(t1)
    # End loop traversing all the time steps. starIDs should now be correct.
    # Save the starIDs for documentation.
    starIDsFile = "mainbranch_" + GALAXY_NAME + "_starIDs.txt"
    numpy.savetxt(OUT_DIR + starIDsFile, starIDs, fmt="%d")
    # Iterate over all the specified time steps generating appropriate files.
    # Notice that we iterate over a copy of the list so we can remove timesteps with no stars.
    for timeStep in listOfTimesteps[:]:
        #
        #Step 2: Load the simulation snapshot and calculate some values for the Sunrise config files.
        #
        t1 = startTimer()
        os.chdir(SIM_DIR)
        snapName = GALAXY_NAME + "." + timeStep
        sys.stdout.write("Loading galaxy simulation file " + SIM_DIR + snapName + "\n")
        sim = pynbody.load(snapName)
        if VERBOSE: sys.stdout.write("ngas = %e, ndark = %e, nstar = %e\n"%(len(sim.g), len(sim.d), len(sim.s)))
        # If there aren't any stars yet, skip the timestep.
        if len(sim.s) == 0:
            sys.stdout.write("Skipping " + snapName + " since there aren't any stars in it.\n")
            # remove the timestep so we don't include it in the controller script later.
            listOfTimesteps.remove(timeStep)
            continue
        if PHYS:
            if VERBOSE: sys.stdout.write("Converting to physical units\n")
            sim.physical_units()
        if VERBOSE: sys.stdout.write("Finding the latest time value\n")
        snaptime = float(sim.s["tform"].max())
        # grab the values of a (1 / 1+z) and z so we can account for redshift
        a = sim.properties["a"]
        z = sim.properties["z"]
        nonzeroRedshift = (abs(z) > 1e-8)
        if VERBOSE: sys.stdout.write("Finding Unit conversion values and m_star_creation\n")
        if PHYS:
            UnitMass = 1.9889e33
            UnitLength = 3.08567758e21
            UnitVelocity = 1e5
            MStarCreation = 63231.4
        else:
            with open(PARAM_FILE) as paramfile:
                for line in paramfile:
                    if line.startswith("dMsolUnit"):
                        UnitMass = 1.98892e33 * float(line.split("=")[1].strip())
                    if line.startswith("dKpcUnit"):
                        UnitLength = (3.086e21 * float(line.split("=")[1].strip())) / a
                    if line.startswith("dInitStarMass"):
                        MStarCreation = float(line.split("=")[1].strip()) 
            UnitVelocity = 1.72756e8 * math.sqrt(a)
        if VERBOSE: stopTimer(t1)
        #
        #Step 3: Cut out the specified radius.
        #
        t1 = startTimer()
        if VERBOSE: sys.stdout.write("Choosing the appropriate halo.\n")
        try:
            h1 = sim.halos()[starIDsDict[timeStep]]
        except KeyError:
            sys.stderr.write("Warning: Defaulting to halo 1 for " + snapName + ".\n")
            h1 = sim.halos()[1]
        if VERBOSE: sys.stdout.write("ngas = %e, ndark = %e, nstar = %e\n"%(len(h1.g), len(h1.d), len(h1.s)))
        if VERBOSE: sys.stdout.write("Centering the simulation around the main halo\n")
        pynbody.analysis.halo.center(h1, mode="hyb")
        if VERBOSE: sys.stdout.write("Rotating to a face-on view\n")
        pynbody.analysis.angmom.faceon(h1, cen=(0,0,0))
        sys.stdout.write("Filtering out a " + CUT_RADIUS + " radius\n")
        if CUBE:
            cut = h1[pynbody.filt.Cuboid("-" + CUT_RADIUS)]
        else:
            cut = h1[pynbody.filt.Sphere(CUT_RADIUS)]
        if VERBOSE: sys.stdout.write("ngas = %e, ndark = %e, nstar = %e\n"%(len(cut.g), len(cut.d), len(cut.s)))
        if VERBOSE: stopTimer(t1)
        #
        #Step 4: Write the snapshot section to the disk in std tipsy and ASCII formats.
        #
        t1 = startTimer()
        os.chdir(WORKING_DIR)
        # <galaxy name>.<time step>.<diameter>kpc.phys|sim.stdtipsy|ascii
        SNAPFILE = GALAXY_NAME + "." + timeStep + "." + cutDiameter + "kpc"
        SNAPFILE += (".phys" if PHYS else ".sim")
        SNAPFILESTD = SNAPFILE + ".stdtipsy"
        SNAPFILEASC = SNAPFILE + ".ascii"
        sys.stdout.write("Writing " + SNAPFILESTD + "\n")
        cut.write(filename=SNAPFILESTD, fmt=pynbody.tipsy.TipsySnap)
        sys.stdout.write("Writing " + SNAPFILEASC + "\n")
        cmd = "cat " + SNAPFILESTD + " | std2ascii > " + SNAPFILEASC
        with open(os.devnull, "w") as FNULL:
            p1 = subprocess.Popen(cmd, shell=True, stdout=FNULL, stderr=FNULL)
            p1.wait()
            if p1.returncode != os.EX_OK:
                sys.stderr.write("Error running std2ascii. Perhaps it's not in your PATH?\n")
                sys.exit(1)
        if VERBOSE: stopTimer(t1)
        #
        #Step 5: Run SMOOTH.
        #
        t1 = startTimer()
        sys.stdout.write("Writing smoothing lengths to smooth.hsm\n")
        cmd = "smooth hsmooth < " + SNAPFILESTD
        with open(os.devnull, "w") as FNULL:
            p2 = subprocess.Popen(cmd, shell=True, stdout=FNULL, stderr=FNULL)
            p2.wait()
            if p2.returncode != os.EX_OK:
                sys.stderr.write("Error running smooth. Perhaps it's not in your PATH?\n")
                sys.exit(1)
        os.remove(SNAPFILESTD) # Sunrise just needs the ASCII version
        if VERBOSE: stopTimer(t1)
        #
        #Step 6: Generate Sunrise config files.
        #
        t1 = startTimer()
        with open(SIM_DIR + PARAM_FILE) as f:
            originalParams = f.readlines()
        # runDirName is the name of the directory in OUT_DIR, and the name of the one in RUN_DIR
        runDirName = GALAXY_NAME + "-" + timeStep + "-sunrise"
        if VERBOSE: sys.stdout.write("The running directory will be " + runDirName + "\n")
        os.chdir(OUT_DIR)
        if os.path.isdir(runDirName):
            shutil.rmtree(runDirName) # in case this isn't the first run
        os.mkdir(runDirName)
        os.chdir(runDirName)
        fullRunDir = RUN_DIR + runDirName + "/"
        sys.stdout.write("Writing config and stub files for sfrhist, mcrx, and broadband.\n")
        # In theory these values should be in a config file not the param, 
        # but Sunrise seems to like it this way.
        with open(PARAM_FILE, "w") as f:
            f.write("snaptime " + str(snaptime) + "\n")
            f.write("UnitMass_in_g " + str(UnitMass) + "\n")
            f.write("UnitLength_in_cm " + str(UnitLength) + "\n")
            f.write("UnitVelocity_in_cm_per_s " + str(UnitVelocity) + "\n")
            f.write("m_star_creation " + str(MStarCreation) + "\n")
            f.write("smoothing_file " + fullRunDir + "smooth.hsm\n")
            f.writelines(originalParams)
        # Write sfrhist config files.
        with open("sfrhist-" + snapName + ".config", "w") as f:
            f.write("include_file " + fullRunDir + "sfrhist.stub\n")
            f.write("snapshot_file " + fullRunDir + SNAPFILEASC + "\n")
            f.write("output_file " + fullRunDir + SNAPFILEASC[:-6] + ".sfrhist.fits\n")
            f.write("simparfile " + fullRunDir + PARAM_FILE + "\n")
        with open("sfrhist.stub", "w") as f:
            f.write("runname " + GALAXY_SET + "-" + GALAXY_NAME + "-" + timeStep + "\n")
            f.write("nbodycod pkdgrav\n")
            f.write("max_level " + MAX_LEVEL + "\n")
            f.write("n_threads " + N_THREADS + "\n")
            f.write("grid_min " + ("-" + CUT_RADIUS.split(" ")[0] + ". ") * 3 + "\n")
            f.write("grid_max " + (CUT_RADIUS.split(" ")[0] + ". ") * 3 + "\n")
            f.write("stellarmodelfile " + SHARE_DIR + "Patrik-imfKroupa-Zmulti-ml.fits\n")
            f.write("mappings_sed_file " + SHARE_DIR + "Smodel.fits\n")
            f.write("mcrx_data_dir " + SRC_DIR + "sunrise/src/\n")
            if not PHYS: 
                f.write("comoving_units true\n")
                f.write("expansion_factor " + str(a) + "\n")
            f.write("use_counters false\n")
            with open(WORKING_DIR + SFRHIST_STUB) as f2:
                f.writelines(f2.readlines())
        # Write mcrx config files.
        with open("mcrx-" + snapName + ".config", "w") as f:
            f.write("include_file " + fullRunDir + "mcrx.stub\n")
            f.write("input_file " + fullRunDir + SNAPFILEASC[:-6] + ".sfrhist.fits\n")
            f.write("output_file " + fullRunDir + SNAPFILEASC[:-6] + ".mcrx.fits\n")
        with open("mcrx.stub", "w") as f:
            f.write("n_threads " + N_THREADS + "\n")
            f.write("mcrx_data_dir " + SRC_DIR + "sunrise/src/\n")
            f.write("grain_data_directory " + SRC_DIR + "crosssections\n")
            f.write("camera_positions " + fullRunDir + CAMPOS_FILE + "\n")
            f.write("camerafov " + str(int(cutDiameter)) + "\n")
            f.write("use_counters false\n")
            with open(WORKING_DIR + MCRX_STUB) as f2:
                f.writelines(f2.readlines())
        # Write broadband config files.
        with open("broadband-" + snapName + ".config", "w") as f:
            f.write("include_file " + fullRunDir + "broadband.stub\n")
            f.write("input_file " + fullRunDir + SNAPFILEASC[:-6] + ".mcrx.fits\n")
            f.write("output_file " + fullRunDir + SNAPFILEASC[:-6] + ".broadband.fits\n")
        # Also make redshifted images if necessary
        if nonzeroRedshift:
            with open("broadband-" + snapName + "-redshift.config", "w") as f:
                f.write("include_file " + fullRunDir + "broadband.stub\n")
                f.write("input_file " + fullRunDir + SNAPFILEASC[:-6] + ".mcrx.fits\n")
                f.write("output_file " + fullRunDir + SNAPFILEASC[:-6] + ".broadband-redshift.fits\n")
                f.write("redshift " + str(z) + "\n")
        with open("broadband.stub", "w") as f:
            f.write("filter_file_directory " + SRC_DIR + "filters/\n")
            f.write("filter_list " + fullRunDir + FILTERS_FILE + "\n")
            f.write("filter_lambda_conversion 1e-10\n")
            f.write("use_counters false\n")
        if VERBOSE: stopTimer(t1)
        #
        #Step 7: Move the files into the final directory, and write out job submission commands.
        #
        t1 = startTimer()
        sys.stdout.write("Moving files from " + WORKING_DIR + " to " + OUT_DIR + runDirName + ".\n")
        for fileName in (SNAPFILEASC, "smooth.hsm"):
            shutil.move(WORKING_DIR + fileName, fileName)
        shutil.copy(WORKING_DIR + CAMPOS_FILE, CAMPOS_FILE)
        shutil.copy(WORKING_DIR + FILTERS_FILE, FILTERS_FILE)
        sys.stdout.write("Writing bsub commands to runsfrhist.sh, runmcrx.sh, and runbroadband.sh.\n")
        with open("runsfrhist.sh", "w") as f:
            f.write("#!/bin/bash\n\n")
            f.write("rm -f " + fullRunDir + "sfrhist.out " + fullRunDir + "sfrhist.err " + fullRunDir + SNAPFILEASC[:-6] + ".sfrhist.fits\n")
            f.write("bsub -q " + QUEUE_NAME + " -n " + N_THREADS + " -R \"span[hosts=1]\" -o " + fullRunDir + "sfrhist.out -e " + fullRunDir + "sfrhist.err")
            if AUTO_RUN: f.write(" -Ep \"bash " + fullRunDir + "runmcrx.sh\"") 
            f.write(" " + BIN_DIR + "sfrhist " + fullRunDir + "sfrhist-" + snapName + ".config\n") 
        os.chmod("runsfrhist.sh", 0744)
        with open("runmcrx.sh", "w") as f:
            f.write("#!/bin/bash\n\n")
            f.write("rm -f " + fullRunDir + "mcrx.out " + fullRunDir + "mcrx.err " + fullRunDir + SNAPFILEASC[:-6] + ".mcrx.fits\n")
            f.write("bsub -q " + QUEUE_NAME + " -n " + N_THREADS + " -R \"span[hosts=1]\" -o " + fullRunDir + "mcrx.out -e " + fullRunDir + "mcrx.err")
            if AUTO_RUN: f.write(" -Ep \"bash " + fullRunDir + "runbroadband.sh\"") 
            f.write(" " + BIN_DIR + "mcrx " + fullRunDir + "mcrx-" + snapName + ".config\n") 
        os.chmod("runmcrx.sh", 0744)
        with open("runbroadband.sh", "w") as f:
            f.write("#!/bin/bash\n\n")
            f.write("rm -f " + fullRunDir + "broadband.out " + fullRunDir + "broadband.err " + fullRunDir + SNAPFILEASC[:-6] + ".broadband.fits\n")
            f.write("bsub -q " + QUEUE_NAME + " -n " + N_THREADS + " -R \"span[hosts=1]\" -o " + fullRunDir + "broadband.out -e " + fullRunDir + "broadband.err " + BIN_DIR + "broadband " + fullRunDir + "broadband-" + snapName + ".config\n") 
            if nonzeroRedshift:
                f.write("rm -f " + fullRunDir + "broadband-redshift.out " + fullRunDir + "broadband-redshift.err " + fullRunDir + SNAPFILEASC[:-6] + ".broadband-redshift.fits\n")
                f.write("bsub -q " + QUEUE_NAME + " -n " + N_THREADS + " -R \"span[hosts=1]\" -o " + fullRunDir + "broadband-redshift.out -e " + fullRunDir + "broadband-redshift.err " + BIN_DIR + "broadband " + fullRunDir + "broadband-" + snapName + "-redshift.config\n")
        os.chmod("runbroadband.sh", 0744)
        if VERBOSE: stopTimer(t1)
    # end primary loop over time steps
    #
    #Step 8: Write a script to start sfrhist for all time steps.
    #
    t1 = startTimer()
    os.chdir(OUT_DIR)
    listOfTimestepInts = [int(x) for x in listOfTimesteps]
    minTimestep = min(listOfTimestepInts)
    maxTimestep = max(listOfTimestepInts)
    # runall-<galaxy name>-t<min timestep>[-<max timestep>].sh
    jobRunner = "runall-" + GALAXY_NAME + "-t" + format(minTimestep, '05')
    if len(listOfTimesteps) > 1: jobRunner += "-" + format(maxTimestep, '05')
    jobRunner += ".sh"
    sys.stdout.write("Writing a script to start jobs for each timestep: " + jobRunner + "\n")
    with open(jobRunner, "w") as f:
        f.write("#!/bin/bash\n\n")
        for timeStep in listOfTimesteps:
            f.write(RUN_DIR + GALAXY_NAME + "-" + timeStep + "-sunrise/runsfrhist.sh\n")
    os.chmod(jobRunner, 0744)
    if VERBOSE: stopTimer(t1)
    #
    #Step 9: If requested, make a tarball of the relevant parts of the output directory.
    #
    t1 = startTimer()
    # All the folders for this galaxy are included, even if they're not the specified time steps.
    # If you want it to only include the right ones, use a different folder for each run.
    if TARBALL:
        sys.stdout.write("Tarring up the output directory.\n")
        # <galaxy name>-t<min timestep>[-<max timestep>].tgz
        cmd = "tar czf " + jobRunner[7:-3] + ".tgz *" + GALAXY_NAME + "*"
        p3 = subprocess.Popen(cmd, shell=True)
        p3.wait()
        if p3.returncode != os.EX_OK:
            sys.stderr.write("Error encountered while tarring output directory!\n")
            sys.exit(1)
    if VERBOSE: stopTimer(t1)
    # Exit.
    sys.stdout.write("Finished. You're ready to run Sunrise!\n\n")
    sys.exit(0)

def startTimer():
    sys.stdout.write("Starting timer.\n")
    return datetime.datetime.now()

def stopTimer(t1):
    sys.stdout.write("Stopping timer. ")
    t2 = datetime.datetime.now()
    sys.stdout.write("Elapsed time: " + str(t2 - t1) + "\n")

# The following are Owain Snaith's modified versions of pynbody functions.
    
def mycenterpot(sim,haloid):
    i = sim["phi"][sim['grp']==haloid].argmin()
    return sim["pos"][sim['grp']==haloid][i].copy()
    
def my_hybrid_center(sim,haloid, r='3 kpc', **kwargs):
    """

    Determine the center of the halo by finding the shrink-sphere
    -center inside the specified distance of the potential minimum

    """

    try:
        cen_a = mycenterpot(sim,haloid)
    except KeyError:
        cen_a = pynbody.analysis.halo.center_of_mass(sim)
    return pynbody.analysis.halo.shrink_sphere_center(sim[pynbody.filt.Sphere(r, cen_a)], **kwargs)

def mycenter(sim,haloid=1, mode=None, retcen=False, vel=True, cen_size="1 kpc", move_all=True, wrap=False, **kwargs):
    """

    Determine the center of mass of the given particles using the
    specified mode, then recenter the particles (of the entire
    ancestor snapshot) accordingly

    Accepted values for *mode* are

      *mypot*: potential minimum including only host particles

      *myhyb*: for sane halos, returns the same as ssc, but works faster by
             starting iteration near potential minimum

    or a function returning the COM.

    **Other keywords:**

    *retcen*: if True only return the center without centering the
     snapshot (default = False)


    *vel*: if True, translate velocities so that the velocity of the
    central 1kpc (default) is zeroed. Other values can be passed with cen_size.

    *move_all*: if True (default), move the entire snapshot. Otherwise only move
    the particles in the halo passed in.

    *wrap*: if True, pre-centre and wrap the simulation so that halos on the edge
    of the box are handled correctly. Default False.
    """

    if mode is None:
        mode = 'myhyb'

    try:
        fn = {'mypot': mycenterpot,
              'myhyb': my_hybrid_center}[mode]
    except KeyError:
        fn = mode

    if move_all:
        target = sim.ancestor
    else:
        target = sim

    if wrap:
        # centre on something within the halo and wrap
        target = pynbody.transformation.inverse_translate(target, sim['pos'][0])
        target.sim.wrap()

    #print fn(sim,haloid, **kwargs)
    if retcen:
        return fn(sim,haloid, **kwargs)
    else:
        cen = fn(sim,haloid, **kwargs)
        tx = pynbody.transformation.inverse_translate(target, cen)

    if vel:
        velc = pynbody.analysis.halo.vel_center(sim, cen_size=cen_size, retcen=True)
        tx = pynbody.transformation.inverse_v_translate(tx, velc)

    return tx

if __name__=="__main__":
    main()

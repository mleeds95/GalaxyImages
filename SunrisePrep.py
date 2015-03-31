#!/usr/bin/env python

###############################################################################
#
# File: SunrisePrep.py
# Author: Matthew Leeds <mwleeds@crimson.ua.edu>
# License: GNU GPL v3 <gnu.org/licenses>
# Last Edit: 2015-03-30
# Purpose: Use PyNBody to extract the galaxy and its surroundings
# for a certain radius from the larger simulation file, run SMOOTH on it,
# generate the Sunrise config files, and write job submission scripts. Please 
# ensure pynbody is installed, std2ascii and smooth are in your PATH, and all 
# necessary parameters are set. Campos, filters, sfrhist/mcrx stub files, and 
# 'config.ini' should be in the current directory. Note that the format 
# <galaxy name>.<time step> for simulation files is hard-coded, as is the 
# format for param files.
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

CONFIG_FILE = "config.ini"
SECTION_NAME = "Sunrise Prep"
WORKING_DIR = os.getcwd() + os.sep

def main():
    #
    #Step 1: Read the config file and initialize listOfTimesteps
    #
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
    for directory in (SIM_DIR, OUT_DIR):
        if not os.path.isdir(directory):
            sys.stderr.write("Error: configured directory does not exist.")
            sys.exit(1)
    GALAXY_NAME = config.get(SECTION_NAME, "GALAXY_NAME")
    listOfTimesteps = []
    try:
        TIME_STEPS = config.get("Sunrise Prep", "TIME_STEPS")
        if len(TIME_STEPS) > 0:
            listOfTimesteps = TIME_STEPS.split(",")
    except ConfigParser.NoOptionError:
        pass
    # If no time steps were specified, look for them.
    if len(listOfTimesteps) == 0:
        for f in os.listdir(SIM_DIR):
            if re.match(r"^" + GALAXY_NAME + "\.\d+$", f) != None:
                listOfTimesteps.append(f.split(".")[1])
    if len(listOfTimesteps) == 0: # we couldn't find any
        sys.stderr.write("Error: No time steps specified or found in " + SIM_DIR + "\n")
        sys.exit(1)
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
    # Iterate over all the time steps generating appropriate files.
    for timeStep in listOfTimesteps:
        #
        #Step 2: Load the simulation snapshot and calculate some values for the Sunrise config files.
        #
        os.chdir(SIM_DIR)
        snapName = GALAXY_NAME + "." + timeStep
        sys.stdout.write("Loading galaxy simulation file " + SIM_DIR + snapName + "\n")
        sim = pynbody.load(snapName)
        if VERBOSE: sys.stdout.write("ngas = %e, ndark = %e, nstar = %e\n"%(len(sim.g), len(sim.d), len(sim.s)))
        # If there aren't any stars yet, skip the timestep.
        if len(sim.s) == 0:
            sys.stdout.write("Skipping " + snapName + " since there aren't any stars in it.\n")
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
        #
        #Step 3: Cut out the specified radius.
        #
        if VERBOSE: sys.stdout.write("Finding the largest halo\n")
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
        #
        #Step 4: Write the snapshot section to the disk in std tipsy and ASCII formats.
        #
        os.chdir(WORKING_DIR)
        # output filename format: <galaxy name>.<time step>.<diameter>kpc.phys|sim.stdtipsy|ascii
        SNAPFILE = GALAXY_NAME + "." + timeStep + "." + cutDiameter + "kpc"
        SNAPFILE += (".phys" if PHYS else ".sim")
        SNAPFILESTD = SNAPFILE + ".stdtipsy"
        SNAPFILEASC = SNAPFILE + ".ascii"
        sys.stdout.write("Writing " + SNAPFILESTD + "\n")
        cut.write(filename=SNAPFILESTD, fmt=pynbody.tipsy.TipsySnap)
        sys.stdout.write("Writing " + SNAPFILEASC + "\n")
        cmd = "cat " + SNAPFILESTD + " | std2ascii > " + SNAPFILEASC
        FNULL = open(os.devnull, "w")
        p1 = subprocess.Popen(cmd, shell=True, stdout=FNULL, stderr=FNULL)
        p1.wait()
        if p1.returncode != os.EX_OK:
            sys.stderr.write("Error running std2ascii. Perhaps it's not in your PATH?\n")
            FNULL.close()
            sys.exit(1)
        #
        #Step 5: Run SMOOTH.
        #
        sys.stdout.write("Writing smoothing lengths to smooth.hsm\n")
        cmd = "smooth hsmooth < " + SNAPFILESTD
        p2 = subprocess.Popen(cmd, shell=True, stdout=FNULL, stderr=FNULL)
        p2.wait()
        if p2.returncode != os.EX_OK:
            sys.stderr.write("Error running smooth. Perhaps it's not in your PATH?\n")
            FNULL.close()
            sys.exit(1)
        os.remove(SNAPFILESTD) # Sunrise just needs the ASCII version
        FNULL.close()
        #
        #Step 6: Generate Sunrise config files.
        #
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
        #
        #Step 7: Move the files into the final directory, and write out job submission commands.
        #
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
    # end for loop over time steps
    # write a file, runall-<galaxy name>-t<min timestep>[-<max timestep>].sh in OUT_DIR that will start sfrhist for all the timesteps
    os.chdir(OUT_DIR)
    listOfTimestepInts = [int(x) for x in listOfTimesteps]
    minTimestep = min(listOfTimestepInts)
    maxTimestep = max(listOfTimestepInts)
    jobRunner = "runall-" + GALAXY_NAME + "-t" + format(minTimestep, '05')
    if len(listOfTimesteps) > 1: jobRunner += "-" + format(maxTimestep, '05')
    jobRunner += ".sh"
    sys.stdout.write("Writing a script to start jobs for each timestep: " + jobRunner + "\n")
    with open(jobRunner, "w") as f:
        f.write("#!/bin/bash\n\n")
        for timeStep in listOfTimesteps:
            f.write(RUN_DIR + GALAXY_NAME + "-" + timeStep + "-sunrise/runsfrhist.sh\n")
    os.chmod(jobRunner, 0744)
    # make a gzipped tarball of the whole output directory
    if TARBALL:
        sys.stdout.write("Tarring up the output directory.\n")
        cmd = "tar czf " + jobRunner[7:-3] + ".tgz *"
        p3 = subprocess.Popen(cmd, shell=True)
        p3.wait()
        if p3.returncode != os.EX_OK:
            sys.stderr.write("Error encountered while tarring output directory!\n")
            sys.exit(1)
    if VERBOSE: sys.stdout.write("Finished. You're ready to run Sunrise!\n\n")
    sys.exit(0)
    
if __name__=="__main__":
    main()

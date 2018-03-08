#!/usr/bin/env python
# -*- coding: utf-8 -*-

# python run script for the event generator fortran executable
# generates a config file and runs the executable using it
# runs locally or submits a job to the lxplus/iridis/qmul batch system
# declan.millar@cern.ch

import os
import StringIO
import argparse
import subprocess
import sys
import glob
import socket
import time
import fnmatch

parser = argparse.ArgumentParser(description="generate ttbar events")
parser.add_argument("-a", "--include_a", help = "include photon mediated interaction", default = True, action = "store_false")
parser.add_argument("-b", "--include_background", help = "include tt background", default = False, action = "store_true")
parser.add_argument("-c", "--cut", help = "apply detector cuts", default = False, action = "store_true")
parser.add_argument("-d", "--include_dd", help = "include down-down initiated interactions", default = False, action = "store_true")
parser.add_argument("-E", "--energy", help = "set collider energy", type = int, default = 13)
parser.add_argument("-e", "--nevents", help = "set number of events", type = int, default = 10000)
parser.add_argument("-F", "--flatten_integrand", help = "flatten resonances in integrand", default = True, action = "store_false")
parser.add_argument("-f", "--final_state", help = "set final state id", type = int, default = 1)
parser.add_argument("-g", "--include_gg", help = "include gluon-gluon initiated interactions", default = False, action = "store_true")
parser.add_argument("-H", "--lhef", help = "write events to an LHEF file", default = True, action = "store_false")
parser.add_argument("-I", "--interference", help = "set interference (default = full)", type = int, default = 1)
parser.add_argument("-i", "--index", help = "append specified file index", default = "")
parser.add_argument("-j", "--job", help = "submit as a batch job", default = True,  action = "store_false")
parser.add_argument("-k", "--walltime", help = "set walltime in 'hh:mm:ss' format", default = "60:00:00")
parser.add_argument("-L", "--energy_low", help = "set an upper limit on the collider energy", type = int, default = 0)
parser.add_argument("-M", "--multichannel", help = "use multichannel integration", default = False, action = "store_true")
parser.add_argument("-m", "--model", help = "set physics model", default = "SM")
parser.add_argument("-N", "--iterations", help = "number of VAMP iterations", type = int, default = 16)
parser.add_argument("-n", "--npoints", help = "number of VAMP calls", type = int, default = 16000000)
parser.add_argument("-o", "--overwrite", help = "overwrite existing grid file", default = False, action = "store_true")
parser.add_argument("-P", "--pdf", help = "structure_functions", type = int, default = 11)
parser.add_argument("-p", "--ppbar", help = "set initial state: 0 = pp, 1 = pp~", type = int, default = 0)
parser.add_argument("-Q", "--queue", help = "lxbatch queue", default = "1nw")
parser.add_argument("-q", "--include_qq", help = "include quark-quark initiated interactions", default = False, action = "store_true")
parser.add_argument("-r", "--use_rambo", help = "use RAMBO for phase space", default = False, action = "store_true")
parser.add_argument("-s", "--include_signal", help = "include tt signal", default = True,  action = "store_false")
parser.add_argument("-T", "--ntuple", help = "write events to ROOT n-tuple", default = False,  action = "store_true")
parser.add_argument("-t", "--tag", help = "add a name tag to output files", default = "")
parser.add_argument("-U", "--energy_up", help = "set an upper limit on the collider energy", type = int, default = 0)
parser.add_argument("-u", "--include_uu", help = "include up-up initiated interactions", default = False, action = "store_true")
parser.add_argument("-v", "--verbose", help = "print extra run information", default = False, action = "store_true")
parser.add_argument("-W", "--use_nwa", help = "use Narrow Width Approximation", default = False, action = "store_true")
parser.add_argument("-w", "--unweighted", help = "generate unweighted events", default = True,  action = "store_false")
parser.add_argument("-x", "--include_x", help = "include Z' boson mediated interactions", default = True,  action = "store_false")
parser.add_argument("-z", "--include_z", help = "include Z boson mediated interaction", default = True,  action = "store_false")
args = parser.parse_args()

hostname = socket.gethostname()
home_directory = "."
data_directory = "."
if "Lorkhan" in hostname:
    home_directory = "/Users/declan/Projects/"
    data_directory = "/Users/declan/Data/"
elif "lxplus" in hostname:
    home_directory = "/afs/cern.ch/user/d/demillar/"
    data_directory = "/afs/cern.ch/work/d/demillar/"
elif "cyan" in hostname:
    home_directory = "/home/dam1g09/"
    data_directory = "/scratch/dam1g09/"
elif "heppc" in hostname:
    home_directory = "/users/millar/"
    data_directory = "/data/millar/"
else:
    exit("ERROR: unknown host")
run_directory = home_directory + "prophet/"
data_directory = data_directory + "zprime/"

# check model file exists
if os.path.isfile(run_directory + "Models/%s.mdl" % args.model) is False:
    print "ERROR: %s is not a valid model" % args.model
    print "available models:"
    for name in glob.glob("Models/*.mdl"):
        print "\t", name[7:-4]
    sys.exit()

if args.energy < 0:
    sys.exit("ERROR: collider energy must be positive definite")

if args.npoints < 2:
    sys.exit("ERROR: must have at least 2 sampling points")

if args.energy_low < 0 or args.energy_up < 0:
    sys.exit("ERROR: collider energy bounds must be positive definite")

if args.energy_low > args.energy or args.energy_up > args.energy:
    sys.exit("ERROR: energy bounds cannot exceed collider energy")

if args.energy_low > 0 and args.energy_up > 0 and args.energy_up <= args.energy_low:
    sys.exit("ERROR: energy upper bound must be greater than lower bound")

if args.interference < 0 or args.interference > 4:
    sys.exit("ERROR: interference must be between 0 and 4")

if args.pdf < 1 or args.pdf > 11:
    sys.exit("ERROR: PDF ID must be between 1 and 11")

if args.include_background == False and args.include_signal == False:
    sys.exit("ERROR: signal and background both off")

initial_states = 0
if args.include_gg: initial_states += 1
if args.include_qq: initial_states += 1
if args.include_dd: initial_states += 1
if args.include_uu: initial_states += 1

if args.lhef and initial_states > 1:
    sys.exit("ERROR: currently when outputting to LHEF, only one initial state can be active")

if initial_states == 0:
    print "WARNING: no initial states specified; will integrate phase space only"

if not ("lxplus" in hostname or "cyan" in hostname):
    args.job = False

if args.model == "SM":
    args.include_x = False

if args.final_state == -1:
    args.include_gg = False
    args.include_qq = False

if args.include_uu == False and args.include_dd == False:
    args.include_a = False
    args.include_z = False
    args.include_x = False

if args.final_state == -1 and args.energy_low == 0:
    # avoid divergence in Drell-Yan
    args.energy_low = 0.5

if args.use_rambo or args.include_background:
    args.flatten_integrand = False

if args.final_state < 1:
    args.use_nwa = False

executable = "generator"
options = ""

if args.ppbar == 1:
    options += "_ppbar"

pdf = ""
if args.pdf ==  1: pdf = "CTEQ6M"
if args.pdf ==  2: pdf = "CTEQ6D"
if args.pdf ==  3: pdf = "CTEQ6L"
if args.pdf ==  4: pdf = "CTEQ6L1"
if args.pdf ==  5: pdf = "MRS9901"
if args.pdf ==  6: pdf = "MRS9902"
if args.pdf ==  7: pdf = "MRS9903"
if args.pdf ==  8: pdf = "MRS9904"
if args.pdf ==  9: pdf = "MRS9905"
if args.pdf == 10: pdf = "CT14LN"
if args.pdf == 11: pdf = "CT14LL"

if not args.unweighted:
    args.lhef = False

if args.final_state < 2:
    args.include_background = False

if args.include_background:
    flatten_integrand = False

# include non-default major options in output file name
if args.interference != 1: options += "_int%i" % args.interference
if args.use_nwa: options += "_nwa"
if args.multichannel: options += "_multi"
if args.cut: options += "_cut"
if args.use_rambo: options += "_rambo"
if args.include_background == False and args.flatten_integrand == False: options += "_unmapped"
if args.energy_low != 0 or args.energy_up != 0: options += "_%s-%s" % (args.energy_low, args.energy_up)
if len(args.tag) > 0: options += "_" + args.tag

initial_partons = ""
if args.include_gg: initial_partons += "gg"
if args.include_qq: initial_partons += "qq"
if args.include_dd: initial_partons += "dd"
if args.include_uu: initial_partons += "uu"
if initial_partons == "": initial_partons = "phase-space"
initial_partons += "-"

intermediates = ""
if args.include_a: intermediates += "A"
if args.include_z: intermediates += "Z"
if args.include_x: intermediates += "X"
if len(intermediates) > 0: intermediates = intermediates + "-"

grid_state = "tt-bbllvv"
gridproc = initial_partons + intermediates + grid_state

final_state = ""
if   args.final_state == -1: final_state = "ll"
elif args.final_state ==  0: final_state = "tt"
elif args.final_state ==  1: final_state = "tt-bbllvv"
elif args.final_state == 11: final_state = "tt-bbeevv"
elif args.final_state == 22: final_state = "tt-bbmumuvv"
elif args.final_state == 12: final_state = "tt-bbemuvv"
elif args.final_state == 21: final_state = "tt-bbmuevv"
elif args.final_state == 33: final_state = "bbtatavv"

if initial_partons == "-" and args.final_state > 0: process = "2-6-phase-space"
elif initial_partons == "-" and args.final_state < 1: process = "2-2-phase-space"
else: process = initial_partons + intermediates + final_state

grid_name = '%s_%s_%sTeV_%s%s' % (gridproc, args.model, str(args.energy), pdf, options)
events_name = '%s_%s_%sTeV_%s%s' % (process, args.model, str(args.energy), pdf, options)

grid_path = data_directory + grid_name

datafiles = [ f for f in os.listdir(data_directory) if os.path.isfile(os.path.join(data_directory, f)) ]
filtered = fnmatch.filter(datafiles, events_name + "_??.lhef")
# print filtered
new_index = int(args.index) if args.index != "" else len(filtered) + 1
events_name = events_name + "_%03d" % new_index
events_path = data_directory + events_name

if os.path.isdir(data_directory) is False:
    sys.exit("ERROR: specified run directory '%s' does not exist" % run_directory)
    sys.exit("ERROR: specified data directory '%s' does not exist" % data_directory)

if args.multichannel: grid = ".grids"
else: grid = ".grid"

if args.unweighted: wgt = ""
else: wgt = ".wgt"

events_path = events_path + wgt
grid_file = "%s%s" % (grid_path, grid)
xsec_file = "%s.txt" % (grid_path)

print "gridfile = ", grid_file

if os.path.isfile(grid_file): new_grid = False
else: new_grid = True

if args.overwrite: new_grid = True

ntuple_file = "%s.root" % (events_path)
lhe_file = "%s.lhef" % (events_path)

if new_grid:
    config_name = '%s.cfg' % (grid_path)
    logfile = "%s.log" % (grid_path)
    handler_name = "%s.sh" % (grid_name)
else:
    config_name = '%s.cfg' % (events_path)
    logfile = "%s.log" % (events_path)
    handler_name = "%s.sh" % (events_name)

config = StringIO.StringIO()
print >> config, '%r    ! ntuple'             % args.ntuple
print >> config, '%r    ! lhef'               % args.lhef
print >> config, '%r    ! new_grid'           % new_grid
print >> config, '%s'                         % grid_file
print >> config, '%s'                         % xsec_file
print >> config, '%s'                         % logfile
print >> config, '%s'                         % lhe_file
print >> config, '%s'                         % ntuple_file
print >> config, '%i    ! initial state'      % args.ppbar
print >> config, '%i    ! final state'        % args.final_state
print >> config, '%s    ! model'              % args.model
print >> config, '%i    ! pdf'                % args.pdf
print >> config, '%r    ! include signal'     % args.include_signal
print >> config, '%r    ! include background' % args.include_background
print >> config, '%r    ! include gg'         % args.include_gg
print >> config, '%r    ! include qq'         % args.include_qq
print >> config, '%r    ! include uu'         % args.include_uu
print >> config, '%r    ! include dd'         % args.include_dd
print >> config, '%r    ! include a'          % args.include_a
print >> config, '%r    ! include z'          % args.include_z
print >> config, '%r    ! include x'          % args.include_x
print >> config, '%i    ! interference'       % args.interference
print >> config, '%r    ! use nwa'            % args.use_nwa
print >> config, '%i.d3 ! energy'             % args.energy
print >> config, '%i    ! iterations'         % args.iterations
print >> config, '%i    ! npoints'            % args.npoints
print >> config, '%i    ! nevents'            % args.nevents
print >> config, '%r    ! unweighted'         % args.unweighted
print >> config, '%r    ! use rambo'          % args.use_rambo
print >> config, '%r    ! map phase space'    % args.flatten_integrand
print >> config, '%r    ! multichannel'       % args.multichannel
print >> config, '%r    ! verbose mode'       % args.verbose
print >> config, '%i.d3 ! energy low'         % args.energy_low
print >> config, '%i.d3 ! energy up'          % args.energy_up
print >> config, '%r    ! batch mode'         % args.job
print >> config, '%r    ! detector cuts'      % args.cut

try:
    with open('%s' % config_name, 'w') as config_file:
        config_file.write(config.getvalue())
        print "config = %s" % config_name
except OSError: sys.exit("ERROR: Cannot write to %s" % config_name)

if args.job:
    handler = StringIO.StringIO()
    if "lxplus" in hostname:
        print >> handler, "#!/bin/bash"
        print >> handler, "source /afs/cern.ch/sw/lcg/external/gcc/6.2/x86_64-slc6/setup.sh"
        print >> handler, "source /cvmfs/sft.cern.ch/lcg/releases/ROOT/6.10.04-22868/x86_64-slc6-gcc62-opt/bin/thisroot.sh"
        # print >> handler, "export PATH=/afs/cern.ch/sw/lcg/external/openmpi/1.8.1/x86_64-slc6-gcc48-opt/bin:$PATH"
        # print >> handler, "export LD_LIBRARY_PATH=/afs/cern.ch/sw/lcg/external/openmpi/1.8.1/x86_64-slc6-gcc48-opt/lib:$LD_LIBRARY_PATH"
        print >> handler, "export PATH=/cvmfs/sft.cern.ch/lcg/releases/Python/2.7.13-597a5/x86_64-slc6-gcc62-opt/bin:$PATH"
        print >> handler, "export LD_LIBRARY_PATH=/cvmfs/sft.cern.ch/lcg/releases/Python/2.7.13-597a5/x86_64-slc6-gcc62-opt/lib:$LD_LIBRARY_PATH"
        print >> handler, "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/afs/cern.ch/user/d/demillar/root-tuple/lib"
        print >> handler, "source /afs/cern.ch/sw/IntelSoftware/linux/setup.sh"
        print >> handler, "source /afs/cern.ch/sw/IntelSoftware/linux/x86_64/xe2017/bin/compilervars.sh intel64"
        print >> handler, "cd %s" % run_directory
        print >> handler, "%s/bin/%s < %s" % (run_directory, executable, config_name)
    elif "cyan" in hostname:
        print "walltime = %s" % args.walltime
        print >> handler, "#!/bin/bash"
        print >> handler, "source /home/dam1g09/.bash_profile"
        print >> handler, "module load gcc/6.1.0"
        print >> handler, "module load openmpi/2.0.2/gcc"
        print >> handler, "module load intel/2017"
        print >> handler, "module load intel/mpi/2017"
        print >> handler, "cd %s" % run_directory
        print >> handler, "%sbin/%s < %s > %s" % (run_directory, executable, config_name, logfile)
    elif "heppc" in hostname:
        print "h_rt = %s" % args.walltime
        print >> handler, "#!/bin/bash"
        print >> handler, "source /users/millar/.bash_profile"
        print >> handler, "cd %s" % run_directory
        print >> handler, "%sbin/%s < %s" % (run_directory, executable, config_name)
    else:
        sys.exit("ERROR: hostname not recognised")

    print >> handler, "gzip -v9 %s >> %s" % (lhe_file, logfile)

    try:
        with open('%s' % handler_name, 'w') as handler_file:
            handler_file.write(handler.getvalue())
        print "handler = %s" % handler_name
    except OSError:
        sys.exit("ERROR: cannot write handler file")

    subprocess.call("chmod a+x %s" % handler_name, shell = True)
    print "submitting batch job ..."
    if "lxplus" in hostname: print "queue = ", args.queue
    if "lxplus" in hostname: subprocess.call('bsub -q %s -o %s %s%s' % (args.queue, logfile, run_directory, handler_name), shell = True)
    elif "cyan" in hostname: subprocess.call('qsub -l walltime=%s %s%s' % (args.walltime, run_directory, handler_name), shell = True)
    elif "heppc" in hostname: subprocess.call('qsub -l h_rt=%s %s%s' % (args.walltime, run_directory, handler_name), shell = True)
    else:
        sys.exit("ERROR: hostname not recognised")
else:
    # if "cyan" in hostname:
        # print "loading openmpi module ..."
        # subprocess.call("module load gcc/6.1.0; module load openmpi/2.0.2/gcc")
    subprocess.call(run_directory + "bin/%s < %s | tee %s" % (executable, config_name, logfile), shell = True)
    # subprocess.call("gzip -v9 %s >> %s" % (lhe_file, logfile), shell = True)

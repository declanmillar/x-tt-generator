#!/usr/bin/env python

# Python run script for the Zprime program.
# Generates a configuration file, then runs the zprime executable using it, either locally or by submission to the LXPLUS or Iridis batch system.
# Do './generate.py -h' for help.
# Author: Declan Millar, d.millar@soton.ac.uk.

import os, StringIO, re, optparse, subprocess, time, sys, random, glob, socket, string

parser = optparse.OptionParser()

parser.add_option("-t", "--tag", default = "", type = "string", help = "add a name tag to output files")
parser.add_option("-b", "--batch", default = True, action = "store_false", help = "run in batch mode")
parser.add_option("-W", "--walltime", default = "60:00:00", action = "store", help = "walltime 'hh:mm:ss'")
parser.add_option("-Q", "--queue", default = "8nh", action = "store", help = "lxbatch queue'")
parser.add_option("-T", "--ntuple_out", default = True, action = "store_false", help = "write events to ROOT ntuple")
parser.add_option("-L", "--lhef_out", default = False, action = "store_true", help = "write events to lhef file")

parser.add_option("-m", "--model", default = "SM", action = "store", help = "set model")

parser.add_option("-p", "--initial_state", default = 0, const = 1, action = "store_const", help = "switch to p-pbar collisions")
parser.add_option("-g", "--include_gg", default = True, action = "store_false", help = "exclude gluon-gluon interactions")
parser.add_option("-q", "--include_qq", default = True, action = "store_false", help = "exclude quark-quark interactions")
parser.add_option("-G", "--include_g", default = False, action = "store_true", help = "exclude gluon mediated interactions")
parser.add_option("-A", "--include_a", default = True, action = "store_false", help = "exclude photon mediated interaction")
parser.add_option("-Z", "--include_z", default = True, action = "store_false", help = "exclude Z boson mediated interaction")
parser.add_option("-X", "--include_x", default = True, action = "store_false", help = "exclude Z' boson mediated interactions")
parser.add_option("-f", "--final_state", default = "tt-bbllvv", action = "store", help = "set final state: ll, tt, tt-bbllvv, bbeevv, bbemuvevm, bbmumuvmvm, bbtatavtvt")

parser.add_option("-E", "--collider_energy", default = 13, action = "store", help = "collider energy")
parser.add_option("-P", "--structure_function", default = 4, type = "int", help = "structure_functions set: 1 = CTEQ6M; 2 = CTEQ6D; 3 = CTEQ6L; 4 = CTEQ6L1; ...")
parser.add_option("-S", "--include_signal", default = 1, const = 0, action = "store_const", help = "include tt signal")
parser.add_option("-B", "--include_background", default = 0, const = 1, action = "store_const", help = "include tt background")

parser.add_option("-i", "--interference", default = 2, type = "int", help = "specify interference")
parser.add_option("-w", "--use_nwa", default = False, action = "store_true", help = "use Narrow Width Approximation")
parser.add_option("-l", "--ecm_low", default = 0, type = "int", help = " Ecm lower limit")
parser.add_option("-u", "--ecm_up", default = 0, type = "int", help = "Ecm upper limit")

# Monte Carlo options
parser.add_option("-F", "--fixed_seed", default = False, action = "store_true", help = "use fixed seed for random number generator")
parser.add_option("-n", "--vegas_points", default = 5000000, type = "int", help = "number of VEGAS points")
parser.add_option("-N", "--itmx", default = 5, type = "int", help = "maximum number of VEGAS iterations")
parser.add_option("-s", "--symmetrise", default = True, action = "store_false", help = "symmetrise phase space x1<->x2")
parser.add_option("-R", "--use_rambo", default = False, action = "store_true", help = "use RAMBO for PS")
parser.add_option("-M", "--map_phase_space", default = True, action = "store_false", help = "flatten Breit-Wigners in integrand for manual phase space")

# Debug option
parser.add_option("-O", "--phase_space_only", default = False, action = "store_true", help = "Set |M|^2 =    1")
parser.add_option("-v", "--verbose", default = False, action = "store_true", help = "Run in verbose mode.")

(option, args) = parser.parse_args()

hostname = socket.gethostname()

final_state = str(option.final_state)
model_name = str(option.model)
try:
    collider_energy = int(option.collider_energy)
except ValueError:
    sys.exit("Error: invalid collider energy '%s [TeV]'.\n%s" % (args[2],usage))
try:
    ncall = int(option.vegas_points)
except ValueError:
    sys.exit("Error: invalid VEGAS points '%s'.\n%s" % (args[3],usage))

# Check arguments are valid
if os.path.isfile("Models/%s.mdl" % model_name) is False:
    sys.exit("%s is not a valid model.\n Available model files: %s" % (model_name, glob.glob("Models/*.mdl")))

if collider_energy < 0:
    sys.exit("Error: collider energy must be positive definite.\n" % usage)

if ncall < 2:
    sys.exit("Error: Must have at least 2 VEGAS points.\n%s" % usage)

# Check options are valid
if "Sunder" in hostname:
    option.batch = False

if option.ecm_low < 0 or option.ecm_up < 0:
    sys.exit("Error: COM energy must be positive definite")

if option.ecm_low > collider_energy or option.ecm_up > collider_energy:
    sys.exit("Error: COM energy cannot exceed collider energy")

if (option.ecm_low > 0 and option.ecm_up > 0 and option.ecm_up <= option.ecm_low):
    sys.exit("Error: E_CM up must be greater than E_CM low")

if option.batch and ("lxplus" in hostname or "iridis" in hostname):
    sys.exit("Error: Must be on lxplus or iridis to submit a batch job.")

if option.interference < 0 or option.interference > 4:
    sys.exit("Error: interference must be from 0-4.")

if option.structure_function < 1 or option.structure_function > 9:
    sys.exit("Error: structure_function ID must be from 1 to 9.")

# Modify configuration for consistency
if model_name == "SM":
    option.include_x = False

if final_state == "ll":
    option.include_g = False
    option.include_gg = False

if option.phase_space_only:
    option.include_g = False
    option.include_a = False
    option.include_z = False
    option.include_x = False

if option.use_rambo or option.include_background:
    option.map_phase_space = False

if final_state == "ll" or final_state == "tt":
    option.use_nwa = False

# Default iseed
seed = 12345 if option.fixed_seed else random.randint(0,100000)

# Strings
executable = "zprime"

options = ""

if option.initial_state == 1:
    options += "p"

if option.structure_function != 4:
    options += "S%s" % option.structure_function

if option.interference == 0:
    options += "i0"
elif option.interference == 1:
    options += "i1"
elif option.interference == 3:
    options += "i3"
elif option.interference == 4:
    options += "i4"

elif option.use_nwa:
    options += "w"

if option.ecm_low != 0 and option.ecm_up != 0:
    options += "%s-%s" % (option.ecm_low, option.ecm_up)

# exclude divergence
if final_state == "ll" and option.ecm_low == 0:
    option.ecm_low = 0.5

if option.fixed_seed:
    options += "f"

if not option.symmetrise:
    options += "s"

if option.use_rambo:
    options += "R"
if not option.map_phase_space:
    options += "M"

if len(options) > 0:
    options = "_" + options

if len(option.tag) > 0:
    options += "_" + option.tag

npoints = str(ncall)
if "000000" in npoints:
    npoints = "M".join(npoints.rsplit("000000", 1))
if "000" in npoints:
    npoints = "k".join(npoints.rsplit("000", 1))

energy_collider = "_" + str(collider_energy) if option.collider_energy != 13 else ""

if option.include_g == False:
    option.include_gg = False
if option.include_gg == False:
    option.include_g = False

initial_partons = ""
if option.include_qq:
    initial_partons += "qq"
if option.include_gg:
    initial_partons += "gg"

intermediates = ""
if option.include_g:
    intermediates += "G"
if option.include_a:
    intermediates += "A"
if option.include_z:
    intermediates += "Z"
if option.include_x:
    intermediates += "X"

filename = '%s_%s-%s-%s%s%s_%sx%s' % (model_name, initial_partons, intermediates, final_state, energy_collider, options, option.itmx, npoints)

# Generate fortran friendly configuration
if final_state == "ll":
    final_state_id = -1
elif final_state == "tt":
    final_state_id = 0
elif final_state == "tt-bbllvv":
    final_state_id = 1
elif final_state == "bbeeveve":
    final_state_id = 2
elif final_state == "bbemuvevm":
    final_state_id = 3
elif final_state == "bbmuvmvm":
    final_state_id = 4
elif final_state == "bbtatavtvt":
    final_state_id = 5
else:
    sys.exit("Error: unavailable final state '%s'." % final_state)

if final_state_id < 2:
    include_background = False

if include_background:
    map_phase_space = False

# Logfile
run_directory = "."
data_directory = "."
if "Sunder" in hostname:
    data_directory = "/Users/declan/Data/Zprime"
elif "lxplus" in hostname:
    run_directory = "/afs/cern.ch/user/d/demillar/zprime-top-generator"
    data_directory = "/afs/cern.ch/work/d/demillar/zprime"
elif "iridis" in hostname:
    run_directory = "/home/dam1g09/zprime-top-generator"
    data_directory = "/scratch/dam1g09/zprime"

if os.path.isdir(data_directory) is False:
    sys.exit("The target data directory '%s' does not exist" % data_directory)

config_name = '%s/%s.cfg' % (data_directory, filename)
logfile = "%s/%s.log" % (data_directory, filename)
handler_name = "%s.sh" % filename
ntuple_file = "%s/%s.root" % (data_directory, filename)
lhe_file = "%s/%s.lhe" % (data_directory, filename)

config = StringIO.StringIO()
print >> config, '%i ! ntuple_out' % option.ntuple_out
print >> config, '%i ! lhef_out' % option.lhef_out
print >> config, '%s' % ntuple_file
print >> config, '%s' % lhe_file
print >> config, '%s' % logfile
print >> config, '%i ! initial_state' % option.initial_state
print >> config, '%i ! final_state' % final_state_id
print >> config, '%s ! model_name' % model_name
print >> config, '%i ! istructure' % option.structure_function
print >> config, '%i ! include_signal' % option.include_signal
print >> config, '%i ! include_background' % option.include_background
print >> config, '%i ! include_g' % option.include_g
print >> config, '%i ! include_a' % option.include_a
print >> config, '%i ! include_z' % option.include_z
print >> config, '%i ! include_x' % option.include_x
print >> config, '%i ! include_gg' % option.include_gg
print >> config, '%i ! include_qq' % option.include_qq
print >> config, '%i ! phase_space_only' % option.phase_space_only
print >> config, '%i ! interference' % option.interference
print >> config, '%i ! use_nwa' % option.use_nwa
print >> config, '%i.d3 ! ecm_col' % collider_energy
print >> config, '%i ! iseed' % seed
print >> config, '%i ! itmx' % option.itmx
print >> config, '%i ! ncall' % ncall
print >> config, '-1.d0 ! acc'
print >> config, '%i ! use rambo' % option.use_rambo
print >> config, '%i ! map phase space' % option.map_phase_space
print >> config, '%i ! symmetrise' % option.symmetrise
print >> config, '%i ! verbose mode' % option.verbose
print >> config, '%f ! ecm_low' % (option.ecm_low*1000)
print >> config, '%f ! ecm_up' % (option.ecm_up*1000)

try:
    with open('%s' % config_name,'w') as config_file:
        config_file.write(config.getvalue())
        print " Config: %s" % config_name
except IOError:
    sys.exit(" Error: cannot write to %s" % config_name)

if option.batch:
    handler = StringIO.StringIO()
    if "lxplus" in hostname:
        print >> handler, "export LD_LIBRARY_PATH=/afs/cern.ch/user/d/demillar/.RootTuple:$LD_LIBRARY_PATH"
        print >> handler, "source /afs/cern.ch/sw/lcg/external/gcc/4.8/x86_64-slc6/setup.sh"
        print >> handler, "cd %s" % run_directory
        print >> handler, '%s/Binary/%s < %s' % (run_directory, executable, config_name)
        # print >> handler, 'mv LSFJOB_* Jobs'
    if "iridis" in hostname:
        print "walltime = %s" % option.walltime
        print >> handler, "#!/bin/bash"
        print >> handler, "module load gcc/4.8.1; source /local/software/cern/root_v5.34.14/bin/thisroot.sh"
        print >> handler, "export LD_LIBRARY_PATH=/home/dam1g09/.RootTuple:$LD_LIBRARY_PATH"
        print >> handler, "cd %s" % run_directory
        # print >> handler, "cd $PBS_O_WORKDIR"
        print >> handler, '%s/Binary/%s < %s' % (run_directory, executable, config_name)
    print >> handler, 'rm -- "$0"'
    try:
        with open('%s' % handler_name, 'w') as handler_file:
            handler_file.write(handler.getvalue())
        print "Handler file written to %s.sh." % filename
    except IOError:
        sys.exit(" Error: cannot write handler file.")

    subprocess.call("chmod a+x %s.sh" % filename, shell = True)
    print "Submitting batch job."
    if "lxplus" in hostname: subprocess.call('bsub -q %s %s/%s.sh' % (option.queue, run_directory, filename), shell = True)
    if "iridis" in hostname: subprocess.call('qsub -l walltime=%s %s/%s.sh' % (option.walltime, run_directory, filename), shell = True)

else:
    if "lxplus" in hostname:
        print "Sourcing ROOT..."
        subprocess.call("source /afs/cern.ch/sw/lcg/external/gcc/4.8/x86_64-slc6/setup.sh", shell = True)
        print "Adding RootTuple libraries to library path..."
        subprocess.call("export LD_LIBRARY_PATH=/afs/cern.ch/user/d/demillar/.RootTuple:$LD_LIBRARY_PATH", shell = True)
    print " Program will run locally..."
    subprocess.call("./Binary/%s < %s" % (executable, config_name), shell = True)

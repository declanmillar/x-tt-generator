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
parser.add_option("-T", "--ntuple", default = True, action = "store_false", help = "write events to ROOT ntuple")
parser.add_option("-H", "--lhef", default = False, action = "store_true", help = "write events to lhef file")
parser.add_option("-m", "--model", default = "SM", action = "store", help = "set model")
parser.add_option("-p", "--initial_state", default = 0, const = 1, action = "store_const", help = "switch to p-pbar collisions")
parser.add_option("-g", "--include_gg", default = True, action = "store_false", help = "gluon-gluon interactions")
parser.add_option("-q", "--include_qq", default = True, action = "store_false", help = "quark-quark interactions")
parser.add_option("-u", "--include_uu", default = True, action = "store_false", help = "up-up interactions")
parser.add_option("-d", "--include_dd", default = True, action = "store_false", help = "down-down interactions")
parser.add_option("-A", "--include_a", default = True, action = "store_false", help = "exclude photon mediated interaction")
parser.add_option("-Z", "--include_z", default = True, action = "store_false", help = "exclude Z boson mediated interaction")
parser.add_option("-X", "--include_x", default = True, action = "store_false", help = "exclude Z' boson mediated interactions")
parser.add_option("-f", "--final_state", default = "1", type = "int", action = "store", help = "set final state")
parser.add_option("-E", "--collider_energy", default = 13, action = "store", help = "collider energy")
parser.add_option("-P", "--structure_function", default = 4, type = "int", help = "structure_functions")
parser.add_option("-S", "--include_signal", default = True, const = 0, action = "store_const", help = "include tt signal")
parser.add_option("-B", "--include_background", default = False, const = 0, action = "store_const", help = "include background")
parser.add_option("-i", "--interference", default = 2, type = "int", help = "specify interference")
parser.add_option("-w", "--use_nwa", default = False, action = "store_true", help = "use Narrow Width Approximation")
parser.add_option("-L", "--ecm_low", default = 0, type = "int", help = " Ecm lower limit")
parser.add_option("-U", "--ecm_up", default = 0, type = "int", help = "Ecm upper limit")
parser.add_option("-F", "--fixed_seed", default = False, action = "store_true", help = "use fixed seed")
parser.add_option("-n", "--vegas_points", default = 10000000, type = "int", help = "number of VEGAS points")
parser.add_option("-N", "--itmx", default = 5, type = "int", help = "maximum number of VEGAS iterations")
parser.add_option("-s", "--symmetrise", default = True, action = "store_false", help = "symmetrise phase space x1<->x2")
parser.add_option("-R", "--use_rambo", default = False, action = "store_true", help = "use RAMBO for PS")
parser.add_option("-M", "--map_phase_space", default = True, action = "store_false", help = "map phase space")
parser.add_option("-O", "--phase_space_only", default = False, action = "store_true", help = "Set |M|^2 =    1")
parser.add_option("-v", "--verbose", default = False, action = "store_true", help = "Run in verbose mode.")
parser.add_option("-c", "--cut", default = False, action = "store_true", help = "Apply fiducial cuts.")

(option, args) = parser.parse_args()
hostname = socket.gethostname()

if os.path.isfile("Models/%s.mdl" % option.model) is False:
    sys.exit("ERROR! %s is not a valid model.\n Available model files: %s" % (option.model, glob.glob("Models/*.mdl")))
if option.collider_energy < 0:
    sys.exit("ERROR! Collider energy must be positive definite.\n" % usage)
if option.vegas_points < 2:
    sys.exit("ERROR! Must have at least 2 VEGAS points.\n%s" % usage)
if option.ecm_low < 0 or option.ecm_up < 0:
    sys.exit("ERROR! Ecm bounds must be positive definite")
if option.ecm_low > option.collider_energy or option.ecm_up > option.collider_energy:
    sys.exit("ERROR! Ecm bounds cannot exceed collider energy")
if (option.ecm_low > 0 and option.ecm_up > 0 and option.ecm_up <= option.ecm_low):
    sys.exit("ERROR! Ecm upper bound must be greater than lower bound")
if option.interference < 0 or option.interference > 4:
    sys.exit("ERROR! Interference must be from 0-4.")
if option.structure_function < 1 or option.structure_function > 11:
    sys.exit("ERROR! Structure_function ID must be from 1 to 11.")
if option.include_background == False and option.include_signal == False:
    sys.exit("ERROR! Signal and background both off.")
if option.final_state < -1 or option.final_state > 3:

    sys.exit("ERROR! invalid final state id." % option.final_state)
initial_states = 0
if option.include_gg: initial_states += 1
if option.include_qq: initial_states += 1
if option.include_uu: initial_states += 1
if option.include_dd: initial_states += 1
if option.lhef and initial_states > 1:
    sys.exit("ERROR! When outputting to LHEF, only one initial state can be active.")


if not ("lxplus" in hostname or "cyan" in hostname):
    option.batch = False
if option.model == "SM":
    option.include_x = False
if option.final_state == -1:
    option.include_gg = False
    option.include_qq = False
if option.include_uu == False and option.include_dd == False:
    option.include_a = False
    option.include_z = False
    option.include_x = False
if option.phase_space_only:
    option.include_gg = False
    option.include_qq = False
    option.include_uu = False
    option.include_dd = False
    option.include_a = False
    option.include_z = False
    option.include_x = False
if option.final_state == -1 and option.ecm_low == 0:
    option.ecm_low = 0.5
if option.use_rambo or option.include_background:
    option.map_phase_space = False
if option.final_state < 1:
    option.use_nwa = False
seed = 12345 if option.fixed_seed else random.randint(0, 100000)

# Strings
executable = "zprime"
options = ""
if option.initial_state == 1:
    options += "_ppbar"
if option.structure_function != 4:
    options += "_pdf%s" % option.structure_function
if option.interference != 2:
    options += "_int%i" % option.interference
if option.use_nwa:
    options += "_nwa"
if option.fixed_seed:
    options += "_fixed"
if not option.symmetrise:
    options += "_nosym"
if option.use_rambo:
    options += "_R"
if  option.include_background == False and option.map_phase_space == False:
    options += "_M"
if option.ecm_low != 0 and option.ecm_up != 0:
    options += "_%s-%s" % (option.ecm_low, option.ecm_up)
if option.cut:
    options += "-fid"
if len(option.tag) > 0:
    options += "_" + option.tag
if option.final_state < 2:
    option.include_background = False
if option.include_background:
    map_phase_space = False

npoints = str(option.vegas_points)
if "000000" in npoints:
    npoints = "M".join(npoints.rsplit("000000", 1))
if "000" in npoints:
    npoints = "k".join(npoints.rsplit("000", 1))

energy_collider = "_" + str(option.collider_energy) if option.collider_energy != 13 else ""

initial_partons = ""
if option.include_qq:
    initial_partons += "qq"
if option.include_gg:
    initial_partons += "gg"
if option.include_uu:
    initial_partons += "uu"
if option.include_dd:
    initial_partons += "dd"

intermediates = ""
if option.final_state < 2:
    if option.include_a:
        intermediates += "A"
    if option.include_z:
        intermediates += "Z"
    if option.include_x:
        intermediates += "X"
else:
    if option.include_background == False and option.include_signal == True:
        intermediates += "tt"
    if option.include_background == True and option.include_signal == False:
        intermediates += "nott"

if len(intermediates) > 0:
    intermediates = intermediates + "-"

final_state = ""
if option.final_state == -1:
    final_state = "ll"
elif option.final_state == 0:
    final_state = "tt"
elif option.final_state == 1:
    final_state = "tt-bbllvv"
elif option.final_state == 2:
    final_state = "tt-blvbqq"
elif option.final_state == 11:
    final_state = "bbtatavtvt"
elif option.final_state == 12:
    final_state = "bbemuvevm"



filename = '%s_%s-%s%s%s%s_%sx%s' % (option.model, initial_partons, intermediates, final_state, energy_collider, options, option.itmx, npoints)

# Logfile
run_directory = "."
data_directory = "."
if "Sunder" in hostname:
    data_directory = "/Users/declan/Data/zprime"
elif "lxplus" in hostname:
    run_directory = "/afs/cern.ch/user/d/demillar/zprime-top-generator"
    data_directory = "/afs/cern.ch/work/d/demillar/zprime"
elif "cyan" in hostname:
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
print >> config, '%i ! ntuple' % option.ntuple
print >> config, '%i ! lhef' % option.lhef
print >> config, '%s' % ntuple_file
print >> config, '%s' % lhe_file
print >> config, '%s' % logfile
print >> config, '%i ! initial_state' % option.initial_state
print >> config, '%i ! final_state' % option.final_state
print >> config, '%s ! model' % option.model
print >> config, '%i ! istructure' % option.structure_function
print >> config, '%i ! include_signal' % option.include_signal
print >> config, '%i ! include_background' % option.include_background
print >> config, '%i ! include_gg' % option.include_gg
print >> config, '%i ! include_qq' % option.include_qq
print >> config, '%i ! include_uu' % option.include_uu
print >> config, '%i ! include_dd' % option.include_dd
print >> config, '%i ! include_a' % option.include_a
print >> config, '%i ! include_z' % option.include_z
print >> config, '%i ! include_x' % option.include_x
print >> config, '%i ! phase_space_only' % option.phase_space_only
print >> config, '%i ! interference' % option.interference
print >> config, '%i ! use_nwa' % option.use_nwa
print >> config, '%i.d3 ! ecm_col' % option.collider_energy
print >> config, '%i ! iseed' % seed
print >> config, '%i ! itmx' % option.itmx
print >> config, '%i ! ncall' % option.vegas_points
print >> config, '-1.d0 ! acc'
print >> config, '%i ! use rambo' % option.use_rambo
print >> config, '%i ! map phase space' % option.map_phase_space
print >> config, '%i ! symmetrise' % option.symmetrise
print >> config, '%i ! verbose mode' % option.verbose
print >> config, '%f ! ecm_low' % (option.ecm_low*1000)
print >> config, '%f ! ecm_up' % (option.ecm_up*1000)
print >> config, '%i ! apply fiducial cuts' % option.cut

try:
    with open('%s' % config_name,'w') as config_file:
        config_file.write(config.getvalue())
        print " Config: %s" % config_name
except IOERROR:
    sys.exit("ERROR! Cannot write to %s" % config_name)

if option.batch:
    handler = StringIO.StringIO()
    if "lxplus" in hostname:
        print >> handler, "export LD_LIBRARY_PATH=/afs/cern.ch/user/d/demillar/.RootTuple:$LD_LIBRARY_PATH"
        # print >> handler, "source /afs/cern.ch/sw/lcg/external/gcc/4.8/x86_64-slc6/setup.sh"
        print >> handler, "cd %s" % run_directory
        print >> handler, '%s/Binary/%s < %s' % (run_directory, executable, config_name)
    elif "cyan" in hostname:
        print "walltime = %s" % option.walltime
        print >> handler, "#!/bin/bash"
        print >> handler, "module load gcc/4.8.1; source /local/software/cern/root_v5.34.14/bin/thisroot.sh"
        print >> handler, "export LD_LIBRARY_PATH=/home/dam1g09/.RootTuple:$LD_LIBRARY_PATH"
        print >> handler, "cd %s" % run_directory
        print >> handler, '%s/Binary/%s < %s' % (run_directory, executable, config_name)
    else:
        sys.exit("Hostname not recognised. No handler file created.")
    # print >> handler, 'rm -- "$0"'

    try:
        with open('%s' % handler_name, 'w') as handler_file:
            handler_file.write(handler.getvalue())
        print "Handler file written to %s.sh." % filename
    except IOERROR:
        sys.exit("ERROR! Cannot write handler file.")

    subprocess.call("chmod a+x %s.sh" % filename, shell = True)
    print "Submitting batch job."
    if "lxplus" in hostname: subprocess.call('bsub -q %s %s/%s.sh' % (option.queue, run_directory, filename), shell = True)
    elif "cyan03" in hostname: subprocess.call('qsub -l walltime=%s %s/%s.sh' % (option.walltime, run_directory, filename), shell = True)
    else: print "Hostname not recognised. No job submitted."
else:
    if "lxplus" in hostname:
        subprocess.call("source /afs/cern.ch/sw/lcg/external/gcc/4.8/x86_64-slc6/setup.sh", shell = True)
        subprocess.call("export LD_LIBRARY_PATH=/afs/cern.ch/user/d/demillar/.RootTuple:$LD_LIBRARY_PATH", shell = True)
    subprocess.call("./Binary/%s < %s" % (executable, config_name), shell = True)

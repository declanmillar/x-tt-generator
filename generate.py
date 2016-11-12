#!/usr/bin/env python

# python run script for the zprime-top-generator program.
# generates a configuration file, then runs the zprime executable using it, either locally or by submission to the LXPLUS or Iridis batch system.
# author: Declan Millar <declan.millar@cern.ch>

import os, StringIO, re, optparse, subprocess, time, sys, random, glob, socket, string

parser = optparse.OptionParser()
parser.add_option("-t", "--tag", default = "", type = "string", help = "add a name tag to output files")
parser.add_option("-B", "--batch", default = True, action = "store_false", help = "run in batch mode")
parser.add_option("-w", "--walltime", default = "60:00:00", action = "store", help = "walltime 'hh:mm:ss'")
parser.add_option("-Q", "--queue", default = "8nh", action = "store", help = "lxbatch queue'")
parser.add_option("-v", "--verbose", default = False, action = "store_true", help = "run in verbose mode.")

parser.add_option("-T", "--ntuple", default = True, action = "store_false", help = "write events to ROOT ntuple")
parser.add_option("-H", "--lhef", default = False, action = "store_true", help = "write events to lhef file")
parser.add_option("-m", "--model", default = "SM", action = "store", help = "set model")

parser.add_option("-E", "--collider_energy", default = 13, action = "store", help = "collider energy")
parser.add_option("-i", "--initial_state", default = 0, const = 1, action = "store_const", help = "initial state: 0 = pp, 1 = ppbar")
parser.add_option("-f", "--final_state", default = "1", type = "int", action = "store", help = "set final state")
parser.add_option("-g", "--include_gg", default = False, action = "store_true", help = "include gluon-gluon initiated interactions")
parser.add_option("-q", "--include_qq", default = False, action = "store_true", help = "include quark-quark initiated interactions")
parser.add_option("-u", "--include_uu", default = False, action = "store_true", help = "include up-up initiated interactions")
parser.add_option("-d", "--include_dd", default = False, action = "store_true", help = "include down-down initiated interactions")
parser.add_option("-A", "--include_a", default = True, action = "store_false", help = "include photon mediated interaction")
parser.add_option("-Z", "--include_z", default = True, action = "store_false", help = "include Z boson mediated interaction")
parser.add_option("-X", "--include_x", default = True, action = "store_false", help = "include Z' boson mediated interactions")
parser.add_option("-s", "--include_signal", default = True, action = "store_false", help = "include tt signal")
parser.add_option("-b", "--include_background", default = False, action = "store_true", help = "include tt background")

parser.add_option("-P", "--structure_function", default = 4, type = "int", help = "structure_functions")
parser.add_option("-I", "--interference", default = 2, type = "int", help = "specify interference")
parser.add_option("-W", "--use_nwa", default = False, action = "store_true", help = "use Narrow Width Approximation")
parser.add_option("-L", "--ecm_low", default = 0, type = "int", help = " Ecm lower limit")
parser.add_option("-U", "--ecm_up", default = 0, type = "int", help = "Ecm upper limit")
parser.add_option("-F", "--fixed_seed", default = False, action = "store_true", help = "use fixed seed")

parser.add_option("-N", "--vegas_iterations", default = 5, type = "int", help = "number of VEGAS iterations")
parser.add_option("-n", "--vegas_points", default = 10000000, type = "int", help = "number of VEGAS points")

parser.add_option("-S", "--symmetrise", default = True, action = "store_false", help = "symmetrise phase space x1<->x2")
parser.add_option("-R", "--use_rambo", default = False, action = "store_true", help = "use RAMBO for phase space")
parser.add_option("-M", "--map_phase_space", default = True, action = "store_false", help = "map phase space")
parser.add_option("-c", "--cut", default = False, action = "store_true", help = "apply fiducial cuts.")

(option, args) = parser.parse_args()
hostname = socket.gethostname()

if os.path.isfile("Models/%s.mdl" % option.model) is False:
    sys.exit("NOPE! %s is not a valid model.\n Available model files: %s" % (option.model, glob.glob("Models/*.mdl")))
if option.collider_energy < 0:
    sys.exit("NOPE! Collider energy must be positive definite.\n" % usage)
if option.vegas_points < 2:
    sys.exit("NOPE! Must have at least 2 VEGAS points.\n%s" % usage)
if option.ecm_low < 0 or option.ecm_up < 0:
    sys.exit("NOPE! Ecm bounds must be positive definite")
if option.ecm_low > option.collider_energy or option.ecm_up > option.collider_energy:
    sys.exit("NOPE! Ecm bounds cannot exceed collider energy")
if (option.ecm_low > 0 and option.ecm_up > 0 and option.ecm_up <= option.ecm_low):
    sys.exit("NOPE! Ecm upper bound must be greater than lower bound")
if option.interference < 0 or option.interference > 4:
    sys.exit("NOPE! Interference must be from 0-4.")
if option.structure_function < 1 or option.structure_function > 11:
    sys.exit("NOPE! Structure_function ID must be from 1 to 11.")
if option.include_background == False and option.include_signal == False:
    sys.exit("NOPE! Signal and background both off.")
if option.final_state < -1 or option.final_state > 3:

    sys.exit("NOPE! invalid final state id." % option.final_state)
initial_states = 0
if option.include_gg: initial_states += 1
if option.include_qq: initial_states += 1
if option.include_uu: initial_states += 1
if option.include_dd: initial_states += 1
if option.lhef and initial_states > 1:
    sys.exit("NOPE! When outputting to LHEF, only one initial state can be active.")


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
if option.ecm_low != 0 or option.ecm_up != 0:
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
if option.include_gg:
    initial_partons += "gg"
if option.include_qq:
    initial_partons += "qq"
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

filename = '%s_%s-%s%s%s%s_%sx%s' % (option.model, initial_partons, intermediates, final_state, energy_collider, options, option.vegas_iterations, npoints)

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
print >> config, '%r ! ntuple' % option.ntuple
print >> config, '%r ! lhef' % option.lhef
print >> config, '%s' % ntuple_file
print >> config, '%s' % lhe_file
print >> config, '%s' % logfile
print >> config, '%i ! initial_state' % option.initial_state
print >> config, '%i ! final_state' % option.final_state
print >> config, '%s ! model' % option.model
print >> config, '%i ! istructure' % option.structure_function
print >> config, '%r ! include_signal' % option.include_signal
print >> config, '%r ! include_background' % option.include_background
print >> config, '%r ! include_gg' % option.include_gg
print >> config, '%r ! include_qq' % option.include_qq
print >> config, '%r ! include_uu' % option.include_uu
print >> config, '%r ! include_dd' % option.include_dd
print >> config, '%r ! include_a' % option.include_a
print >> config, '%r ! include_z' % option.include_z
print >> config, '%r ! include_x' % option.include_x
print >> config, '%i ! interference' % option.interference
print >> config, '%r ! use_nwa' % option.use_nwa
print >> config, '%i.d3 ! ecm_col' % option.collider_energy
print >> config, '%i ! iseed' % seed
print >> config, '%i ! vegas_iterations' % option.vegas_iterations
print >> config, '%i ! ncall' % option.vegas_points
# print >> config, '-1.d0 ! acc'
print >> config, '%r ! use rambo' % option.use_rambo
print >> config, '%r ! map phase space' % option.map_phase_space
print >> config, '%r ! symmetrise' % option.symmetrise
print >> config, '%r ! verbose mode' % option.verbose
print >> config, '%f ! ecm_low' % (option.ecm_low*1000)
print >> config, '%f ! ecm_up' % (option.ecm_up*1000)
print >> config, '%r ! apply fiducial cuts' % option.cut

try:
    with open('%s' % config_name,'w') as config_file:
        config_file.write(config.getvalue())
        print " config:  %s" % config_name
except IOERROR:
    sys.exit("NOPE! Cannot write to %s" % config_name)

if option.batch:
    handler = StringIO.StringIO()
    if "lxplus" in hostname:
        print >> handler, "export LD_LIBRARY_PATH=/afs/cern.ch/user/d/demillar/.RootTuple:$LD_LIBRARY_PATH"
        # print >> handler, "source /afs/cern.ch/sw/lcg/external/gcc/4.8/x86_64-slc6/setup.sh"
        print >> handler, "cd %s" % run_directory
        print >> handler, '%s/bin/%s < %s' % (run_directory, executable, config_name)
    elif "cyan" in hostname:
        print "walltime = %s" % option.walltime
        print >> handler, "#!/bin/bash"
        print >> handler, "module load gcc/4.8.1; source /local/software/cern/root_v5.34.14/bin/thisroot.sh"
        print >> handler, "export LD_LIBRARY_PATH=/home/dam1g09/.RootTuple:$LD_LIBRARY_PATH"
        print >> handler, "cd %s" % run_directory
        print >> handler, '%s/bin/%s < %s' % (run_directory, executable, config_name)
    else:
        sys.exit("Hostname not recognised. No handler file created.")
    # print >> handler, 'rm -- "$0"'

    try:
        with open('%s' % handler_name, 'w') as handler_file:
            handler_file.write(handler.getvalue())
        print "Handler file written to %s.sh." % filename
    except IOERROR:
        sys.exit("NOPE! Cannot write handler file.")

    subprocess.call("chmod a+x %s.sh" % filename, shell = True)
    print "Submitting batch job."
    if "lxplus" in hostname: subprocess.call('bsub -q %s %s/%s.sh' % (option.queue, run_directory, filename), shell = True)
    elif "cyan03" in hostname: subprocess.call('qsub -l walltime=%s %s/%s.sh' % (option.walltime, run_directory, filename), shell = True)
    else: print "Hostname not recognised. No job submitted."
else:
    if "lxplus" in hostname:
        subprocess.call("source /afs/cern.ch/sw/lcg/external/gcc/4.8/x86_64-slc6/setup.sh", shell = True)
        subprocess.call("export LD_LIBRARY_PATH=/afs/cern.ch/user/d/demillar/.RootTuple:$LD_LIBRARY_PATH", shell = True)
    subprocess.call("./bin/%s < %s" % (executable, config_name), shell = True)

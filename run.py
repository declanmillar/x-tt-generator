#!/usr/bin/env python

# Python run script for the Zprime program.
# Generates a configuration file, then runs the zprime executable using it, either locally or by submission to the LXPLUS or Iridis batch system.
# Do './run.py -h' for help.
# Author: Declan Millar, d.millar@soton.ac.uk.

import os, StringIO, re, optparse, subprocess, time, sys, random, glob, socket

usage = "Usage: ./run.py final_state model_name collider_energy vegas_points [options]"
parser = optparse.OptionParser(usage)

# Execution options
parser.add_option("-L", "--write_logfile", default = False, action = "store_true" , help = "output to logfile")
parser.add_option("-t", "--tag", default = "", type = "string", help = "add a name tag to output files")
parser.add_option("-B", "--batch", default = False, action = "store_true", help = "run in batch mode")
parser.add_option("-W", "--walltime", default = "60:00:00", action = "store", help = "walltime 'hh:mm:ss'")

# Physics options
parser.add_option("-p", "--initial_state", default = 0, const = 1, action = "store_const", help = "switch to p-pbar collisions")
parser.add_option("-S", "--structure_function", default = 4, type = "int", help = "structure_functions set: 1 = CTEQ6M; 2 = CTEQ6D; 3 = CTEQ6L; 4 = CTEQ6L1; ...")
parser.add_option("-C", "--include_qcd", default = True, action = "store_false", help = "turn off QCD")
parser.add_option("-F", "--include_ew", default = True, action = "store_false", help = "turn off EW")
parser.add_option("-Z", "--include_bsm", default = True, action = "store_false", help = "turn off Z'")
parser.add_option("-g", "--include_gg", default = True, action = "store_false", help = "turn off gg")
parser.add_option("-q", "--include_qq", default = True, action = "store_false", help = "turn off qq")
parser.add_option("-i", "--interference", default = 2, type = "int", help = "specify interference: 0 = none, 1 = SM, 2 = full, 3 = full-SM")
parser.add_option("-w", "--use_nwa", default = False, action = "store_true", help = "use NWA")
parser.add_option("-l", "--ecm_low", default = 0, type = "int", help = "ecm lower limit")
parser.add_option("-u", "--ecm_up", default = 0, type = "int", help = "ecm upper limit")

# Monte Carlo options
parser.add_option("-s", "--fixed_seed", default = False, action = "store_true", help = "use fixed seed for random number generator")
parser.add_option("-m", "--itmx", default = 5, type = "int", help = "maximum number of VEGAS iterations")
parser.add_option("-x", "--symmetrise_x1x2", default = False, action = "store_true", help = "symmetrise phase space over x1 and x2")
parser.add_option("-c", "--symmetrise_costheta_t", default = False, action = "store_true", help = "symmetrise phase space over costheta_t")
parser.add_option("-5", "--symmetrise_costheta_5", default = False, action = "store_true", help = "symmetrise phase space over costheta_5")
parser.add_option("-7", "--symmetrise_costheta_7", default = False, action = "store_true", help = "symmetrise phase space over costheta_7")
parser.add_option("-R", "--use_rambo", default = False, action = "store_true", help = "use RAMBO for PS")
parser.add_option("-M", "--map_phase_space", default = True, action = "store_false", help = "flatten Breit-Wigners in integrand for manual phase space")

# Debug option
parser.add_option("-P", "--phase_space_only", default = False, action = "store_true", help = "Set |M|^2  =  1")
parser.add_option("-v", "--verbose", default = False, action = "store_true", help = "Run in verbose mode.")

(option, args) = parser.parse_args()

hostname = socket.gethostname()

if len(args) != 4:
  sys.exit("Error: incorrect number of arguments %i/4.\n%s" % (len(args),usage))

final_state = str(args[0])
model_name = str(args[1])
try:
  collider_energy = int(args[2])
except ValueError:
  sys.exit("Error: invalid collider energy '%s [TeV]'.\n%s" % (args[2],usage))
try:
  ncall = int(args[3])
except ValueError:
  sys.exit("Error: invalid VEGAS points '%s'.\n%s" % (args[3],usage))

# Check arguments are valid
if os.path.isfile("Models/%s.mdl" % model_name) is False:
  sys.exit("'Model/%s.mdl' does not exist.\n%s\nAvailable model files: %s" % (model_name,usage,glob.glob("Models/*.mdl")))

if collider_energy  < 0:
  sys.exit("Error: collider energy must be positive definite.\n" % usage)

if final_state != "ll" and final_state != "tt" and final_state != "bbllnn":
  sys.exit("Error: unavailable final state '%s'.\n%s\nPossible final states: ll, tt, bbllnn" % (final_state,usage))

if ncall < 2:
  sys.exit("Error: Must have at least 2 VEGAS points.\n%s" % usage)

# Check options are valid
if option.batch and "lxplus" in hostname:
  option.write_logfile = False

if option.ecm_low < 0 or option.ecm_up < 0:
  sys.exit("Error: COM energy must be positive definite")

if option.ecm_low > collider_energy or option.ecm_up > collider_energy:
  sys.exit("Error: COM energy cannot exceed collider energy")

if (option.ecm_low > 0 and option.ecm_up > 0 and option.ecm_up <= option.ecm_low):
  sys.exit("Error: E_CM up must be greater than E_CM low")

if option.batch:
  if ("lxplus" not in hostname) or ("iridis" not in hostname)):
    sys.exit("Error: Must be on lxplus or iridis to submit a batch job.")

if option.interference < 0 or option.interference > 4:
  sys.exit("Error: interference must be from 0-4.")

if option.structure_function < 1 or option.structure_function > 9:
  sys.exit("Error: structure_function ID must be from 1 to 9.")

if option.itmx > 20:
   sys.exit('Error: itmx does not have to exceed 20.')

# Modify configuration for consistency
if model_name == "SM":
  option.include_bsm = False

if option.include_bsm is False:
  model_name = "SM"

if final_state == "ll":
  option.include_qcd = False
  option.include_gg = False

if option.phase_space_only:
  option.include_qcd = False
  option.include_ew = False
  option.include_bsm = False

if option.interference == 4 and option.include_ew is False:
  print "EW sector must be active to calculate interference with Zprimes. Switching to default interference."
  option.interference = 2

if option.use_rambo:
  option.map_phase_space = False

if final_state == "ll" or final_state == "tt":
  option.use_nwa = False
  option.symmetrise_costheta_5 = False
  option.symmetrise_costheta_7 = False

# Default iseed
seed = 12345 if option.fixed_seed else random.randint(0,100000)

# Strings
executable = "zprime"

options = ""

if option.initial_state == 1:
  options += "p"

if option.structure_function != 4:
  options += "S%s" % option.structure_function

if option.include_qcd is False and final_state != "ll":
  options += "C"
if option.include_ew is False:
  options += "F"
if option.include_gg is False and final_state != "ll":
  options += "g"
if option.include_qq is False:
  options += "q"

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

if option.ecm_low != 0:
  options += "l%s" % option.ecm_low
if option.ecm_up != 0:
  options += "u%s" % option.ecm_up

# exclude divergence
if final_state == "ll" and option.ecm_low == 0:
  option.ecm_low = 20

if option.fixed_seed:
  options += "s"

# symmetrization
if option.symmetrise_x1x2:
  options += "x"
if option.symmetrise_costheta_t:
  options += "c"
if option.symmetrise_costheta_5:
  options += "5"
if option.symmetrise_costheta_7:
  options += "7"

if option.use_rambo:
  options += "R"
if not option.map_phase_space:
  options += "M"

if len(options) > 0:
  options = "_" + options

if len(option.tag) > 0:
  options = "_" + option.tag

# filename
filename = '%s_%s_%s%s_%sx%s' % (final_state, model_name, collider_energy, options, option.itmx, ncall)

# Generate fortran friendly configuration
if final_state == "ll":
  final_state_id = -1
elif final_state == "tt":
  final_state_id = 0
elif final_state == "bbllnn":
  final_state_id = 1

# Logfile
run_directory = "."
data_directory = "."
if "Lorkhan" in hostname:
  data_directory = "/Users/declan/Data/Zprime"
elif "lxplus" in hostname:
  run_directory = "/afs/cern.ch/user/d/demillar/zprime-top-generator"
  data_directory = "/afs/cern.ch/work/d/demillar/zprime"
elif "iridis" in hostname:
  run_directory = "/home/dam1g09/zprime-top-generator"
  data_directory = "/scratch/dam1g09/zprime"

if os.path.isdir("%s" % data_directory) is False:
  sys.exit("The target data directory '%s' does not exist" % data_directory)

config_name = "%s.cfg" % filename
logfile = "%s.log" % filename
handler_name = "%s.sh" % filename
data_file = '%s/%s.txt' % (data_directory, filename)
ntuple_file = "%s/%s.root" % (data_directory, filename)
logfile_command = "> %s/%s &" % (data_directory, logfile) if option.write_logfile else ""

# Print config file
config = StringIO.StringIO()

print >> config, '%s' % ntuple_file
print >> config, '%s' % data_file
print >> config, '%i ! initial_state' % option.initial_state
print >> config, '%i ! final_state' % final_state_id
print >> config, '%s ! model_name' % model_name
print >> config, '%i ! istructure' % option.structure_function
print >> config, '%i ! include_qcd' % option.include_qcd
print >> config, '%i ! include_ew' % option.include_ew
print >> config, '%i ! include_bsm' % option.include_bsm
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
print >> config, '%i ! symmetrise_x1x2' % option.symmetrise_x1x2
print >> config, '%i ! symmetrise_costheta_t' % option.symmetrise_costheta_t
print >> config, '%i ! symmetrise_costheta_5' % option.symmetrise_costheta_5
print >> config, '%i ! symmetrise_costheta_7' % option.symmetrise_costheta_7
print >> config, '%i ! verbose mode' % option.verbose
print >> config, '%i.d3 ! ecm_low' % option.ecm_low
print >> config, '%i.d3 ! ecm_up' % option.ecm_up

try:
  with open('Config/%s' % config_name,'w') as config_file:
    config_file.write(config.getvalue())
    print " Config: %s." % config_name
except IOError:
  sys.exit(" Error: cannot write Config/%s. Are you running in the right directory?" % config_name)

if option.batch:
  handler = StringIO.StringIO()
  if "lxplus" in hostname:
    print >> handler, "export LD_LIBRARY_PATH=/afs/cern.ch/user/d/demillar/.RootTuple:$LD_LIBRARY_PATH"
    print >> handler, "source /afs/cern.ch/sw/lcg/external/gcc/4.8/x86_64-slc6/setup.sh"
    print >> handler, "cd /afs/cern.ch/user/d/demillar/Zp-tt_pheno/Generation/"
    print >> handler, '%s/Binary/%s < %s/Config/%s' % (run_directory, executable, run_directory, config_name)
    print >> handler, 'mv LSFJOB_* Jobs'
  if "iridis" in hostname:
    print "walltime = %s" % option.walltime
    print >> handler, "#!/bin/bash"
    print >> handler, "module load gcc/4.8.1; source /local/software/cern/root_v5.34.14/bin/thisroot.sh"
    print >> handler, "export LD_LIBRARY_PATH=/home/dam1g09/.RootTuple:$LD_LIBRARY_PATH"
    print >> handler, "cd %s" % run_directory
    # print >> handler, "cd $PBS_O_WORKDIR"
    print >> handler, '%s/Binary/%s < %s/Config/%s %s' % (run_directory, executable, run_directory, config_name, logfile_command)
  print >> handler, 'rm -- "$0"'
  try:
    with open('%s' % handler_name, 'w') as handler_file:
      handler_file.write(handler.getvalue())
    print " Handler file written to %s.sh." % filename
  except IOError:
    sys.exit(" Error: cannot write handler file.")

  subprocess.call("chmod a+x %s.sh" % filename, shell = True)
  print " Submitting batch job."
  if "lxplus" in hostname: subprocess.call('bsub -q 1nh /afs/cern.ch/user/d/demillar/Zp-tt_pheno/Generation/%s.sh' % filename, shell = True)
  if "iridis" in hostname: subprocess.call('qsub -l walltime=%s %s/%s.sh' % (option.walltime, run_directory, filename), shell = True)

else:
  if "lxplus" in hostname:
    print " Sourcing ROOT..."
    subprocess.call("source /afs/cern.ch/sw/lcg/external/gcc/4.8/x86_64-slc6/setup.sh", shell = True)
    print " Adding RootTuple libraries to library path..."
    subprocess.call("export LD_LIBRARY_PATH=/afs/cern.ch/user/d/demillar/.RootTuple:$LD_LIBRARY_PATH", shell = True)
  print " Executing locally."
  if option.write_logfile:
    print " Logfile: /%s" % logfile
  subprocess.call("./Binary/%s < Config/%s %s" % (executable, config_name, logfile_command), shell = True)

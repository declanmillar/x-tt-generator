#!/usr/bin/env python 
# ------------------------------------------------------------------------------------------
# Command line executor for Zp-tt code.
# Generates a .com file then runs the zprime executable using it.
# Do ./run.py -h for help.

# Arguments:
#   0 = mdl:   name of model_name file in Models directory
#   2 = ecm:   centre of mass energy (input:TeV, output: GeV)
#   3 = final_state: final_state state (0: no decay, 1: dilepton)
#   4 = ncall: number of vegas calls (monte carlo/phase space points)
# ------------------------------------------------------------------------------------------

import os, StringIO, re, optparse, subprocess, time, sys, random

# Usage help
parser = optparse.OptionParser("usage: ./%prog [options] model_name ecm_col o_final ncall")

# Execution options
parser.add_option("-l", "--output", default=False, action="store_true" , help="output to logfile")
parser.add_option("-t", "--tag", default="", type="string", help="add a name tag to logfile")
parser.add_option("-b", "--batch", default=False, action="store_true", help = "run in batch mode")

# Physics options
parser.add_option("-p", "--initial_state", default = 0, const = 1, action = "store_const", help = "switch to p-pbar collisions")
parser.add_option("-S", "--structure_function", default = 4, type = "int", help = "structure_functions set: 1 = CTEQ6M; 2 = CTEQ6D; 3 = CTEQ6L; 4 = CTEQ6L1; ...")
parser.add_option("-C", "--include_qcd", default = 1, const = 0, action = "store_const", help = "turn off QCD")
parser.add_option("-F", "--include_ew", default = 1, const = 0, action = "store_const", help = "turn off EW")
parser.add_option("-Z", "--include_bsm", default = 1, const = 0, action = "store_const", help = "turn off Z'")
parser.add_option("-g", "--include_gg", default = 1, const = 0, action = "store_const", help = "turn off gg")
parser.add_option("-q", "--include_qq", default = 1, const = 0, action = "store_const", help = "turn off qq")
parser.add_option("-i", "--interference", default = 2, type = "int", help = "specify interference: 0 = none, 1 = SM, 2 = full, 3 = full-SM")
parser.add_option("-w", "--use_nwa", default = 0, const = 1, action = "store_const", help = "turn on use_nwa")

# Monte Carlo options
parser.add_option("-s", "--iseed", default = False, action = "store_true", help = "used fixed iseed for random number generator")
parser.add_option("-m", "--itmx", default = 5, type = "int", help = "maximum number of VEGAS iterations")
parser.add_option("-x", "--symmetrise_x1x2", default = 0, const = 1, action = "store_const", help = "symmatrise phase space over x1 and x2")
parser.add_option("-c", "--symmetrise_costheta_t", default = 0, const = 1, action = "store_const", help = "symmatrise phase space over costheta_t")
parser.add_option("-5", "--symmetrise_costheta_5", default = 0, const = 1, action = "store_const", help = "symmatrise phase space over costheta_5")
parser.add_option("-7", "--symmetrise_costheta_7", default = 0, const = 1, action = "store_const", help = "symmatrise phase space over costheta_7")
parser.add_option("-R", "--use_rambo", default = 0, const = 1, action = "store_const", help = "Use RAMBO for PS. Default is manual.")
parser.add_option("-M", "--map_phase_space", default = 1, const = 0, action = "store_const", help = "Flatten Breit-Wigners in integrand for manual phase space.")

# Debug options
parser.add_option("-P", "--phase_space_only", default = 0, const = 1, action = "store_const", help = "Set |M|^2  =  1")
parser.add_option("-v", "--verbose", default = 0, const = 1, action = "store_const", help = "Run in verbose mode.")

(options, args) = parser.parse_args()
print "\n Generating config file..."
  
# Collect arguments
model_name = args[0]
collider_energy = args[1]
final_state = args[2]
ncall = args[3]

if options.phase_space_only == 1:
  options.include_qcd = 0
  options.include_ew = 0
  options.include_bsm = 0

if options.include_bsm == 1:
  smodel = args[0]
elif options.include_bsm == 0:
  smodel = ""

# Default iseed
seed = 12345 if options.iseed else random.randint(0,100000)

# Strings
executable="zprime"
sfinal = "2to2" if (final_state == "0") else "2to6"
config = StringIO.StringIO()
handler = StringIO.StringIO()

# Gauge sectors
if   options.include_qcd == 1 and options.include_ew == 1 and options.include_bsm == 1:
  sector = ""
elif options.include_qcd == 1 and options.include_ew == 1 and options.include_bsm == 0:
  sector = "SM"
elif options.include_qcd == 1 and options.include_ew == 0 and options.include_bsm == 1:
  sector = "_QCD-Zp"  
elif options.include_qcd == 0 and options.include_ew == 1 and options.include_bsm == 1:
  sector = "_EW-Zp" 
elif options.include_qcd == 1 and options.include_ew == 0 and options.include_bsm == 0:
  sector = "QCD"    
elif options.include_qcd == 0 and options.include_ew == 1 and options.include_bsm == 0:
  sector = "EW"
elif options.include_qcd == 0 and options.include_ew == 0 and options.include_bsm == 1:
  sector = "_Zp"
elif options.include_qcd == 0 and options.include_ew == 0 and options.include_bsm == 0:
  sector = "PS"

all_options = ""

# Interference
if (options.interference == 4) and (options.include_ew == 0):
  options.interference = 2
  print 'EW sector must be active to calculate interference with Zprimes.'
  print 'Switching to default interference.'

if options.interference == 0:
  all_options += "i0"
elif options.interference == 1:
  all_options += "i1"
elif options.interference == 3:
  all_options += "i3"
elif options.interference == 4:
  all_options += "i4"

# symmetrization
if options.symmetrise_x1x2 == 1:
  all_options += "x"
if options.symmetrise_costheta_t == 1:
  all_options += "c"
if options.symmetrise_costheta_5 == 1:
  all_options += "5"
if options.symmetrise_costheta_7 == 1:
  all_options += "7"

if options.use_rambo == 1:
  all_options += "R"
if options.map_phase_space == 0:
  all_options += "M"

if options.include_gg == 0:
  all_options += "g"
if options.include_qq == 0:
  all_options += "q"

if len(all_options) > 0:
  all_options = "_" + all_options

# filename
filename = '%s_%s%s_%s%s%s_%sx%s' % (sfinal, smodel, sector, collider_energy, all_options, options.tag, options.itmx, ncall)

# logfile 
os = sys.platform
if (os == "darwin"):
  ntuple_directory = "/Users/declan/Data/Ntuples_Zprime/"
elif (os == "linux2"):
  if options.batch == False:
    subprocess.call("export LD_LIBRARY_PATH=/afs/cern.ch/user/d/demillar/.RootTuple:$LD_LIBRARY_PATH", shell=True)
    subprocess.call("source /afs/cern.ch/sw/lcg/external/gcc/4.8/x86_64-slc6/setup.sh", shell=True)
  ntuple_directory = "/afs/cern.ch/work/d/demillar/Ntuples_Zprime/"
  
logfile = '> Logs/%s.log &' % (filename) if options.output else ''
output_file = "%s.out" % filename
ntuple_file = ntuple_directory + '%s.root' % filename

# print config file
print >> config, '%s' % ntuple_file

print >> config, '%s ! output file' % output_file

print >> config, '%s ! initial_state' % options.initial_state

print >> config, '%s ! final_state' % final_state

print >> config, '%s ! model_name' % model_name

print >> config, '%s ! istructure' % options.structure_function

print >> config, '%s ! include_qcd' % options.include_qcd

print >> config, '%s ! include_ew' % options.include_ew

print >> config, '%s ! include_bsm' % options.include_bsm

print >> config, '%s ! include_gg' % options.include_gg

print >> config, '%s ! include_qq' % options.include_qq

print >> config, '%s ! phase_space_only' % options.phase_space_only

print >> config, '%s ! interference' % options.interference

print >> config, '%s ! use_nwa' % options.use_nwa

print >> config, '%s.d3 ! ecm_col' % collider_energy

print >> config, '%s ! iseed' % seed

print >> config, '%s ! itmx' % options.itmx

print >> config, '%s ! ncall' % ncall

print >> config, '-1.d0 ! acc'

print >> config, '%s ! use rambo' % options.use_rambo

print >> config, '%s ! map phase space' % options.map_phase_space

print >> config, '%s ! symmetrise_x1x2' % options.symmetrise_x1x2

print >> config, '%s ! symmetrise_costheta_t' % options.symmetrise_costheta_t

print >> config, '%s ! symmetrise_costheta_5' % options.symmetrise_costheta_5

print >> config, '%s ! symmetrise_costheta_7' % options.symmetrise_costheta_7

print >> config, '%s ! verbose mode' % options.verbose

try:
  with open('Config/%s.com' % filename,'w') as cfile1:
    cfile1.write(config.getvalue())
except IOError:
  print "Not in right directory?"
  sys.exit()

if options.batch == True:
  print >> handler, "export LD_LIBRARY_PATH=/afs/cern.ch/user/d/demillar/.RootTuple:$LD_LIBRARY_PATH"
  print >> handler, "source /afs/cern.ch/sw/lcg/external/gcc/4.8/x86_64-slc6/setup.sh"
  print >> handler, "cd /afs/cern.ch/user/d/demillar/Zp-tt_pheno/Generation/"
  print >> handler, '/afs/cern.ch/user/d/demillar/Zp-tt_pheno/Generation/Binary/%s < /afs/cern.ch/user/d/demillar/Zp-tt_pheno/Generation/Config/%s.com %s' % (executable,filename,logfile)

  try:
    with open('%s.sh' % filename,'w') as hfile:
      hfile.write(handler.getvalue())
  except IOError:
    print "Not in right directory?"
    sys.exit()
else:
  subprocess.call("export LD_LIBRARY_PATH=/afs/cern.ch/user/d/demillar/.RootTuple:$LD_LIBRARY_PATH", shell=True)
  subprocess.call("source /afs/cern.ch/sw/lcg/external/gcc/4.8/x86_64-slc6/setup.sh", shell=True)

# Command
if options.batch == True:
 permission = "chmod a+x %s.sh" % filename
 subprocess.call(permission, shell=True)
 command = 'bsub -q 1nw /afs/cern.ch/user/d/demillar/Zp-tt_pheno/Generation/%s.sh' % filename
else:
 command = './Binary/%s < Config/%s.com %s' % (executable,filename,logfile)
print " ...complete."
print ' Executing ', command, '...'
subprocess.call(command, shell=True)

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

# Local options
parser.add_option("-l", "--output", default=False, action="store_true" , help="output to logfile")
parser.add_option("-t", "--tag", default="",  type="string", help="add a name tag to logfile")

# Physics options
parser.add_option("-p", "--initial_state", default=0, const=1, action="store_const", help="switch to p-pbar collisions")
parser.add_option("-P", "--structure_function", default=4, type="int", help="structure_functions set: 1=CTEQ6M; 2=CTEQ6D; 3=CTEQ6L; 4=CTEQ6L1; ...")
parser.add_option("-q", "--include_qcd", default=1, const=0, action="store_const", help="turn off QCD")
parser.add_option("-e", "--include_ew", default=1, const=0, action="store_const", help="turn off EW")
parser.add_option("-z", "--include_bsm", default=1, const=0, action="store_const", help="turn off Z'")
parser.add_option("-i", "--interference", default=2, type="int", help="specify interference: 0=none, 1=SM, 2=full, 3=full-SM")
parser.add_option("-w", "--use_NWA", default=0, const=1, action="store_const", help="turn on use_NWA")
parser.add_option("-K", "--additional_kinematics", default=1, const=0, action="store_const", help="Turn off additional kinematics")
parser.add_option("-T", "--include_transverse", default=1, const=0, action="store_const", help="switch off include transverse mass variables")
parser.add_option("-A", "--include_asymmetries", default=1, const=0, action="store_const", help="switch off asymmetry variables")

# Monte Carlo options
parser.add_option("-s", "--iseed", default=False, action="store_true", help="used fixed iseed for random number generator")
parser.add_option("-m", "--itmx", default=1, type="int", help="maximum number of VEGAS iterations")
parser.add_option("-x", "--symmetrise_x1x2", default=0, const=1, action="store_const", help="symmatrise phase space over x1 and x2")
parser.add_option("-c", "--symmetrise_costheta_t", default=0, const=1, action="store_const", help="symmatrise phase space over costheta_t")
parser.add_option("-D", "--print_distributions", default=0, const=1, action="store_const", help="turn on built in histogams")
parser.add_option("-E", "--include_errors", default=0, const=1, action="store_const", help="turn on distribution errors")

# Debug options
parser.add_option("-M", "--phase_space_only", default=0, const=1, action="store_const", help="Set |M|^2 = 1")
parser.add_option("-v", "--verbose", default=0, const=1, action="store_const", help="Run in verbose mode.")

(options, args) = parser.parse_args()
print "\n Generating config file..."

if options.phase_space_only==1:
  options.include_qcd=0
  options.include_ew=0
  options.include_bsm=0
  
# Collect arguments
model_name = args[0]

if   options.include_bsm == 1:
  smodel = args[0]
elif options.include_bsm == 0:
  smodel = ""

collider_energy = args[1]
final_state = args[2]
ncall = args[3]

# Default iseed
seed = 12345 if options.iseed else random.randint(0,100000)

# Strings
executable="Binary/"+"zprime"
sfinal = "2to2" if (final_state=="0") else "2to6"
config=StringIO.StringIO()

if (final_state == "1"):
  sfinal += "_DL"
elif (final_state == "2"):
  sfinal += "_SL"
elif (final_state == "3"):
  sfinal += "_FH"

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

if len(all_options) > 0:
  all_options = "_" + all_options

# filename
filename = '%s%s_%s_%s%s%s_%sx%s' % (smodel, sector, collider_energy, sfinal, all_options, options.tag, options.itmx, ncall)

# logfile 
os = sys.platform
if (os == "darwin"):
  ntuple_directory = "/Users/declan/Data/Ntuples_Zprime/"
elif (os == "linux2"):
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

print >> config, '%s ! phase_space_only' % options.phase_space_only

print >> config, '%s ! interference' % options.interference

print >> config, '%s ! use_NWA' % options.use_NWA

print >> config, '%s ! additional_kinematics' % options.additional_kinematics

print >> config, '%s.d3 ! ecm_col' % collider_energy

print >> config, '%s ! iseed' % seed

print >> config, '%s ! itmx' % options.itmx

print >> config, '%s ! ncall' % ncall

print >> config, '-1.d0 ! acc'

print >> config, '%s ! symmetrise_x1x2' % options.symmetrise_x1x2

print >> config, '%s ! symmetrise_costheta_t' % options.symmetrise_costheta_t

print >> config, '%s ! print_distributions' % options.print_distributions

print >> config, '%s ! include distribution errors' % options.include_errors

print >> config, '%s ! verbose mode' % options.verbose

try:
      with open('Config/%s.com' % filename,'w') as cfile1:
            cfile1.write(config.getvalue())
except IOError:
      print "Not in right directory?"
      sys.exit()

# Command
command = './%s < Config/%s.com %s' % (executable,filename,logfile)
print " ...complete."
print ' Executing ', command, '...'
subprocess.call(command, shell=True)

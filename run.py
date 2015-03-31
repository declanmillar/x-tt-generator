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

import os,StringIO,re,optparse,subprocess,time,sys,random

# Usage help
parser = optparse.OptionParser("usage: ./%prog [options] model_name ecm_col o_final ncall")

# Local options
parser.add_option("-l", "--output" , default=False, action="store_true"               , help="output to logfile")
parser.add_option("-t", "--tag"    , default=""   ,   type="string"                   , help="add a name tag to logfile")

# Physics options
parser.add_option("-p", "--initial_state"   ,  default=0   , const=1      , action="store_const", help="switch to p-pbar collisions")
parser.add_option("-P", "--structure_function"   ,  default=4   ,  type="int"                        , help="structure_functions set: 1=CTEQ6M; 2=CTEQ6D; 3=CTEQ6L; 4=CTEQ6L1; ...")
parser.add_option("-q", "--include_qcd"   ,  default=1   , const=0      , action="store_const", help="turn off QCD")
parser.add_option("-e", "--include_ew"    ,  default=1   , const=0      , action="store_const", help="turn off EW")
parser.add_option("-z", "--include_bsm"   ,  default=1   , const=0      , action="store_const", help="turn off Z'")
parser.add_option("-i", "--interference"   ,  default=2   ,  type="int"                        , help="specify interference: 0=none, 1=SM, 2=full, 3=full-SM")
parser.add_option("-w", "--use_NWA"   ,  default=0   , const=1      , action="store_const", help="turn on use_NWA")
parser.add_option("-B", "--use_branching_ratio"    ,  default=0   , const=1      , action="store_const", help="multiply 2to2 process by use_branching_ratio")
parser.add_option("-T", "--include_transverse"  ,  default=1   , const=0      , action="store_const", help="switch off include_transversesverse mass variables")
parser.add_option("-A", "--include_asymmetries"  ,  default=1   , const=0      , action="store_const", help="switch off asymmetry variables")
parser.add_option("-C", "--cuts"  ,  default=1   , const=0      , action="store_const", help="turn off all cuts")
parser.add_option("-y", "--ytmax" ,  default=100 ,  type="float",                       help="rapidity cut")
parser.add_option("-Y", "--yttmin",  default=0   ,  type="float",                       help="top pair rapidity cut")

# Monte Carlo options
parser.add_option("-s", "--iseed"  , default=False,               action="store_true" , help="used fixed iseed for random number generator")
parser.add_option("-m", "--itmx"   , default=5    ,  type="int" ,                       help="maximum number of VEGAS iterations")
parser.add_option("-x", "--symmetries_x1x2", default=0    , const=1     , action="store_const", help="symmatrise phase space over x1 and x2")
parser.add_option("-c", "--symmetrise_costheta_t", default=0    , const=1     , action="store_const", help="symmatrise phase space over costheta_t")
parser.add_option("-D", "--print_all_distributions", default=1    , const=0     , action="store_const", help="turn off distributions")
parser.add_option("-2", "--print_2d_distributions" , default=1    , const=0     , action="store_const", help="turn off 2d-distributions")

# Debug options
parser.add_option("-M", "--phase_space_only",  default=0    , const=1     , action="store_const", help="Set |M|^2 = 1")

(options, args) = parser.parse_args()

if options.phase_space_only==1:
  options.include_qcd=0
  options.include_ew=0
  options.include_bsm=0
  
# Collect arguments
model_name = args[0]
if   options.include_bsm==1:
  smodel = args[0]
elif options.include_bsm==0:
  smodel = ""
collider_energy=args[1]
final_state=args[2]
ncall=args[3]

# Default iseed
seed = 12345 if options.iseed else random.randint(0,100000)

# Strings
executable="Binary/"+"zprime"
sfinal = "2to2" if (final_state=="0") else "2to6"
config=StringIO.StringIO()

# Gauge sectors
if   options.include_qcd==1 and options.include_ew==1 and options.include_bsm==1:
  sector = ""
elif options.include_qcd==1 and options.include_ew==1 and options.include_bsm==0:
  sector = "SM"
elif options.include_qcd==1 and options.include_ew==0 and options.include_bsm==1:
  sector = "_QCD-Zp"  
elif options.include_qcd==0 and options.include_ew==1 and options.include_bsm==1:
  sector = "_EW-Zp" 
elif options.include_qcd==1 and options.include_ew==0 and options.include_bsm==0:
  sector = "QCD"    
elif options.include_qcd==0 and options.include_ew==1 and options.include_bsm==0:
  sector = "EW"
elif options.include_qcd==0 and options.include_ew==0 and options.include_bsm==1:
  sector = "_Zp"
elif options.include_qcd==0 and options.include_ew==0 and options.include_bsm==0:
  sector = "PS"

# Interference  
if   options.interference==0:
  interference="_int0"
elif options.interference==1:
  interference="_int1"
elif options.interference==2:
  interference=""
elif options.interference==3:
  interference="_int3"

# Symmetrization
if   (options.symmetries_x1x2==1) and (options.symmetrise_costheta_t==0):
  symmetrization = "_x"
elif (options.symmetries_x1x2==0) and (options.symmetrise_costheta_t==1):
  symmetrization = "_c"
elif (options.symmetries_x1x2==1) and (options.symmetrise_costheta_t==1):
  symmetrization = "_xc"
else:
  symmetrization = ""

# Filename
filename = '%s%s_%s%s_%s%s%s_%sx%s' % (smodel,sector,sfinal,interference,collider_energy,symmetrization,options.tag,options.itmx,ncall)

# Logfile      
logfile = '> Output/%s.log &' % (filename) if options.output else ''
output_file = '%s.out' % filename

# Print config file

print >> config, '%s  ! output file' % output_file

print >> config, '%s ! initial_state' % options.initial_state

print >> config, '%s ! final_state' % final_state

print >> config, '%s ! model_name' % model_name

print >> config, '%s ! ISTRUCTURE' % options.structure_function

print >> config, '%s ! include_qcd' % options.include_qcd

print >> config, '%s ! include_ew' % options.include_ew

print >> config, '%s ! include_bsm' % options.include_bsm

print >> config, '%s ! phase_space_only' % options.phase_space_only

print >> config, '%s ! interference' % options.interference

print >> config, '%s ! use_branching_ratio' % options.use_branching_ratio

print >> config, '%s ! use_NWA' % options.use_NWA

print >> config, '%s ! include_transverse' % options.include_transverse

print >> config, '%s ! include_asymmetries' % options.include_asymmetries

print >> config, '%s.d3 ! ecm_col' % collider_energy

print >> config, '%sd0 ! ytmax' % options.ytmax

print >> config, '%sd0 ! yttmin' % options.yttmin

print >> config, '%s ! iseed' % seed

print >> config, '%s ! itmx' % options.itmx

print >> config, '%s ! ncall' % ncall

print >> config, '-1.d0 ! acc'

print >> config, '%s ! symmetries_x1x2' % options.symmetries_x1x2

print >> config, '%s ! symmetrise_costheta_t' % options.symmetrise_costheta_t

print >> config, '%s ! print_all_distributions' % options.print_all_distributions

print >> config, '%s ! print_2d_distributions' % options.print_2d_distributions

try:
      with open('Config/%s.com' % filename,'w') as cfile1:
            cfile1.write(config.getvalue())
except IOError:
      print "Not in right directory?"
      sys.exit()

# Command
command = './%s < Config/%s.com %s' % (executable,filename,logfile)
print 'executing: ',command
subprocess.call(command, shell=True)
   


#!/usr/bin/env python 
# ------------------------------------------------------------------------------------------
# Command line executor for Zp-tt code.
# Generates a .com file then runs the ttbar_BSM executable using it.
# Do ./run.py -h for help.

# Arguments:
#   0 = mdl:   name of model file in Models directory
#   2 = ecm:   centre of mass energy (input:TeV, output: GeV)
#   3 = final: final state (0: no decay, 1: dilepton)
#   4 = ncall: number of vegas calls (monte carlo/phase space points)
# ------------------------------------------------------------------------------------------

import os,StringIO,re,optparse,subprocess,time,sys,random

# Usage help
parser = optparse.OptionParser("usage: ./%prog [options] model ecm_col o_final ncall")

# Local options
parser.add_option("-o", "--output" , default=False, action="store_true"               , help="output to terminal")
parser.add_option("-t", "--tag"    , default=""   ,   type="string"                   , help="add a name tag to logfile")

# Physics options
parser.add_option("-p", "--col"   ,  default=0   , const=1      , action="store_const", help="switch to p-pbar collisions")
parser.add_option("-P", "--pdf"   ,  default=4   ,  type="int"                        , help="PDF set: 1=CTEQ6M; 2=CTEQ6D; 3=CTEQ6L; 4=CTEQ6L1; ...")
parser.add_option("-q", "--qcd"   ,  default=1   , const=0      , action="store_const", help="turn off QCD")
parser.add_option("-e", "--ew"    ,  default=1   , const=0      , action="store_const", help="turn off EW")
parser.add_option("-z", "--bsm"   ,  default=1   , const=0      , action="store_const", help="turn off Z'")
parser.add_option("-i", "--int"   ,  default=2   ,  type="int"                        , help="specify interference: 0=none, 1=SM, 2=full, 3=full-SM")
parser.add_option("-n", "--NWA"   ,  default=0   , const=1      , action="store_const", help="turn on NWA")
parser.add_option("-r", "--BR"    ,  default=0   , const=1      , action="store_const", help="multiply 2to2 process by BR")
parser.add_option("-T", "--tran"  ,  default=1   , const=0      , action="store_const", help="switch off transverse mass variables")
parser.add_option("-A", "--asym"  ,  default=1   , const=0      , action="store_const", help="switch off asymmetry variables")
parser.add_option("-C", "--cuts"  ,  default=1   , const=0      , action="store_const", help="turn off all cuts")
parser.add_option("-y", "--ytmax" ,  default=100 ,  type="float",                       help="rapidity cut")
parser.add_option("-Y", "--yttmin",  default=0   ,  type="float",                       help="top pair rapidity cut")

# Monte Carlo options
parser.add_option("-s", "--iseed"  , default=False,               action="store_true" , help="used fixed iseed for random number generator")
parser.add_option("-m", "--itmx"   , default=5    ,  type="int" ,                       help="maximum number of VEGAS iterations")
parser.add_option("-x", "--symx1x2", default=0    , const=1     , action="store_const", help="symmatrise phase space over x1 and x2")
parser.add_option("-c", "--symcost", default=0    , const=1     , action="store_const", help="symmatrise phase space over costheta_t")
parser.add_option("-D", "--distros", default=1    , const=0     , action="store_const", help="turn off distributions")

# Debug options
parser.add_option("-M", "--M_eq_1",  default=0    , const=1     , action="store_const", help="Set |M|^2 = 1")

(options, args) = parser.parse_args()

if options.M_eq_1==1:
  options.qcd=0
  options.ew=0
  options.bsm=0
  
# Collect arguments
model = args[0]
if   options.bsm==1:
  smodel = args[0]
elif options.bsm==0:
  smodel = ""
emc_col=args[1]
final=args[2]
ncall=args[3]

# Default iseed
seed = 12345 if options.iseed else random.randint(0,100000)

# Strings
executable="ttbar_BSM"
sfinal = "2to2" if (final=="0") else "2to6"
config=StringIO.StringIO()

# Gauge sectors
if   options.qcd==1 and options.ew==1 and options.bsm==1:
  sector = ""
elif options.qcd==1 and options.ew==1 and options.bsm==0:
  sector = "SM"
elif options.qcd==1 and options.ew==0 and options.bsm==1:
  sector = "_QCD-Zp"  
elif options.qcd==0 and options.ew==1 and options.bsm==1:
  sector = "_EW-Zp" 
elif options.qcd==1 and options.ew==0 and options.bsm==0:
  sector = "QCD"    
elif options.qcd==0 and options.ew==1 and options.bsm==0:
  sector = "EW"
elif options.qcd==0 and options.ew==0 and options.bsm==1:
  sector = "_Zp"
elif options.qcd==0 and options.ew==0 and options.bsm==0:
  sector = "PS"

# Interference  
if   options.int==0:
  interference="_int0"
elif options.int==1:
  interference="_int1"
elif options.int==2:
  interference=""
elif options.int==3:
  interference="_int3"

# Symmetrization
if   (options.symx1x2==1) and (options.symcost==0):
  symmetrization = "_symx"
elif (options.symx1x2==0) and (options.symcost==1):
  symmetrization = "_symc"
elif (options.symx1x2==1) and (options.symcost==1):
  symmetrization = "_symxc"
else:
  symmetrization = ""

# Print config file
print >> config, '%s ! col' % options.col
print >> config, '%s.d3 ! ecm_col' % emc_col
print >> config, '%s ! ISTRUCTURE' % options.pdf
print >> config, '%s ! model' % model
print >> config, '%s ! QCD' % options.qcd
print >> config, '%s ! EW' % options.ew
print >> config, '%s ! BSM' % options.bsm
print >> config, '%s ! int' % options.int
print >> config, '%s ! sfinal' % final
print >> config, '%s ! NWA' % options.NWA
print >> config, '%s ! BR' % options.BR
print >> config, '%s ! tran' % options.tran
print >> config, '%s ! asym' % options.asym
print >> config, '%s.d0 ! ytmax' % options.ytmax
print >> config, '%s.d0 ! yttmin' % options.yttmin

print >> config, '%s ! iseed' % seed
print >> config, '%s ! itmx' % options.itmx
print >> config, '%s ! ncall' % ncall
print >> config, '-1.d0 ! acc'
print >> config, '%s ! symx1x2' % options.symx1x2
print >> config, '%s ! symcost' % options.symcost
print >> config, '%s ! distros' % options.distros

print >> config, '%s ! M_eq_1' % options.M_eq_1


# Filename
filename = '%s_%s%s%s_%s_%sx%s%s%s' % (sfinal,smodel,sector,interference,emc_col,options.itmx,ncall,symmetrization,options.tag)
try:
      with open('Config/%s.com' % filename,'w') as cfile1:
            cfile1.write(config.getvalue())
except IOError:
      print "Not in right directory?"
      sys.exit()

# Logfile      
write = '' if options.output else '> Output/%s.log &' % (filename)

# Command
command = './%s < Config/%s.com %s' % (executable,filename,write)
print 'executing: ',command
subprocess.call(command, shell=True)
   


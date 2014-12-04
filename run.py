#!/usr/bin/env python 
# ------------------------------------------------------------------------------------------
# Command line executor for Zp-tt code. Use run.py -h for help.
# Generates a .com file to feed into fortran Zp-tt executable and runs it.

# Arguments:
# 	0 = executable e.g. ttbar_BSM. see makefile
#   1 = mdl: Name of model file in Models directory
# 	2 = ecm: centre of mass energy input:TeV, output: GeV
#   3 = final: final state
# 	4 = ncall: number of vegas calls (monte carlo points)
# ------------------------------------------------------------------------------------------

import os,StringIO,re,optparse,subprocess,time,sys,random

# Options
parser = optparse.OptionParser("usage: %prog [options] executable model ecm_col ifs ncall")

# External options
parser.add_option("-o", "--output", action="store_true", default=False, help="output to terminal")
parser.add_option("-t", "--tag", default="",type="string",help="add specific tag to logfile name")

# Internal options
parser.add_option("-q", "--iqcd",   default=1,action="store_const",const=0, help="Turn off QCD contribution for ttbar")
parser.add_option("-e", "--iew",    default=1,action="store_const",const=0, help="Turn off EW (gamma/Z) contribution")
parser.add_option("-b", "--ibsm",   default=1,action="store_const",const=0, help="Turn off Z' contribution")
parser.add_option("-f", "--iint",   type="int",default=2, help="Specify interference: 0=none; 1=SM interference only; 2=all interference; 3=Zp interference only.")
parser.add_option("-n", "--iNWA",   action="store_const", const=1, default=0, help="NWA")
parser.add_option("-c", "--icol",   action="store_const", const=1, default=0, help="collider")
parser.add_option("-P", "--pdf",    default=4,type="int", help="PDF set: 1=CTEQ6M; 2=CTEQ6D; 3=CTEQ6L; 4=CTEQ6L1; ...")
parser.add_option("-y", "--ytcut",  default=100,type="float", help="rapidity cut")
parser.add_option("-Y", "--yttcut", default=0,type="float", help="top pair rapidity cut")
parser.add_option("-s", "--iseed",  action="store_true", default=False, help="used fixed seed for random number generator")
parser.add_option("-D", "--idist",  action="store_const", const=0, default=1, help="Do not make standard distributions.")
parser.add_option("-T", "--itdist",  action="store_const", const=0, default=1, help="Do not make transverse mass distributions.")
parser.add_option("-A", "--iadist",  action="store_const", const=0, default=1, help="Do not make asymmetry distributions.")
parser.add_option("-m", "--itmx",   type="int",default=5, help="Maximum number of VEGAS iterations.")
parser.add_option("-l", "--ilhe",   default=0,action="store_const",const=1, help="Output in lhe format.")

(options, args) = parser.parse_args()

# Collect arguments
model=args[0]
emc_col=args[1]
ifs=args[2]
ncall=args[3]

seed = 12345 if options.iseed else random.randint(0,100000)

# Strings
executable="ttbar_BSM"
name = "2to2" if (ifs=="0") else "2to6"
config=StringIO.StringIO()
sector = "_EW" if (options.iqcd==0) else ""
if   options.iint==0:
	interference="_int0"
elif options.iint==1:
	interference="_int1"
elif options.iint==2:
	interference=""
elif options.iint==3:
	interference="_int3"

# Printing
print >> config, '%s ! iQCD' % options.iqcd
print >> config, '%s ! iEW' % options.iew
print >> config, '%s ! iBSM' % options.ibsm
print >> config, '%s ! iint' % options.iint
print >> config, '%s ! ifs' % ifs
print >> config, '%s ! iNWA' % options.iNWA
print >> config, '%s ! model' % model
print >> config, '%s.d3 ! ecm_col' % emc_col
print >> config, '%s ! icol' % options.icol
print >> config, '%s ! ISTRUCTURE' % options.pdf
print >> config, '%s.d0 ! ytcut' % options.ytcut
print >> config, '%sd0 ! yttcut' % options.yttcut
print >> config, '%s ! ncall' % ncall
print >> config, '%s ! itmx' % options.itmx
print >> config, '-1.d0 ! acc'
print >> config, '%s ! iseed' % seed
print >> config, '%s ! idist' % options.idist
print >> config, '%s ! itdist' % options.itdist
print >> config, '%s ! iadist' % options.iadist
print >> config, '%s ! ilhe' % options.ilhe


# Filename
filename = '%s_%s%s%s_%s_%s%s' % (name,model,sector,interference,emc_col,ncall,options.tag)
try:
      with open('Config/%s.com' % filename,'w') as cfile1:
            cfile1.write(config.getvalue())
except IOError:
      print "Not in right directory?"
      sys.exit()

# Logfile      
write = '' if options.output else '> Output/%s.log' % (filename)

# Command
command = './%s < Config/%s.com %s &' % (executable,filename,write)
print 'executing: ',command
subprocess.call(command, shell=True)
   


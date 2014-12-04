#!/usr/bin/env python 
# ------------------------------------------------------------------------------------------
# Command line executor for Zp-tt code. Use run.py -h for help
# generates a .com file to feed into fortran Zp-tt executable and runs it

# Arguments:
# 	0 executable e.g. ttbar_BSM. see makefile
#   1 mdl: Name of model file in Models directory
# 	2 ecm: centre of mass energy (TeV)
#   3 ifs: final state
# 	4 ncall: number of vegas calls (monte carlo points ~ 500k should be ok for good statistics)
# ------------------------------------------------------------------------------------------
import os,StringIO,re,optparse,subprocess,time,sys,random

# Options
parser = optparse.OptionParser("usage: %prog [options] executable model ecm_col ifs ncall")
# External
parser.add_option("-o", "--output", action="store_true", default=False, help="output to terminal")
parser.add_option("-t", "--tag", default="",type="string",help="add specific tag to logfile name")

# Internal
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
parser.add_option("-d", "--dists",  action="store_const", const=0, default=1, help="Remove distributions to log file")
parser.add_option("-m", "--itmx",   type="int",default=5, help="Maximum number of VEGAS iterations.")
parser.add_option("-l", "--ilhe",   default=0,action="store_const",const=1, help="Output in lhe format.")
# parser.add_option("-t", "--isycost", default=0,action="store_const",const=1, help="Manually sum over costheta")
# parser.add_option("-i", "--it", action="store_const", default=5, help="Iterations")
# parser.add_option("-p", "--path", default="",type="string",help="specify where to save log files")
# parser.add_option("-x", "--ptcut", default=0,type="float", help="pT cut")
# parser.add_option("-Q", "--qfact", default=0,type="float", help="factorisation scale for PDFs, default = Z mass")
# parser.add_option("-c", "--cosx", action="store_const", const=1,default=0, help="symmetric cos theta and x sums")
# parser.add_option("-r", "--range", default=None, help="specify tuple of min and max energy range over which to integrate as 'min:max'")
########################################################################
(options, args) = parser.parse_args()

executable="ttbar_BSM"
model=args[0]
emc_col=args[1]
ifs=args[2]
ncall=args[3]
name = "2to2" if (ifs=="0") else "2to6"
config1=StringIO.StringIO()
seed = 12345 if options.iseed else random.randint(0,100000)
sector = "_EW" if (options.iqcd==0) else ""
if   options.iint==0:
	interference="_int0"
elif options.iint==1:
	interference="_int1"
elif options.iint==2:
	interference=""
elif options.iint==3:
	interference="_int3"

# invmass=(0,args[2]) if options.range==None else re.findall(r"\d+",options.range)
########################################################################
print >> config1, '%s ! iQCD' % options.iqcd
print >> config1, '%s ! iEW' % options.iew
print >> config1, '%s ! iBSM' % options.ibsm
print >> config1, '%s ! iint' % options.iint
print >> config1, '%s ! ifs' % ifs
print >> config1, '%s ! iNWA' % options.iNWA
print >> config1, '%s ! model' % model
print >> config1, '%s.d3 ! ecm_col' % emc_col
print >> config1, '%s ! icol' % options.icol
print >> config1, '%s ! ISTRUCTURE' % options.pdf
print >> config1, '%s.d0 ! ytcut' % options.ytcut
print >> config1, '%sd0 ! yttcut' % options.yttcut
# print >> config1, '%sd0 ! pTcut' % options.ptcut
# print >> config1, '%s,111,333 ! iseed,jdummy,kdummy' % seed
print >> config1, '%s ! ncall' % ncall
print >> config1, '%s ! itmx' % options.itmx
print >> config1, '-1.d0 ! acc'
print >> config1, '%s ! iseed' % seed
print >> config1, '%s ! idist' % options.dists
print >> config1, '%s ! ilhe' % options.ilhe
# print >> config1, '%s ! isycost' % options.isycost
# print >> config1, '100000.d0 ! alumpb'
# print >> config1, '%sd0 ! min invariant mass' % invmass[0]
# print >> config1, '%sd0 ! max invariant mass' % invmass[1]
# print >> config1, '%s ! dist switch' % options.dists
# print >> config1, '%s ! cosx switch' % options.cosx
# print >> config1, '%s ! select CONT' % options.cont
# print >> config1, '%s ! QCD switch' % options.qcd
# print >> config1, '%s ! EW switch' % options.ew
# print >> config1, '%s ! BSM switch' % options.bsm
# print >> config1, '%s ! final state switch' % options.final
# print >> config1, '%s ! name' % name
# print >> config1, '%s ! model' % args[4]
# print >> config1, '%sd0 ! factorisation scale' % options.qfact
########################################################################
# cosstr=""
# if options.cosx==1:cosstr="_cos_x"
# iecm=int(float(args[2])/1000)
filename = '%s_%s%s%s_%s_%s%s' % (name,model,sector,interference,emc_col,ncall,options.tag)
try:
      with open('Config/%s.com' % filename,'w') as cfile1:
            cfile1.write(config1.getvalue())
except IOError:
      print "Not in right directory?"
      sys.exit()
write = '' if options.output else '> Output/%s.log' % (filename)
command1 = './%s < Config/%s.com %s &' % (executable,filename,write)
print 'executing: ',command1
subprocess.call(command1, shell=True)
   


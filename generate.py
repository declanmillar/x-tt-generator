#!/usr/bin/env python
# -*- coding: utf-8 -*-

# python run script for the generator
# generates a config file, runs generator using it 
# locally or submit to the lxplus/iridis/qmul batch system
# declan.millar@cern.ch

import os
import StringIO
import optparse
import subprocess
import sys
import random
import glob
import socket
import time

parser = optparse.OptionParser()
parser.add_option("-D", "--overwrite",          default = False, action = "store_true",  help = "overwrite existing grid if present.")
parser.add_option("-v", "--verbose",            default = False, action = "store_true",  help = "run in verbose mode.")
parser.add_option("-B", "--batch",              default = True,  action = "store_false", help = "run in batch mode")
parser.add_option("-T", "--ntuple",             default = True,  action = "store_false", help = "write events to ROOT n-tuple")
parser.add_option("-H", "--lhef",               default = False, action = "store_true",  help = "write events to lhe file")
parser.add_option("-g", "--include_gg",         default = False, action = "store_true",  help = "include gluon-gluon initiated interactions")
parser.add_option("-q", "--include_qq",         default = False, action = "store_true",  help = "include quark-quark initiated interactions")
parser.add_option("-d", "--include_dd",         default = False, action = "store_true",  help = "include down-down initiated interactions")
parser.add_option("-u", "--include_uu",         default = False, action = "store_true",  help = "include up-up initiated interactions")
parser.add_option("-A", "--include_a",          default = True,  action = "store_false", help = "include photon mediated interaction")
parser.add_option("-Z", "--include_z",          default = True,  action = "store_false", help = "include Z boson mediated interaction")
parser.add_option("-X", "--include_x",          default = True,  action = "store_false", help = "include Z' boson mediated interactions")
parser.add_option("-s", "--include_signal",     default = True,  action = "store_false", help = "include tt signal")
parser.add_option("-b", "--include_background", default = False, action = "store_true",  help = "include tt background")
parser.add_option("-M", "--multichannel",       default = False, action = "store_true",  help = "use multichannel integration")
parser.add_option("-x", "--symmetrise",         default = False, action = "store_true",  help = "symmetrise phase space around x1<->x2 in integral")
parser.add_option("-R", "--use_rambo",          default = False, action = "store_true",  help = "use RAMBO for phase space")
parser.add_option("-F", "--flatten_integrand",  default = True,  action = "store_false", help = "flatten resonances")
parser.add_option("-W", "--use_nwa",            default = False, action = "store_true",  help = "use Narrow Width Approximation")
parser.add_option("-w", "--unweighted",         default = True,  action = "store_false", help = "unweighted events")
parser.add_option("-f", "--final_state",        default = 1,         type = int,         help = "set final state")
parser.add_option("-i", "--initial_state",      default = 0,         type = int,         help = "initial state: 0 = pp, 1 = pp~")
parser.add_option("-N", "--iterations",         default = 10,        type = int,         help = "number of VAMP iterations")
parser.add_option("-n", "--ncall",              default = 1000000,   type = int,         help = "number of VAMP calls")
parser.add_option("-e", "--nevents",            default = 100000,    type = int,         help = "number of events")
parser.add_option("-P", "--pdf",                default = 11,        type = int,         help = "structure_functions")
parser.add_option("-I", "--interference",       default = 1,         type = int,         help = "specify interference")
parser.add_option("-E", "--energy",             default = 13,        type = int,         help = "collider energy")
parser.add_option("-L", "--energy_low",         default = 0,         type = int,         help = "Ecm lower limit")
parser.add_option("-U", "--energy_up",          default = 0,         type = int,         help = "Ecm upper limit")
parser.add_option("-t", "--tag",                default = "",                            help = "add a name tag to output files")
parser.add_option("-k", "--walltime",           default = "60:00:00",                    help = "walltime 'hh:mm:ss'")
parser.add_option("-Q", "--queue",              default = "1nw",                         help = "lxbatch queue'")
parser.add_option("-m", "--model",              default = "SM",                          help = "set model")
parser.add_option("-c", "--cut",                default = False, action = "store_true",  help = "apply fiducial cuts")
(option, args) = parser.parse_args()

if os.path.isfile("Models/%s.mdl" % option.model) is False: sys.exit("error: %s is not a valid model.\n Available model files: %s" % (option.model, glob.glob("Models/*.mdl")))
if option.energy < 0: sys.exit("error: collider energy must be positive definite\n" % usage)
if option.ncall < 2: sys.exit("error: must have at least 2 VEGAS points\n%s" % usage)
if option.energy_low < 0 or option.energy_up < 0: sys.exit("error: energy bounds must be positive definite")
if option.energy_low > option.energy or option.energy_up > option.energy: sys.exit("error: energy bounds cannot exceed collider energy")
if (option.energy_low > 0 and option.energy_up > 0 and option.energy_up <= option.energy_low): sys.exit("error: energy upper bound must be greater than lower bound")
if option.interference < 0 or option.interference > 4: sys.exit("error: interference must be 0 - 4")
if option.pdf < 1 or option.pdf > 11: sys.exit("error: pdf id must be 1 - 11")
if option.include_background == False and option.include_signal == False: sys.exit("error: signal and background both off")
if option.final_state < -1 or option.final_state > 3: sys.exit("error: invalid final state id" % option.final_state)

initial_states = 0
if option.include_gg: initial_states += 1
if option.include_qq: initial_states += 1
if option.include_dd: initial_states += 1
if option.include_uu: initial_states += 1

if option.lhef and initial_states > 1: sys.exit("error: currently when outputting to LHEF, only one initial state can be active.")

hostname = socket.gethostname()
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

if option.final_state == -1 and option.energy_low == 0:
    option.energy_low = 0.5

if option.use_rambo or option.include_background:
    option.flatten_integrand = False

if option.final_state < 1:
    option.use_nwa = False

executable = "generator"
options = ""

if option.initial_state == 1:
    options += ".ppbar"

pdf = ""
if option.pdf ==  1: pdf = "CTEQ6M"
if option.pdf ==  2: pdf = "CTEQ6D"
if option.pdf ==  3: pdf = "CTEQ6L"
if option.pdf ==  4: pdf = "CTEQ6L1"
if option.pdf ==  5: pdf = "MRS9901"
if option.pdf ==  6: pdf = "MRS9902"
if option.pdf ==  7: pdf = "MRS9903"
if option.pdf ==  8: pdf = "MRS9904"
if option.pdf ==  9: pdf = "MRS9905"
if option.pdf == 10: pdf = "CT14LN"
if option.pdf == 11: pdf = "CT14LL"

if option.interference != 1: options += ".int%i" % option.interference
if option.use_nwa: options += ".nwa"
if option.multichannel: options += ".multi"
if option.cut: options += ".cut"
if option.symmetrise: options += ".symmetrised"
if option.use_rambo: options += ".rambo"
if option.include_background == False and option.flatten_integrand == False: options += ".unmapped"
if option.energy_low != 0 or option.energy_up != 0: options += ".%s-%s" % (option.energy_low, option.energy_up)
if len(option.tag) > 0: options += "." + option.tag
if option.final_state < 2: option.include_background = False
if option.include_background: flatten_integrand = False

initial_partons = ""
if option.include_gg:
    initial_partons += "gg"
if option.include_qq:
    initial_partons += "qq"
if option.include_dd:
    initial_partons += "dd"
if option.include_uu:
    initial_partons += "uu"
initial_partons += "-"

intermediates = ""
if option.final_state < 2:
    if option.include_a: intermediates += "A"
    if option.include_z: intermediates += "Z"
    if option.include_x: intermediates += "X"
else:
    if option.include_background == False and option.include_signal == True:
        intermediates += "tt"
    if option.include_background == True and option.include_signal == False:
        intermediates += ".bkg-only"

if len(intermediates) > 0:
    intermediates = intermediates + "-"

final_state = ""
if   option.final_state == -1: final_state = "ll"
elif option.final_state == 0:  final_state = "tt"
elif option.final_state == 1:  final_state = "tt-bbllvv"
elif option.final_state == 2:  final_state = "tt-blvbqq"
elif option.final_state == 11: final_state = "bbtatavtvt"
elif option.final_state == 12: final_state = "bbemuvevm"

process = initial_partons + intermediates + final_state

filename = '%s.%s.%sTeV.%s%s' % (process, option.model, str(option.energy), pdf, options)

home_directory = "."
data_directory = "."
if "Sunder" in hostname:
    home_directory = "/Users/declan/Code/"
    data_directory = "/Users/declan/Data/"
elif "lxplus" in hostname:
    home_directory = "/afs/cern.ch/user/d/demillar/"
    data_directory = "/afs/cern.ch/work/d/demillar/"
elif "cyan" in hostname:
    home_directory = "/home/dam1g09/"
    data_directory = "/scratch/dam1g09/"
elif "heppc" in hostname:
    home_directory = "/users/millar/"
    data_directory = "/data/millar/"
else:
    exit("error: unknown host")
run_directory = home_directory + "zprime-top-generator"
data_directory = data_directory + "zprime"

if os.path.isdir(data_directory) is False:
    sys.exit("error: specified run directory '%s' does not exist" % run_directory)
    sys.exit("error: specified data directory '%s' does not exist" % data_directory)

if option.multichannel:
    grid = "grids"
else:
    grid = "grid"

if option.unweighted:
    weight = "unweighted"
else:
    weight = "weighted"

config_name = '%s/%s.cfg' % (data_directory, filename)
logfile = "%s/%s.log" % (data_directory, filename)
handler_name = "%s.sh" % filename
ntuple_file = "%s/%s.%s.root" % (data_directory, filename, weight)
lhe_file = "%s/%s.%s.lhe" % (data_directory, filename, weight)
grid_file = "%s/%s.%s" % (data_directory, filename, grid)

if os.path.isfile(grid_file): 
    new_grid = False
else:
    new_grid = True
if option.overwrite:
    new_grid = True

config = StringIO.StringIO()
print >> config, '%r    ! ntuple'             % option.ntuple
print >> config, '%r    ! lhef'               % option.lhef
print >> config, '%r    ! new_grid'           % new_grid
print >> config, '%s'                         % ntuple_file
print >> config, '%s'                         % lhe_file
print >> config, '%s'                         % logfile
print >> config, '%s'                         % grid_file
print >> config, '%i    ! initial state'      % option.initial_state
print >> config, '%i    ! final state'        % option.final_state
print >> config, '%s    ! model'              % option.model
print >> config, '%i    ! pdf'                % option.pdf
print >> config, '%r    ! include signal'     % option.include_signal
print >> config, '%r    ! include background' % option.include_background
print >> config, '%r    ! include gg'         % option.include_gg
print >> config, '%r    ! include qq'         % option.include_qq
print >> config, '%r    ! include uu'         % option.include_uu
print >> config, '%r    ! include dd'         % option.include_dd
print >> config, '%r    ! include a'          % option.include_a
print >> config, '%r    ! include z'          % option.include_z
print >> config, '%r    ! include x'          % option.include_x
print >> config, '%i    ! interference'       % option.interference
print >> config, '%r    ! use nwa'            % option.use_nwa
print >> config, '%i.d3 ! energy'             % option.energy
print >> config, '%i    ! iterations'         % option.iterations
print >> config, '%i    ! ncall'              % option.ncall
print >> config, '%i    ! nevents'            % option.nevents
print >> config, '%r    ! unweighted'         % option.unweighted
print >> config, '%r    ! use rambo'          % option.use_rambo
print >> config, '%r    ! map phase space'    % option.flatten_integrand
print >> config, '%r    ! multichannel'       % option.multichannel
print >> config, '%r    ! symmetrise'         % option.symmetrise
print >> config, '%r    ! verbose mode'       % option.verbose
print >> config, '%i.d3 ! energy low'         % option.energy_low
print >> config, '%i.d3 ! energy up'          % option.energy_up
print >> config, '%r    ! batch mode'         % option.batch
print >> config, '%r    ! fiducial cuts'      % option.cut

try:
    with open('%s' % config_name,'w') as config_file:
        config_file.write(config.getvalue())
        print " config:  %s" % config_name
except IOERROR:
    sys.exit("error: Cannot write to %s" % config_name)

if option.batch:
    handler = StringIO.StringIO()
    if "lxplus" in hostname:
        print >> handler, "#!/bin/bash"
        print >> handler, "source /afs/cern.ch/user/d/demillar/.bash_profile"
        print >> handler, "cd %s" % run_directory
        print >> handler, "%s/bin/%s < %s" % (run_directory, executable, config_name)
    elif "cyan" in hostname:
        print "walltime = %s" % option.walltime
        print >> handler, "#!/bin/bash"
        print >> handler, "source /home/dam1g09/.bash_profile"
        print >> handler, "cd %s" % run_directory
        print >> handler, '%s/bin/%s < %s > %s/%s.log' % (run_directory, executable, config_name, data_directory, filename)
    elif "heppc" in hostname:
        print "h_rt = %s" % option.walltime
        print >> handler, "#!/bin/bash"
        print >> handler, "source /users/millar/.bash_profile"
        print >> handler, "cd %s" % run_directory
        print >> handler, '%s/bin/%s < %s' % (run_directory, executable, config_name)
    else:
        sys.exit("error: hostname not recognised")

    try:
        with open('%s' % handler_name, 'w') as handler_file:
            handler_file.write(handler.getvalue())
        print "Handler file written to %s.sh." % filename
    except IOERROR:
        sys.exit("error: Cannot write handler file.")

    subprocess.call("chmod a+x %s.sh" % filename, shell = True)
    print "submitting batch job ..."
    if "lxplus" in hostname: subprocess.call('bsub -q %s -o %s %s/%s.sh' % (option.queue, logfile, run_directory, filename), shell = True)
    elif "cyan03" in hostname: subprocess.call('qsub -l walltime=%s %s/%s.sh' % (option.walltime, run_directory, filename), shell = True)
    elif "heppc" in hostname: subprocess.call('qsub -l h_rt=%s %s/%s.sh' % (option.walltime, run_directory, filename), shell = True)
    else: 
        print "error: hostname not recognised"
else:
    subprocess.call("./bin/%s < %s | tee %s" % (executable, config_name, logfile), shell = True)

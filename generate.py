#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Python run script for the event generator Fortran executable.
Generates a config file and runs the executable using it.
Runs locally or submits a job to the lxplus/iridis/qmul batch system.
"""

__author__ = "Declan Millar"
__version__ = "dev"
__license__ = "MIT"

import os
import argparse
import subprocess
import fnmatch
import enum


def main():
    """
    Generates a config file and runs the executable using it.
    Runs locally or submits a job to the lxplus/iridis/qmul batch system.
    """

    args = parse_args()

    # Check executable exists
    if not os.path.isfile(f"{args.filename}"):
        raise FileNotFoundError(f"Executable '{args.filename}' does not exist")

    if not os.path.isdir(args.data_dir):
        print(f"Warning: '{args.data_dir}' does not exist; it will be created.")
        os.makedirs(args.data_dir)

    # Check model_name file exists
    model_path = f"./models/{args.model_name}.mdl"
    if not os.path.isfile(model_path):
        raise FileNotFoundError(f"{model_path} does not exist.")

    # Drell-Yan has no gluon-mediated channels
    if args.final_state == -1 and (args.initial_state == 1
                                   or args.initial_state == 2):
        raise ValueError("Drell-Yan has no gluon-mediated channels.")

    # No Z' boson in the Standard Model
    if args.model_name == "SM":
        args.include_x = False

    # No electroweak mediators for gluon-interactions
    if args.initial_state == 1 or args.initial_state == 2:
        args.a, args.z, args.x = False, False, False

    # Avoid divergence in Drell-Yan
    if args.final_state == -1 and args.ecm_low == 0:
        args.ecm_low = 0.5

    # No change of variables to flatten Breit-Wigner peak
    if args.rambo or args.background:
        args.flatten_integrand = False

    # Do not use the Narrow Width Approximation for stable-tops Drell-Yan
    if args.final_state < 1:
        args.nwa = False

    if args.final_state < 2:
        args.background = False

    intermediates = ""
    if args.a:
        intermediates += "A"
    if args.z:
        intermediates += "Z"
    if args.x:
        intermediates += "X"
    if len(intermediates) > 0:
        intermediates += "-"

    initial_states = ""
    if args.gg:
        initial_states += "gg"
    if args.qq:
        initial_states += "qq"
    if args.uu:
        initial_states += "uu"
    if args.dd:
        initial_states += "dd"
    if len(initial_states) > 0:
        initial_states += "-"

    final_state_names = {
        -1: "ll",
        0: "tt",
        1: "tt-bbllvv",
        11: "tt-bbeevv",
        22: "tt-bbmmvv",
        12: "tt-bbemvv",
        21: "tt-bbmevv",
        33: "bbtatavv",
    }

    final_state_name = final_state_names[args.final_state]

    if args.initial_state == 0 and args.final_state > 0:
        process = "2-6-phase-space"
    elif args.initial_state == 0 and args.final_state < 1:
        process = "2-2-phase-space"
    else:
        process = initial_states + intermediates + final_state_name

    gridproc = initial_states + intermediates + "tt-bbllvv"

    pdf_name = PDF(args.pdf).name

    # Include non-default major options_label in output file name
    options_label = ""
    if args.ppbar == 1:
        options_label += "_ppbar"
    if args.interference != 1:
        options_label += f"_int{args.interference}"
    if args.nwa:
        options_label += "_nwa"
    if args.cut:
        options_label += "_cut"
    if args.rambo:
        options_label += "_rambo"
    if args.background is False and args.flatten_integrand is False:
        options_label += "_unmapped"
    if args.ecm_low is not None or args.ecm_up is not None:
        options_label += f"_{args.ecm_low}-{args.ecm_up}"
    if args.tag is not None:
        options_label += "_" + args.tag

    events_name = f"{process}_{args.model_name}_{str(args.sqrts)}TeV_{pdf_name}{options_label}"
    grid_name = f"{gridproc}_{args.model_name}_{str(args.sqrts)}TeV_{pdf_name}{options_label}"

    # Get new file index
    datafiles = [
        f for f in os.listdir(args.data_dir)
        if os.path.isfile(os.path.join(args.data_dir, f))
    ]
    filtered = fnmatch.filter(datafiles, events_name + "_??.lhef")
    new_index = int(
        args.index) if args.index is not None else len(filtered) + 1
    events_name = events_name + f"_{new_index:03d}"

    events_path = args.data_dir + events_name
    grid_path = args.data_dir + grid_name

    wgt = "" if args.unweighted else "_wgt"

    grid_file = f"{grid_path}.grid"
    xsec_file = f"{grid_path}.txt"

    new_grid = not os.path.isfile(grid_file) or args.overwrite

    events_path = events_path + wgt
    lhe_file = events_path + ".lhef"

    if new_grid:
        config_name = grid_path + ".cfg"
        logfile = grid_path + ".log"
        handler_name = grid_path + ".sh"
    else:
        config_name = events_path + ".cfg"
        logfile = events_path + ".log"
        handler_name = events_path + ".sh"

    config = (f"{new_grid} ! new_grid\n"
              f"{grid_file}\n"
              f"{xsec_file}\n"
              f"{logfile}\n"
              f"{lhe_file}\n"
              f"{args.ppbar} ! initial_state\n"
              f"{args.final_state} ! final_state\n"
              f"{args.model_name} ! model_name\n"
              f"{args.pdf} ! pdf\n"
              f"{args.signal} ! include_signal\n"
              f"{args.background} ! include_background\n"
              f"{args.gg} ! include_gg\n"
              f"{args.qq} ! include_qq\n"
              f"{args.uu} ! include_uu\n"
              f"{args.dd} ! include_dd\n"
              f"{args.a} ! include_a\n"
              f"{args.z} ! include_z\n"
              f"{args.x} ! include_x\n"
              f"{args.interference} ! interference\n"
              f"{args.nwa} ! nwa\n"
              f"{args.sqrts}.d3 ! sqrts\n"
              f"{args.itmx} ! itmx\n"
              f"{args.ncall} ! ncall\n"
              f"{args.nevents} ! nevents\n"
              f"{args.unweighted} ! unweighted\n"
              f"{args.rambo} ! use_rambo\n"
              f"{args.flatten_integrand} ! flatten_integrand\n"
              f"{args.verbose} ! verbose\n"
              f"{0 if args.ecm_low is None else args.ecm_low}.d3 ! ecm_low\n"
              f"{0 if args.ecm_up is None else args.ecm_up}.d3 ! ecm_up\n"
              f"{args.job} ! batch\n"
              f"{args.cut} ! cut\n")

    with open(config_name, "w") as config_file:
        config_file.write(config)

    print(f"Config: {config_name}")
    print(config)

    if args.job:
        print(f"walltime: {args.walltime}")
        handler = ("#!/bin/bash\n"
                   "source /home/dam1g09/.bash_profile\n"
                   "module load gcc/6.1.0\n"
                   "module load openmpi/2.0.2/gcc\n"
                   "module load intel/2017\n"
                   "module load intel/mpi/2017\n"
                #    f"cd {run_directory}\n"
                #    f"{run_directory}{executable} < {config_name} > {logfile}\n"
                   f"{args.filename} < {config_name} > {logfile}\n"
                   f"gzip -v9 {lhe_file} >> {logfile}\n")

        with open(handler_name, "w") as handler_file:
            handler_file.write(handler)

        print(f"handler: {handler_name}")

        subprocess.call(f"chmod a+x {handler_name}", shell=True)
        subprocess.call(
            f"qsub -l walltime={args.walltime} {handler_name}",
            shell=True)
    else:
        command = (
            f"{args.filename} < {config_name} | tee {logfile}")
            # f" && gzip -v9 {lhe_file} >> {logfile}")
        subprocess.call(command, shell=True)


def parse_args():
    """
    Process command line arguments.

    Interference
    0: (gamma) + (Z) + (Z')
    1: (gamma + Z + Z')
    2: (gamma + Z) + (Z')
    3: (gamma + Z + Z') - (gamma) - (Z)
    4: (gamma + Z + Z') - (gamma) - (Z) - (Z')
    """

    parser = argparse.ArgumentParser(description="Generate X-tt events.")

    # File and directory options
    parser.add_argument("-f", "--filename", help="Executable file name.",
                        default="./bin/generator")
    parser.add_argument("-d", "--data_dir", help="Executable file name.",
                        default="./data/")

    # Job options
    parser.add_argument("-w", "--walltime",
                        help="Set Iridis walltime in 'hh:mm:ss' format",
                        default="60:00:00")
    parser.add_argument("-j", "--job", help="Submit as a batch job.",
                        action="store_true")
    parser.add_argument("-o", "--overwrite",
                        help="Overwrite existing grid file.",
                        action="store_true")
    parser.add_argument("-v", "--verbose", help="Print extra run information.",
                        action="store_true")
    parser.add_argument("--tag", help="Add a name tag to output files.",
                        default=None)
    parser.add_argument("--index", help="Overwrite filename index tag.",
                        default=None)
    # parser.add_argument("-q",
    #                     "--queue",
    #                     help="lxbatch queue",
    #                     default="1nw") # No access to lxplus

    # Collider options
    parser.add_argument("-s", "--sqrts", help="Set collider energy (s hat).",
                        type=int, default=13)
    parser.add_argument("--ecm_low",
                        help="Set a lower limit on the collider energy.",
                        type=int, default=None)
    parser.add_argument("--ecm_up",
                        help="Set an upper limit on the collider energy",
                        type=int, default=None)
    parser.add_argument("--ppbar",
                        help="Set initial collider state. 0: pp, 1: pp~",
                        default=0)

    # Physics model
    parser.add_argument("-m", "--model_name",
                        help="Set the physics model name.", default="SM")

    # Final/initial state
    parser.add_argument("-F", "--final_state", help="Set final state ID.",
                        type=int, default=1)
    parser.add_argument("-I", "--initial_state", help="Set ititial state ID.",
                        type=int, default=1)
    parser.add_argument("-G", "--gg",
                        help="Include gluon-gluon initial state.",
                        action="store_true")
    parser.add_argument("-Q", "--qq",
                        help="Include quark-quark initial state.",
                        action="store_true")
    parser.add_argument("-U", "--uu", help="Include up-up initial state.",
                        action="store_true")
    parser.add_argument("-D", "--dd", help="Include down-down initial state.",
                        action="store_true")

    # Mediators
    parser.add_argument("-A", "--a",
                        help="Include photon mediated interaction.",
                        action="store_false")
    parser.add_argument("-Z", "--z",
                        help="Include Z boson mediated interaction.",
                        action="store_false")
    parser.add_argument("-X", "--x",
                        help="Include Z' boson mediated interactions.",
                        action="store_false")

    # QFT interference
    parser.add_argument("--interference", help="Set interference", type=int,
                        default=1)

    # 2to6 signal and background
    parser.add_argument("--signal", help="Include ttbar 2to6 signal.",
                        action="store_false")
    parser.add_argument("--background", help="Include ttbar 2to6 background.",
                        action="store_true")

    # Parton Distribution Functions
    parser.add_argument("--pdf", help="Parton Density Function ID.", type=int,
                        default=11)

    # Integration options
    parser.add_argument("-N", "--itmx", help="Number of VAMP iterations.",
                        type=int, default=16)
    parser.add_argument("-n", "--ncall", help="Number of VAMP calls.",
                        type=int, default=16000000)
    parser.add_argument("-e", "--nevents", help="Number of unweighted events.",
                        default=10000)
    parser.add_argument("--flatten_integrand", help="Flatten integrand.",
                        action="store_false")
    parser.add_argument("--unweighted", help="Generate unweighted events.",
                        action="store_false")
    parser.add_argument("--rambo", help="Use RAMBO for phase space.",
                        action="store_true")
    parser.add_argument("--nwa", help="Use the Narrow Width Approximation.",
                        action="store_true")
    parser.add_argument("--cut", help="Apply detector cuts.",
                        action="store_true")

    return parser.parse_args()


class PDF(enum.Enum):
    """
    Parton distribution function enum.
    """

    CTEQ6M = 1
    CTEQ6D = 2
    CTEQ6L = 3
    CTEQ6L1 = 4
    MRS9901 = 5
    MRS9902 = 6
    MRS9903 = 7
    MRS9904 = 8
    MRS9905 = 9
    CT14LN = 10
    CT14LL = 11


if __name__ == "__main__":
    main()

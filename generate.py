#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Python script to run the event generator executable.
Generates a config file and runs the executable using it.
Runs locally or submits a job to the lxplus/iridis/qmul batch system.
"""

import os
import argparse
import subprocess
import fnmatch
import enum


def check_config(config) -> None:
    """
    Check that config is valid
    """

    # Check executable exists
    if not os.path.isfile(config["filename"]):
        raise FileNotFoundError(f"Executable {config['filename']} does not exist")

    if not os.path.isdir(config["data_dir"]):
        print(f"Warning: {config['data_dir']} does not exist; it will be created.")
        os.makedirs(config["data_dir"])

    # Check model_name file exists
    model_path = f"./models/{config['model_name']}.mdl"
    if not os.path.isfile(model_path):
        raise FileNotFoundError(f"{model_path} does not exist.")

    # Drell-Yan has no gluon-mediated channels
    if config["final_state"] == -1 and (
        config["initial_state"] == 1 or config["initial_state"] == 2
    ):
        raise ValueError("Drell-Yan has no gluon-mediated channels.")


def attune_config(config) -> None:
    """
    Ensure config is consistent
    """

    # No Z' boson in the Standard Model
    if config["model_name"] == "SM":
        config["include_x"] = False

    # No electroweak mediators for gluon interactions
    if config["initial_state"] == 1 or config["initial_state"] == 2:
        config["A"], config["Z"], config["X"] = False, False, False

    # Avoid divergence in Drell-Yan
    if config["final_state"] == -1 and config["ecm_low"] == 0:
        config["ecm_low"] = 0.5

    # No change of variables to flatten Breit-Wigner peak
    if config["rambo"] or config["background"]:
        config["flatten_integrand"] = False

    # Do not use the Narrow Width Approximation for stable-tops Drell-Yan
    if config["final_state"] < 1:
        config["nwa"] = False

    if config["final_state"] < 2:
        config["background"] = False


def create_output_names(config):
    """
    Create names for all output files.
    """

    initial_states = ["gg", "qq", "uu", "dd"]
    process = ""

    for state in initial_states:
        if config[state]:
            process += state

    process += "-"

    intermediates = ["A", "Z", "X"]

    for state in initial_states:
        if config[state]:
            process += state

    for state in intermediates:
        if config[state]:
            process += state

    process += "-"

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

    if config["final_state"] > 0:
        grid_proc = process + final_state_names[1]

    process = process + final_state_names[config["final_state"]]

    pdf_name = PDF(config["pdf"]).name

    # Include non-default major tag in output file name
    tag = ""
    if config["ppbar"] == 1:
        tag += "_ppbar"
    if config["interference"] != 1:
        tag += f"_int{config['interference']}"
    if config["nwa"]:
        tag += "_nwa"
    if config["cut"]:
        tag += "_cut"
    if config["rambo"]:
        tag += "_rambo"
    if config["background"] is False and config["flatten_integrand"] is False:
        tag += "_unmapped"
    if config["ecm_low"] is not None or config["ecm_up"] is not None:
        tag += f"_{config['ecm_low']}-{config['ecm_up']}"
    if config["tag"] is not None:
        tag += "_" + config["tag"]

    events_name = f"{process}_{config['model_name']}_{str(config['sqrts'])}TeV_{pdf_name}{tag}"
    grid_name = f"{grid_proc}_{config['model_name']}_{str(config['sqrts'])}TeV_{pdf_name}{tag}"

    # Get new file index
    datafiles = [
        f
        for f in os.listdir(config["data_dir"])
        if os.path.isfile(os.path.join(config["data_dir"], f))
    ]
    filtered = fnmatch.filter(datafiles, events_name + "_??.lhef")
    new_index = int(config["index"]) if config["index"] is not None else len(filtered) + 1
    events_name = events_name + f"_{new_index:03d}"

    events_path = config["data_dir"] + events_name
    grid_path = config["data_dir"] + grid_name

    wgt = "" if config["unweighted"] else "_wgt"

    grid_name = f"{grid_path}.grid"
    xsec_file = f"{grid_path}.txt"

    new_grid = not os.path.isfile(grid_name) or config["overwrite"]

    events_path = events_path + wgt
    lhef = events_path + ".lhef"

    if new_grid:
        config_name = grid_path + ".cfg"
        log = grid_path + ".log"
        handler = grid_path + ".sh"
    else:
        config_name = events_path + ".cfg"
        log = events_path + ".log"
        handler = events_path + ".sh"

    return config_name, handler, new_grid, grid_name, xsec_file, log, lhef


def main():
    """
    Generates a config file and runs the executable using it.
    Runs locally or submits a job to the lxplus/iridis/qmul batch system.
    """

    args = parse_args()

    config = vars(args)
    print("initial config:\n", config)

    check_config(config)
    attune_config(config)

    print("checked config:\n", config)

    config_name, handler, new_grid, grid_name, xsec_file, log, lhef = create_output_names(config)

    config_string = (
        f"{new_grid} ! new_grid\n"
        f"{grid_name}\n"
        f"{xsec_file}\n"
        f"{log}\n"
        f"{lhef}\n"
        f"{config['ppbar']} ! initial_state\n"
        f"{config['final_state']} ! final_state\n"
        f"{config['model_name']} ! model_name\n"
        f"{config['pdf']} ! pdf\n"
        f"{config['signal']} ! include_signal\n"
        f"{config['background']} ! include_background\n"
        f"{config['gg']} ! include_gg\n"
        f"{config['qq']} ! include_qq\n"
        f"{config['uu']} ! include_uu\n"
        f"{config['dd']} ! include_dd\n"
        f"{config['A']} ! include_a\n"
        f"{config['Z']} ! include_z\n"
        f"{config['X']} ! include_x\n"
        f"{config['interference']} ! interference\n"
        f"{config['nwa']} ! nwa\n"
        f"{config['sqrts']}.d3 ! sqrts\n"
        f"{config['itmx']} ! itmx\n"
        f"{config['ncall']} ! ncall\n"
        f"{config['nevents']} ! nevents\n"
        f"{config['unweighted']} ! unweighted\n"
        f"{config['rambo']} ! use_rambo\n"
        f"{config['flatten_integrand']} ! flatten_integrand\n"
        f"{config['verbose']} ! verbose\n"
        f"{0 if config['ecm_low'] is None else config['ecm_low']}.d3 ! ecm_low\n"
        f"{0 if config['ecm_up'] is None else config['ecm_up']}.d3 ! ecm_up\n"
        f"{config['job']} ! batch\n"
        f"{config['cut']} ! cut\n"
    )

    with open(config_name, "w") as config_file:
        config_file.write(config_string)

    print(f"Config: {config_name}")
    print(config)

    if args.job:
        print(f"walltime: {args.walltime}")
        handler = (
            "#!/bin/bash\n"
            "source /home/dam1g09/.bash_profile\n"
            "module load gcc/6.1.0\n"
            "module load openmpi/2.0.2/gcc\n"
            "module load intel/2017\n"
            "module load intel/mpi/2017\n"
            #    f"cd {run_directory}\n"
            #    f"{run_directory}{executable} < {config_name} > {log}\n"
            f"{args.filename} < {config_name} > {log}\n"
            f"gzip -v9 {lhef} >> {log}\n"
        )

        with open(handler, "w") as handler_file:
            handler_file.write(handler)

        print(f"handler: {handler}")

        subprocess.call(f"chmod a+x {handler}", shell=True)
        subprocess.call(f"qsub -l walltime={args.walltime} {handler}", shell=True)
    else:
        command = f"{args.filename} < {config_name} | tee {log}"
        # f" && gzip -v9 {lhef} >> {log}")
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
    parser.add_argument("-f", "--filename", help="Executable file name.", default="./bin/generator")
    parser.add_argument("-d", "--data_dir", help="Executable file name.", default="./data/")

    # Job options
    parser.add_argument(
        "-w", "--walltime", help="Set Iridis walltime in 'hh:mm:ss' format", default="60:00:00"
    )
    parser.add_argument("-j", "--job", help="Submit as a batch job.", action="store_true")
    parser.add_argument(
        "-o", "--overwrite", help="Overwrite existing grid file.", action="store_true"
    )
    parser.add_argument("-v", "--verbose", help="Print extra run information.", action="store_true")
    parser.add_argument("--tag", help="Add a name tag to output files.", default=None)
    parser.add_argument("-i", "--index", help="Overwrite filename index tag.", default=None)
    # parser.add_argument("-q",
    #                     "--queue",
    #                     help="lxbatch queue",
    #                     default="1nw") # No access to lxplus

    # Collider options
    parser.add_argument("-s", "--sqrts", help="Set collider energy (s hat).", type=int, default=13)
    parser.add_argument(
        "--ecm_low", help="Set a lower limit on the collider energy.", type=int, default=None
    )
    parser.add_argument(
        "--ecm_up", help="Set an upper limit on the collider energy", type=int, default=None
    )
    parser.add_argument("--ppbar", help="Set initial collider state. 0: pp, 1: pp~", default=0)

    # Physics model
    parser.add_argument("-m", "--model_name", help="Set the physics model name.", default="SM")

    # Final/initial state
    parser.add_argument("-F", "--final_state", help="Set final state ID.", type=int, default=1)
    parser.add_argument("-I", "--initial_state", help="Set ititial state ID.", type=int, default=1)
    parser.add_argument(
        "-G", "--gg", help="Include gluon-gluon initial state.", action="store_true"
    )
    parser.add_argument(
        "-Q", "--qq", help="Include quark-quark initial state.", action="store_true"
    )
    parser.add_argument("-U", "--uu", help="Include up-up initial state.", action="store_true")
    parser.add_argument("-D", "--dd", help="Include down-down initial state.", action="store_true")

    # Mediators
    parser.add_argument("-A", help="Include photon mediated interaction.", action="store_false")
    parser.add_argument("-Z", help="Include Z boson mediated interaction.", action="store_false")
    parser.add_argument("-X", help="Include Z' boson mediated interactions.", action="store_false")

    # QFT interference
    parser.add_argument("--interference", help="Set interference", type=int, default=1)

    # 2to6 signal and background
    parser.add_argument("--signal", help="Include ttbar 2to6 signal.", action="store_false")
    parser.add_argument("--background", help="Include ttbar 2to6 background.", action="store_true")

    # Parton Distribution Functions
    parser.add_argument("--pdf", help="Parton Density Function ID.", type=int, default=11)

    # Integration options
    parser.add_argument("-N", "--itmx", help="Number of VAMP iterations.", type=int, default=16)
    parser.add_argument("-n", "--ncall", help="Number of VAMP calls.", type=int, default=16000000)
    parser.add_argument("-e", "--nevents", help="Number of unweighted events.", default=10000)
    parser.add_argument("--flatten_integrand", help="Flatten integrand.", action="store_false")
    parser.add_argument("--unweighted", help="Generate unweighted events.", action="store_false")
    parser.add_argument("--rambo", help="Use RAMBO for phase space.", action="store_true")
    parser.add_argument("--nwa", help="Use the Narrow Width Approximation.", action="store_true")
    parser.add_argument("--cut", help="Apply detector cuts.", action="store_true")

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

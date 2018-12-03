# Prophet

Produce ResOnant Particles to High Energy Tops

The generation tool employed for our study is a custom Monte Carlo (MC) program. The matrix element calculations are based on helicity amplitudes using [HELAS](http://inspirehep.net/record/336604?ln=en) subroutines, with Standard Model square matrix elements built up using [MagGraph](http://madgraph.physics.illinois.edu). Beyond the Standard Model amplitudes are then constructed by modifying these as required. [Vegas AMPlified (VAMP)](http://www.sciencedirect.com/science/article/pii/S001046559900209X?via%3Dihub), an enhanced version of the popular [VEGAS](https://en.wikipedia.org/wiki/VEGAS_algorithm) program, is used for the multi-dimensional numerical phase-space integration, and the generation of unweighted events.

A number of different PDF sets are available. The most recent of these are the [CT14](http://hep.pa.msu.edu/cteq/public/index.html) leading order tables, with [CTEQ6](http://hep.pa.msu.edu/cteq/public/cteq6.htmlmrs) and MRS99 available for comparison. Generally we select the CT14LO(LL) table for our simulations, with a factorisation/renormalisation scale at twice the top mass. The bottom and top quarks are assigned masses of 4.18 GeV and 172.5 GeV, respectively, while the lighter quarks are treated in the massless limit.

The program can write the minimal event information (event weight, PDG particle IDs, and 4-vectors) directly to a [ROOT](https://root.cern.ch) n-tuple, in binary format, which minimises storage space and eases a parton-level root analysis. Alternatively, if one wishes to further process the events with a parton shower/hadronisation tool (e.g. [pythia](http://home.thep.lu.se/~torbjorn/Pythia.html)), an output text file in the standard [Les Houches Event Format (LHEF)](https://arxiv.org/abs/hep-ph/0609017) can be produced.

The full source code repository is stored [here](https://gitlab.cern.ch/demillar/zprime-top-generator), but currently accessible to collaborators only.

## Running the program

The program should be built by running `make util` followed by 'make'.
The util path should be added to the library path:

```bash
    export LD_LIBRARY_PATH="<path to repo>/util/build/src:$LD_LIBRARY_PATH"
```

The program is executed via the `generate.py` run file. Do `generate.py -h` for the available options.
If run on `lxplus` or `iridis` the program will create a steering file and submit a batch job.

Example for submitting multiple jobs
```bash
    for i in `seq 400 499`; do ./macros/generate.py -i $i -u -f 11; done
```

## Directory Structure

* `bin/`: compiled and linked `generator` binary executable
* `Diagrams/`: Feynman diagrams for the generation processes
* `lib/`: compiled library files
* `macros/`: scripts for running program and more
* `Models/`: the input files for each BSM model
* `PDFs/`: the tables for the available Parton Distribution functions
* `src/`: the fortran (freeform `.f90`) source files
* `util/`: contains dependencies necessary to directly write events to ROOT files

## Important Files

* `generate.py`: The run file for execution.
* `generator.f90`: The main source file where the Fortran program lives.
* `scattering.f90`: Contains the code for calculating the differential cross section.
* `README.md`: This file!


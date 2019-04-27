# Zprime ttbar phenomenology/Generator

The generation tool employed for our study is a custom Monte Carlo (MC) program.
The matrix element calculations are based on helicity amplitudes using [HELAS](http://inspirehep.net/record/336604?ln=en) subroutines, with Standard Model square matrix elements built up using [MagGraph](http://madgraph.physics.illinois.edu). Beyond the Standard Model amplitudes are then constructed by modifying these as required. [Vegas AMPlified (VAMP)](http://www.sciencedirect.com/science/article/pii/S001046559900209X?via%3Dihub), an enhanced version of the popular [VEGAS](https://en.wikipedia.org/wiki/VEGAS_algorithm) program, is used for the multi-dimensional numerical phase-space integration, and the generation of unweighted events.

A number of different PDF sets are available.
The most recent of these are the [CT14](http://hep.pa.msu.edu/cteq/public/index.html) leading order tables, with [CTEQ6](https://hep.pa.msu.edu/cteq/public/cteq6.html) and [MRS99](https://arxiv.org/abs/hep-ph/9906231) available for comparison. Generally we select the CT14LO(LL) table for our simulations, with a factorisation/renormalisation scale at twice the top mass. The bottom and top quarks are assigned masses of 4.18 GeV and 172.5 GeV, respectively, while the lighter quarks are treated in the massless limit.

The program can write the minimal event information (event weight, PDG particle IDs, and 4-vectors) directly to a [ROOT](https://root.cern.ch) n-tuple, in binary format, which minimises storage space and eases a parton-level root analysis. Alternatively, if one wishes to further process the events with a parton shower/hadronisation tool (e.g. [pythia](http://home.thep.lu.se/~torbjorn/Pythia.html)), an output text file in the standard [Les Houches Event Format (LHEF)](https://arxiv.org/abs/hep-ph/0609017) can be produced.

## Research pipeline

Generator -> 
[Delphes](https://gitlab.com/zprime-ttbar-phenomenology/delphes) -> 
[Analysis](https://gitlab.com/zprime-ttbar-phenomenology/analysis) -> 
[Statistics](https://gitlab.com/zprime-ttbar-phenomenology/statistics)

## Running the program

The program is executed via the `generate.py` run file.
Do `generate.py -h` for the available options.
If run on [LXPLUS](http://information-technology.web.cern.ch/services/lxplus-service) or [Iridis](https://www.southampton.ac.uk/isolutions/staff/iridis.page) this will create a script and submit a batch job.

Example for submitting multiple jobs:

```sh
for i in `seq 400 499`
  do
    ./generate.py -i $i -u -f 11
  done
```

## Directory structure

```
generator/
├── Diagrams/            ─ Feynman diagrams for the generation processes
├── Models/              ─ input model files for each BSM model
├── PDFs/                ─ data for the Parton Distribution Functions
├── lib/                 ─ compiled object and Fortran module files
├── scripts/             ─ shell scripts for checking and modifying output
├── src/                 ─ Fortran source files
│   ├── generator.f90    ─ main source file 
│   └── scattering.f90   ─ code for calculating the differential cross section
└── generate.py          ─ run file for execution
```

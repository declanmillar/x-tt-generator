# Z' ttbar Generator

*NOTE: Currently undergoing refactoring and cleanup. Not guaranteed to be in a working state.*

A custom Monte Carlo (MC) program using matrix element calculations based on helicity amplitudes
that use [HELAS](https://inspirehep.net/record/336604?ln=en) subroutines. The Standard Model (SM)
square matrix elements were constructed using [MadGraph](https://madgraph.physics.illinois.edu).

Beyond the Standard Model (BSM) amplitudes were constructed by modifying SM process files as
required. Vegas AMPlified
([VAMP](https://www.sciencedirect.com/science/article/pii/S001046559900209X?via%3Dihub), an enhanced
version of the popular [VEGAS](https://en.wikipedia.org/wiki/VEGAS_algorithm) program, is used for
the multi-dimensional numerical phase-space integration, and the generation of unweighted events.

## Installation

Compilation requires `gfortran`. On macOS, one can install this using [homebrew](https://brew.sh/)'s
`gcc` formulae.

```sh
brew install gcc
```

One can then compile the executable file using the Makefile.

```sh
make
```

## Usage

The executable can be run directly.

```sh
./bin/generator < <your-config-file>.cfg
```

For convenience one may also execute it using the Python `generate.py` file, which generates a
config file and runs the executable with it.

If run on [LXPLUS](https://information-technology.web.cern.ch/services/lxplus-service) or
[Iridis](https://www.southampton.ac.uk/isolutions/staff/iridis.page) this file will create a script
and submit a batch job.

An example to submit multiple jobs is given below.

```sh
for i in `seq 00 99`
  do
    ./generate.py -i $i -U -F 11
  done
```

Run `generate.py -h` to see all the available options.

## Directory structure

```txt
x-tt-generator/
├── diagrams/            ─ Feynman diagrams for the generation processes
├── models/              ─ Input model files for each BSM model
├── pdfs/                ─ Data for the Parton Distribution Functions
├── shell-scripts/       ─ Shell scripts for checking and modifying output
├── src/                 ─ Fortran source files
│   ├── generator.f90    ─ Main source file
│   └── scattering.f90   ─ Code for calculating the differential cross section
└── generate.py          ─ Runfile for executable
```

## Parton distribution functions (PDFs)

A number of different parton distribution function (PDF) sets are available. The most recent of
these are the [CT14](https://hep.pa.msu.edu/cteq/public/index.html) leading-order tables, with
[CTEQ6](https://hep.pa.msu.edu/cteq/public/cteq6.html) and
[MRS99](https://arxiv.org/abs/hep-ph/9906231) available for comparison.

By default, the CT14LO(LL) table is selected with a factorization and renormalization scale at
double the top mass. The bottom and top quarks are assigned masses of 4.18 GeV and 172.5 GeV,
respectively, while the lighter quarks are treated in the massless limit.

## Output

Output can be written to a plain text file in the standard Les Houches Event Format
([LHEF](https://arxiv.org/abs/hep-ph/0609017). This allows further processing of the events with a
parton shower and hadronization tool, such as
[pythia](http://home.thep.lu.se/~torbjorn/Pythia.html).

Output can also be written using the minimal event information (event weight, PDG particle IDs, and
4-vectors) directly to a [ROOT](https://root.cern.ch) n-tuple, in binary format, which minimizes
storage space and eases a parton-level root analysis. I'm planning to remove this functionality.

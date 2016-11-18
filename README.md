# Z' -> tt Generator

---

## Overview

The generation tool employed for our study is a custom Monte Carlo (MC) program. The matrix element calculations are based on helicity amplitudes using HELAS subroutines, with Standard Model square matrix elements built up using MagGraph [hagiwara2000, Stelzer1994]. Beyond the Standard Model amplitudes are then constructed by modifying these as required. Vegas AMPliﬁed (VAMP), an enhanced version of the popular VEGAS program, is used for the multi-dimensional numerical phase-space integration, and the generation of unweighted events [Lepage1980].

A number of diﬀerent PDF sets are available. The most recent of these are the CT14 leading order tables, with CTEQ6 and MRS99 available for comparison [Pumplin2002]. Generally we select the CT14LO(LL) table for our simulations, with a factorisation/renormalisation scale of Q=μ=2mt. The b and t quarks are assigned masses of 4.18 GeV and 172.5 GeV, respectively, while the lighter quarks are treated in the massless limit.

The program can write the minimal event information (event weight, PDG particle IDs, and 4-vectors) directly to a ROOT n-tuple, in binary format, using RootTuple, which minimises storage space and eases a parton-level root analysis. Alternatively, if one wishes to further process the events with a parton shower/hadronisation tool, an output text ﬁle in the standard Les Houches Event can be produced.

The full source code repository is stored here: https://gitlab.cern.ch/demillar/zprime-top-generator (accessible to collaborators only).
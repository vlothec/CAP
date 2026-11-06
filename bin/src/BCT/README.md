# BCT C++ implementation

This directory contains C++ source for the Bayesian Context Trees (BCT)
algorithm, adapted from the implementation by Ioannis Papageorgiou:

- C++ repository: https://github.com/IoannisPapageorgiou/Bayesian-Context-Trees
- related R package: https://cran.r-project.org/package=BCT
- associated paper: https://arxiv.org/abs/2007.14900
- CAP introduction commit: `57f08aa5d72e44ee33b90b62da6b06fe6cc8ce20`

The upstream algorithms include Context Tree Weighting (CTW), Bayesian Context
Trees (BCT), and k-BCT. CAP adds `main_ctw.cpp` and this directory's Makefile to
build the CTW-specific command-line executable `bin/ctw-calc`. CAP invokes that
executable as a separate process.

## Build

From the CAP repository root, run:

```bash
make -C bin/src/BCT
```

The Makefile writes the executable to `bin/ctw-calc`.

## License and redistribution

The standalone C++ repository does not contain a separate explicit license
notice. The related BCT R package published by the same author declares
`GPL (>= 2)` in its CRAN metadata. On that basis, CAP applies this SPDX license
expression to the adapted source in this directory and the resulting
`bin/ctw-calc` executable:

```text
GPL-2.0-or-later
```

The complete GPL version 2 and version 3 texts are in the CAP repository at:

- `licenses/GPL-2.0.txt`
- `licenses/GPL-3.0.txt`

CAP's remaining original code is licensed under MIT. Distribution of a compiled
`ctw-calc` must include the complete corresponding source used to build it,
including CAP's modifications and build instructions. This tracked
`bin/src/BCT` directory is that source for CAP builds.

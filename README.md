# Disjunctive Sum of Squares (DiSOS)

This repository implements the disjunctive sum of squares (DiSOS) method for certifying nonnegativity of polynomials. Unlike traditional sum of squares methods that use a single algebraic identity, DiSOS method uses multiple algebraic identities with the same degree as the polynomial being verified. The implementation integrates branch-and-bound scheme and includes other applications like certifying copositivity of matrices and computing the clique number of graphs.

The `julia/` folder is a git submodule of the Julia package [`DisjunctiveSOS.jl`](https://github.com/stellatogrp/DisjunctiveSOS.jl), which bundles two complementary workflows from the paper:

- **Rational disjunctive sos certificates** — a symbolic, exact-arithmetic pipeline that produces rational `(s_1, s_2, s_3, s_4)` certificates for the polynomials that admit one (`julia/src/`, `julia/scripts/`).
- **Spatial branch-and-bound** — a Julia reimplementation of the BnB framework with the symmetry-aware deduplication and global projected-gradient upper bound enhancements, covering polynomial minimisation on the sphere, copositive optimisation, and the Motzkin–Straus clique problem; reproduces the paper's numerical tables and convergence figures (`julia/bnb/`).

To pull the submodule contents, clone recursively:

```bash
git clone --recurse-submodules https://github.com/YH7422/DisjunctiveSOS.git
```

or, if you already cloned without `--recurse-submodules`, run `git submodule update --init` inside the checkout.

## Repository Structure

```
.
├── data_structure/
│   ├── MinHeap_BnB.m    # Implementation of Min Heap for Branch and Bound
│   └── Node.m           # Node class implementation
│
├── certifying_nonnegativity/
│   ├── DiSOS_BnB_SD.m          # Main DiSOS Branch-and-Bound implementation
│   ├── test_nonsos.m           # Test script for classic non-SOS polynomials
│   ├── test_nonsos_family.m    # Test script for a familiy of non-SOS polynomials
│   └── plot_nonnegativity.m    # Script for producing plots from the paper
│
├── certifying_copositivity/
│   ├── DiSOS_copositive_BnB.m  # DiSOS implementation for matrix copositivity
│   ├── test_standard_qp.m      # Test script for standard quadratic programming
│   └── plot_copositivity.m     # Script for producing plots from the paper
│
└── clique_number/
    ├── DiSOS_clique_BnB.m      # DiSOS implementation for clique number computation
    ├── test_random_graph.m     # Test script for Erdos-Renyi random graphs
    └── Random_graph/           # Subfolder containing graph data
```

## Features

- Implementation of Disjunctive Sum of Squares (DiSOS) method
- Branch and Bound algorithm integration
- Applications in:
  - Polynomial minimization
  - Copositive programming
  - Graph clique number computation

## Requirements

The following versions were used in our experiments. Other versions might work but haven't been tested.

- MATLAB: 9.12 (R2022a)
- CVX: 2.2
- YALMIP: [R20230622](https://github.com/yalmip/YALMIP/releases/tag/R20230622)
- Mosek: 10.2.1 (recommended SDP solver)

## Usage

This section explains how to use the provided scripts to test the algorithms and reproduce the corresponding plots from the paper. You are free to adjust the hyperparameters directly in the testing scripts to explore different configurations.

#### Certifying Nonnegativity of Polynomials

1. Navigate to the Folder:
   - Open the folder `certifying_nonnegativity`.
2. Run the Testing Scripts:
   - Execute the scripts `test_nonsos.m` and `test_nonsos_family.m` to test the algorithm.
   - The scripts will save the sequence of lower and upper bounds produced by the algorithm.
3. Reproduce Plots from the Paper:
   - Open the script `plot_nonnegativity.m`, specify the name of the polynomial you tested, the number of variables in the polynomial, the type of initialization used in the testing algorithm (`init = 1` represents the first type, `init = 0` represents the second type).
   - Execute the script. The plot will be saved as a PDF file in the same folder.

#### Certifying Copositivity of Matrices

1. Navigate to the Folder:
   - Open the folder `certifying_copositivity`.
2. Run the Testing Scripts:
   - Execute the script `test_standard_qp.m` to test the algorithm.
   - The scripts will save the sequence of lower and upper bounds produced by the algorithm.
3. Reproduce Plots from the Paper:
   - Open the script `plot_copositivity.m`, specify the test index of the matrix.
   - Execute the script. The plot will be saved as a PDF file in the same folder.

#### Computing Clique Number of Graphs

1. Navigate to the Folder:
   - Open the folder `clique_number`.
2. Run the Testing Scripts:
   - Execute the script `test_random_graph.m` to test the algorithm.
   - The scripts will print out the final lower and upper bounds obtained by the algorithm.

## Citation

If you use this code in your research, please cite:

```bibtex
[Citation information to be added]
```

## License

This project is licensed under the MIT License.

## Contact

Yixuan Hua - yh7422@princeton.edu

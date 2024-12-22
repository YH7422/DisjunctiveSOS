# Disjunctive Sum of Squares (DiSOS)

This repository implements the disjunctive sum of squares (DiSOS) method for certifying nonnegativity of polynomials. Unlike traditional sum of squares methods that use a single algebraic identity, DiSOS method uses multiple algebraic identities with the same degree as the polynomial being verified. The implementation integrates branch-and-bound scheme and includes other applications like certifying copositivity of matrices and computing the clique number of graphs.

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
│   ├── runtime_nonneg.csv      # Summary of running times for testing scripts
│   └── plot_nonnegativity.m    # Script for producing plots from the paper
│
├── certifying_copositivity/
│   ├── DiSOS_copositive_BnB.m  # DiSOS implementation for matrix copositivity
│   ├── test_standard_qp.m      # Test script for standard quadratic programming
│   ├── runtime_copos.csv       # Summary of running times for testing scripts
│   └── plot_copositivity.m     # Script for producing plots from the paper
│
└── clique_number/
    ├── DiSOS_clique_BnB.m      # DiSOS implementation for clique number computation
    ├── test_Paley_graph.m      # Test script for Paley graphs
    ├── runtime_clique.csv      # Summary of running times for testing scripts
    └── Paley_graph/              # Subfolder containing graph data
```

## Features

- Implementation of Disjunctive Sum of Squares (DiSOS) method
- Branch and Bound algorithm integration
- Applications in:
  - Polynomial nonnegativity certification
  - Matrix copositivity certification
  - Graph clique number computation

## Requirements

The following versions were used in our experiments. Other versions might work but haven't been tested.

### MATLAB Dependencies

- MATLAB: 9.12 (R2022a)
- CVX: 2.2
- YALMIP: [R20230622](https://github.com/yalmip/YALMIP/releases/tag/R20230622)
- Mosek: 10.2.1 (recommended SDP solver)

### Optional Python Dependencies

For graph generation (optional, as graph files are already included):

- Python: 3.9.13
- Numpy: 1.25.1
- Networkx: 3.1

### Usage

This section explains how to use the provided scripts to test the algorithms and reproduce the corresponding plots from the paper. You are free to adjust the hyperparameters directly in the testing scripts to explore different configurations.

------

#### Certifying Nonnegativity of Polynomials

1. **Navigate to the Folder**:
   - Open the folder `certifying_nonnegativity`.
2. **Run the Testing Scripts**:
   - Execute the scripts `test_nonsos.m` and `test_nonsos_family.m` to test the algorithm.
   - The scripts will save the sequence of lower and upper bounds produced by the algorithm.
   - Running times have been integrated into a summary CSV file located at `certifying_nonnegativity/runtime_nonneg.csv`.

3. **Reproduce Plots from the Paper**:

   - Open the script `plot_nonnegativity.m`, specify the name of the polynomial you tested, the number of variables in the polynomial, the type of initialization used in the testing algorithm (`init = 1` represents the first type, `init = 0` represents the second type).

   - Execute the script. The plot will be saved as a PDF file in the same folder.

#### Certifying Copositivity of Matrices

1. **Navigate to the Folder**:
   - Open the folder `certifying_copositivity`.
2. **Run the Testing Scripts**:
   - Execute the script `test_standard_qp.m` to test the algorithm.
   - The scripts will save the sequence of lower and upper bounds produced by the algorithm.
   - Running times have been integrated into a summary CSV file located at `certifying_copositivity/runtime_copos.csv`.

3. **Reproduce Plots from the Paper**:
   - Open the script `plot_copositivity.m`, specify the test index of the matrix.
   - Execute the script. The plot will be saved as a PDF file in the same folder.

#### Computing Clique Number of Graphs

1. **Navigate to the Folder**:
   - Open the folder `clique_number`.
2. **Run the Testing Scripts**:
   - Execute the script `test_Paley_graph.m` to test the algorithm.
   - The scripts will print out the final lower and upper bounds obtained by the algorithm.
   - Running times have been integrated into a summary CSV file located at `clique_number/runtime_clique.csv`.

## Citation

If you use this code in your research, please cite:

```bibtex
[Citation information to be added]
```

## License

This project is licensed under the MIT License.

## Contact

Yixuan Hua - yh7422@princeton.edu

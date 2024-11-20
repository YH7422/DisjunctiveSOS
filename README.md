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
│   └── test_nonsos_family.m    # Test script for a familiy of non-SOS polynomials
│
├── certifying_copositivity/
│   ├── DiSOS_copositive_BnB.m  # DiSOS implementation for matrix copositivity
│   └── test_standard_qp.m      # Test script for standard quadratic programming
│
└── clique_number/
    ├── DiSOS_clique_BnB.m      # DiSOS implementation for clique number computation
    ├── test_Paley_graph.m      # Test script for Paley graphs
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
- numpy: 1.25.1
- networkx: 3.1

## Citation

If you use this code in your research, please cite:

```bibtex
[Citation information to be added]
```

## License

This project is licensed under the MIT License.

## Contact

Yixuan Hua - yh7422@princeton.edu

Project Link: https://github.com/YH7422/DisjunctiveSOS
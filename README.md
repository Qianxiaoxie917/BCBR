# BCBR: BiConvex Blockwise Regularization for Precision Matrix Estimation

This repository contains the code implementations and related materials for the paper:

> **Regularized estimation of precision matrix for high-dimensional multivariate longitudinal data.**  
> Fang, Qian, Chen Yu, and Zhang Weiping.  
> Journal of Multivariate Analysis 176 (2020): 104580.

## Table of Contents

- [Installation](#installation)
- [File Descriptions](#file-descriptions)
- [Usage](#usage)
- [Contributing](#contributing)
- [Citation](#citation)
- [Contact](#contact)

## Installation

To use the R scripts, ensure you have R installed and set up on your machine. If you wish to utilize the C++ implementation, make sure you have the Rcpp package available in your R environment.

## File Descriptions

1. **mainfun.cpp**: Primary C++ implementation of the BCBR algorithm.
2. **mainfun.R**: R adaptation of the BCBR algorithm from `mainfun.cpp`.
3. **otherfun.R**: Auxiliary functions supporting the BCBR algorithm.
4. **test.R**: Demonstrative runs and test cases for the BCBR algorithm.
5. **CV_CP.R**: Parallelized cross-validation script for the BCBR algorithm.

## Usage

1. Source the `otherfun.R` before using `mainfun.R`.
2. For direct C++ performance, compile `mainfun.cpp` and link to R using the Rcpp package.
3. Run `test.R` for test cases.
4. Utilize `CV_CP.R` for optimal hyperparameter tuning through cross-validation.

## Contributing

Feedback, bug reports, and pull requests are welcome on this repository. For major changes, please open an issue first to discuss the proposed change.

## Citation

If you utilize this codebase in your work, please cite the original paper:

```bibtex
@article{fang2020regularized,
  title={Regularized estimation of precision matrix for high-dimensional multivariate longitudinal data},
  author={Fang, Qian and Yu, Chen and Weiping, Zhang},
  journal={Journal of Multivariate Analysis},
  volume={176},
  pages={104580},
  year={2020}
}

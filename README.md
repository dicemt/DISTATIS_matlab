# DISTATIS_matlab
Author: Daisuke Matsuyoshi (National Institute for Radiological Sciences ([QST-NIRS](https://www.qst.go.jp/site/qst-english/)) and [Araya, Inc.](https://www.araya.org/))

Implementation of DISTATIS written in MATLAB. DISTATIS is a generalization of classial multidimensional scaling and handles *m* sets of data tables collected on the same set of *n* observations, combining datasets into a common consensus structure, called "compromise."


# Usage

## Basic usage
```matlab
[eigval,eigvector,fscore,eigval3,eigvector3,fscore3] = distatis(data);
```

- `data` should be a 3D distance (dissimilarity) matrix [n x n x m].
    - *n x n* distance matrices stacking along Z-direction *m* times
- For example usage, see `distatis_demo.m` 

## Outputs
- eigval         - Eigenvalues for compromise
- eigvector      - Eigenvector for compromise
- fscore         - Factor score for compromise
- eigval3        - Eigenvalues for dim 3
- eigvector3     - Eigenvector for dim 3
- fscore3        - Factor score for dim 3

# Reference
- Abdi H, Toole AJO, Valentin D, Edelman B (2005) [DISTATIS: The analysis of multiple distance matrices.](https://ieeexplore.ieee.org/document/1565340) IEEE CVPR 2005, 42-47. doi: 10.1109/CVPR.2005.445
- Abdi H, Valentin D, Chollet S, Chrea C (2007) [Analyzing assessors and products in sorting tasks: DISTATIS, theory and applications.](https://www.sciencedirect.com/science/article/abs/pii/S0950329306001236) Food Qual Prefer 18:627-640. doi: 10.1016/j.foodqual.2006.09.003

# License
The DISTATIS_matlab is free but copyright software, distributed under the terms of the [BSD 3-Clause "New" or "Revised" License](https://choosealicense.com/licenses/bsd-3-clause/).

Copyright (C) 2020-2022 Daisuke Matsuyoshi

This program is free software: you can redistribute it and/or modify it under the terms of the BSD 3-Clause "New" or "Revised" License.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the BSD 3-Clause "New" or "Revised" License for more details.

You should have received a copy of the BSD 3-Clause "New" or "Revised" License along with this program. If not, see <https://opensource.org/licenses/BSD-3-Clause>.

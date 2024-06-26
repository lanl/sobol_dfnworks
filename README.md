Code for Experiment in Publication “Sensitivity Analysis in the Presence of Intrinsic Stochasticity for Discrete Fracture Network Simulations”
========
[![arXiv](https://img.shields.io/badge/arXiv-2312.04722-GREEN.svg)](https://arxiv.org/abs/2312.04722)

All code used to perform the experiment in the cited publication.  For questions, issues, or clarifications please reach out to Murph: <murph@lanl.gov>.

This repository is organized according to the different directories needed to perform the experiment in the paper.  These directories are:

* dfnworks_drivers: code to perform the DFN simulations to create the data used in the experiment;

* test_data: code mostly identical to the code in dfnworks_drivers.  Used to created a separate testing data set for validation of the methods;

* model_comparisons: compares the joint GP emulator against a joint GAM emulator;

* eda: scripts to create data visualizations on the data created in dfnworks_drivers, sequential_design_10th_percentile, and test_data;

* gp_analysis: code to analyze the Gaussian Process (GP) fit during the Sequential Design in sequential_design_10th_percentile;

* sobol_indices: code to estimate Sobol’ Indices using the GP fit during the Sequential Design in sequential_design_10th_percentile.

## Note on Reproducibility

We suggest cloning this repository into a folder entitled GitHub in your home directory.  Filepaths within this directory structure should then work out of the box.

The figures in the paper can be recreated by running the corresponding scripts in these directory (for instance, figure_4.R creates the visual in Figure 4).  All the data from dfnWorks is pre-compiled into .Rdata files.  This is due to the very large size of a typical dfnWorks output; we compiled this data using scripts in this repository, removing extraneous output for space considerations.  The raw output can be produced using the dfnWorks software suite, using the files in dfnworks_drivers and test_data.  The corresponding seeds a specific inputs can be found in these directories.

To recreate the entire experiment from these files, first run all the figure_##.R files in eda, then these same files in gp_analysis, then these same files in sobol_indices.  We suggest using an HPC environment for the files in this last directory.

If you are attempting to run this experiment from our Github repository, we suggest getting the data files from our Zenodo repository.  Github has insufficient space for the data, so these data are not present.


## Software Required
To create underground particle transport simulation data, one will need access to the dfnWorks simulation suite, available for download [here](https://dfnworks.lanl.gov/). The dfnWorks data files are included in this repository.

## Citation
Alexander C. Murph, Justin D. Strait, Kelly R. Moran, Jeffrey D. Hyman, Hari S. Viswanathan, & Philip H. Stauffer. (2023). Sensitivity Analysis in the Presence of Intrinsic Stochasticity for Discrete Fracture Network Simulations.  _In Review._ [Preprint.](https://arxiv.org/abs/2312.04722)

## Release

This software has been approved for open source release and has been assigned **O4712** 

## Copyright

© 2023. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.

## License

This code repository is distributed under a MIT License:

Copyright 2023
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


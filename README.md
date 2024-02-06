# Project Name

Solving Euler Equations for a Quasi 1D converging diverging nozzle. 

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Contributing](#contributing)
- [License](#license)

## Installation

1. The [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) matrix library needs a variable called `EIGEN3_ROOT` to be set to the top-level Eigen directory
2. [PETSc](http://www.mcs.anl.gov/petsc/) for sparse linear solvers; needs `PETSC_DIR` and `PETSC_ARCH` set.

To build the Doxygen documentation, type the following command in the doc/ directory:

		doxygen Doxyfile.cfg

This requires you to have [Doxygen](http://www.stack.nl/~dimitri/doxygen/index.html) installed. It will generate HTML documentation, which can be accessed through the generated documentation folder
## Usage

Issue the following commands from the top level directory of this repo: 

        make
        ./exec

The solutions are present in `./data` once the solution converges. 



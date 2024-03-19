# Hankel Projection

A package for denosing imaginary time Green's functions with Hankel projections.

## Prerequisites

* Eigen3

## Installation

1. Clone the repository:
```
git clone git@github.com:CQMP/HankelProjection.git
```
2. Build the package:
```
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=PATH_TO_INSTALL/HankelProj && cmake --build build
```
3. Install the package (optional):
```
cmake --build build --target install
```

## Usage

After building the package, refer to the examples directory for sample usage. Execute the command `bash run.sh` within the examples directory to see HankelProj in action.

## Contact

If you have questions, bug reports, or feature suggestions, please feel free to email us at: umyangyu@umich.edu or egull@umich.edu


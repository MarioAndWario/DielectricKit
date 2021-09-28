# DielectricKit

High-performance computing toolkit to calculate and visualize polarizability and
dielectric response function within the random-phase approximation.

## Introduction

`DielectricKit` includes three Fortran programs: `Chi.x`, `EpsInv.x`, and `RealSpace.x`. Some of the libraries and modules are
incorporated from the open-source [`BerkeleyGW`](https://berkeleygw.org) package. The input and output formats are also compatible
with `BerkeleyGW`.

`Chi.x` calculates the polarizability function (`chimat.h5` and/or `chi0mat.h5`) in reciprocal space using Kohn-Sham eigenstates and eigenvalues from a density-functional
theory calculation. Current, we only support [`Quantum ESPRESSO`](https://www.quantum-espresso.org), which is an open-source DFT code.

`EpsInv.x` calculates use `chimat.h5` and/or `chi0mat.h5` as input to calculate the inverse dielectric response function (`epsmat.h5` and/or `eps0mat.h5`) in reciprocal space.

`RealSpace.x` perform fast Fourier transform to calculate polarizability or inverse dielectric response functions in real space. The result with one fixed coordinate is output in the Xcrysden format (.xsf).

## Theoretical formalism

See the pdf file within the `doc` folder.

## Libraries

* Parallelization schemes: hybrid MPI / OpenMP
* Parallel IO: HDF5
* Math libraries:
  * Linear algebra: BLAS, LAPACK, SCALAPACK
  * Fast Fourier transform: FFTW

## Installation

* Copy one `arch.mk` file from `config` folder to `src` folder. Make necessary modifications for the compiler, compilation flags, and library paths.

* Use `make` to compile the source code.  Use `make -j` to enable parallel compilation to same time.

## Usage

### Generate wavefunction from Quantum ESPRESSO

### Calculate polarizability with Chi.x

### Calculate inverse dielectric response function with EpsInv.x

### Plot polarizability or dielectric response function in real space with RealSpace.x

## References

1. Giannozzi, P., et al.QUANTUM ESPRESSO: A modular and open-source software project for quantumsimulations of materials.J. Phys.: Condens. Matter 21, 395502 (2009).

2. Hybertsen, M. S. & Louie, S. G. Electron correlation in semiconductors and insulators: Band gaps andquasiparticles energies.Phys. Rev. B 34, 5390 (1986).

3. Deslippe, J. et al. BerkeleyGW: A massively parallel computer package for the calculation of the quasiparticle and optical properties of materials and nanostructures. Comput. Phys. Commun. 183, 1269-1289 (2012).

4. Adler, S. L. Quantum Theory of the Dielectric Constant in Real Solids. Phys. Rev. 126, 413-420 (1962).

5. Wiser, N. Dielectric Constant with Local Field Effects Included. Phys. Rev. 129, 62-69 (1963).

## Contact

Meng Wu

wu1meng2@berkeley.edu
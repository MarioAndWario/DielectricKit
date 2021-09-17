# EpsInv
  Version 0.1  (2020)

## Introduction


### Dielectric response function


### Installation


### Usage

#### Command line

`EpsInv.x chimat.h5 > EpsInv.out`

or 

`EpsInv.x chi0mat.h5 > Eps0Inv.out`

#### Input:

`chimat.h5`: irreducible polarizability matrix calculated on a uniform k-grid

`chi0mat.h5`: irreducible polarizability matrix calculated on a shifted k-grid with $\Gamma$ shifted to $q_0$

#### Output:

`eps0mat.h5`: inverse dielectric response function matrix calculated on a uniform k-grid

`epsmat.h5`: inverse dielectric response function matrix calculated on a shifted k-grid with $\Gamma$ shifted to $q_0$

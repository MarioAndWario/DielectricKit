# Monolayer MoS2 with spin-orbit coupling (SOC)
============

This example shows how to calculate RPA polarizability and inverse dielectric
response functions of monolayer MoS2 with SOC.

1. Enter the QE folder.

a. Run script_1 to link Quantum ESPRESSO input and output files
b. Run script_2 (with interactive run or batch job) to run Quantum ESPRESSO scf calculations.
   Estimated time: ~20 seconds for 1 Haswell node (32 cpus) at NERSC Cori.
b. Run script_3 (with interactive run or batch job) to run Quantum ESPRESSO non-scf calculations.
   Estimated time: ~4 minutes for 2 Haswell nodes (64 cpus) at NERSC Cori.

   If you want to clean the output, run script_2_3_clean.

2. Go back to the directory containing this file

   Run script_4 to do the following:
a. Calculate polarizability function matrix within Chi folder.
   Estimated time: ~100 seconds for 1 Haswell node (32 cpus) at NERSC Cori.

b. Calculate inverse dielectric response function matrix within EpsInv folder.
   Estimated time: ~13 second for 1 Haswell node (32 cpus) at NERSC Cori.

c. Plot real-space with 2 different r2 position2 within RealSpace folder.
   Estimated time: ~10 seconds for 1 Haswell node (32 cpus) at NERSC Cori.

   If you want to clean the output, run script_4_clean.

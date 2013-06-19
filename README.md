nsolve
======

nsolve is a small Fortran code for finding energy eigenvalues and eigenfunctions of the 1-d Schrodinger equation using the Numerov method. This includes solving the radial equation for 3-d problems.

Contents:

  solve.f90: Main code for doing the Numerov method and finding eigenvalues & wavefunctions
  
  ho.f90: Example program for the 3d harmonic oscillator (radial equation)
  
  pt.f90: Example program for two particles in a harmonic trap interacting via a Poschl-Teller (inverse cosh^2) potential
  
Notes:

  - This code enhances the "standard" method (that I've seen implemented elsewhere) through the use of a high-order numerical derivatve. This works much better.
  
  - Includes an additional code, findpt.f90, which determines the interaction strength of the Poschl-Teller potential that reproduces a given desired ground-state eigenvalue.
  
  - Has been written to easily find eigenfunctions with exactly "k" nodes, where k is a given integer.

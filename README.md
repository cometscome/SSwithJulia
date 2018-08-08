# SSwithJulia
[08/08/18] The code works with Julia 0.7.0-rc3.0. 
Sakurai-Sugiura method to obtain the eigenvalues located in a given domain with Julia 0.6.4. See, T. Sakurai and H. Sugiura: J. Comput. Appl. Math. 159 (2003) 119. and "Numerical Construction of a Low-Energy Effective Hamiltonian in a Self-Consistent Bogoliubov–de Gennes Approach of Superconductivity", Yuki Nagai et al., J. Phys. Soc. Jpn. 82, 094701 (2013) or arXiv:1303.3683 

In SStest.ipynb, the eigenvalues are obtained in the 2D tight-binding model with the superconducting-normal-superconducting π-junction. The shiftedCG method is used.

The calculation speed is not optimized. There might be other good solvers for solving linear equations. 

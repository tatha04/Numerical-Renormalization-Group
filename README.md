# Numerical-Renormalization-Group

## Introduction

The Numerical Renormalization Group (NRG) Technique is a non-perturbative numerical method that was primarily developed to solve the Kondo problem. The Kondo problem deals with a situation that arises in metals due to the interaction of the conduction electrons with spin 1/2 magnetic impurities (arising from unpaired electrons in the d or f orbitals). Any perturbative treatment of this problem exhibits logarithmic divergences at characteristic temperatures. The NRG was able to solve this problem to low temperatures and capture the formation of the singlet ground state. Since then, the NRG has been used in a range of applications including magnetic hosts in metallic, semi-metallic and superconducting hosts, quantum dots, heavy fermion systems and quantum phase transitions. 

For references see
* [Kondo effect](http://www.scholarpedia.org/article/Kondo_effect) (Scholarpedia)
* [Reviews of Modern Physics](https://doi.org/10.1103/RevModPhys.80.395) paper by Bulla et. al.
* [The Kondo Problem to Heavy Fermions](https://doi.org/10.1017/CBO9780511470752) by A. C. Hewson.
* [Computational Physics Blog](https://compphys.go.ro/the-numerical-renormalization-group/) by Adrian Roman.

## Method
The NRG consists of the following key steps:
* Partition of the conduction band into logarithmic bins.
* Mapping the conduction band into a semi-infinite tight binding fermionic chain with nearest neighbor hopping, known as the Wilson chain (WC). The impurity is attached to one end of the WC.
* The hopping coeffecients drop off exponentially, ensuring convergence of the ground state(s).
* Iterative diagonalization of the WC, where an additional site of the WC is added in each iteration. The basis states of the current iteration are formed using the eigenstates of the previous iteration and the basis states of the additional site of the WC.

## About the code
This is a relatively simple code in MATLAB that implements the NRG for a flat band (constant density of states). The main function that runs the code is NRG(). For usage details on the code and usage, please see the file 'Readme.txt'.



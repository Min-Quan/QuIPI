# Inverse iteration quantum eigensolvers assisted with a continuous variable
This repository contains the demonstration code of the article "[Inverse iteration quantum eigensolvers assisted with a continuous variable](https://iopscience.iop.org/article/10.1088/2058-9565/ac5b30)" ([arXiv:2010.03236](https://arxiv.org/abs/2010.03236)). This simulation is based on package [*QuTip*](https://github.com/qutip).

## Introduction
Our quantum algorithm is proposed to solve the eigenstate problem of a given local Hamiltonian. The principle of this algorithm is iteratively performing the inverse Hamiltonian onto a quantum state, so that it converges to an approximate ground state of the Hamiltonian. Here, the inverse Hamiltonian is constructed by integrate of unitaries with the help of continuous-variable qumode.

## Dependecies
- Python 3.7
- Numpy
- Qutip 4.6
- Scipy
- matplotlib

## File description
- [H_generator.py](https://github.com/Min-Quan/QuIPI/blob/main/H_generator.py): a package that can generator Hamiltonian of Hydrogen molecular, transverse field Ising model, and Kitaev model. And it also can generate the local operators of the Ising model, which can be used to apply the Hamiltonian with the first order Trotter decomposition.
- [QuIPI_circuit.py](https://github.com/Min-Quan/QuIPI/blob/main/QuIPI_circuit.py): a package that performs the proposed algorithm. It gives the function that 1. solve the approximate ground state by applying ideal inverse Hamiltonian; 2. evolve the state by our quantum algorithm, where the unitary operator is directly applied; 3. evolve the state by our quantum algorithm with the first order Trotter decomposition under the noiseless or noisy environment.
- [Demon_of_H2_bond_dissociation.ipynb](https://github.com/Min-Quan/QuIPI/blob/main/Demon_of_H2_bond_dissociation.ipynb): the demonstration of solving the ground state energy of Hydrogen molecular with different bond distance.
- [Demon_of_kitaev_ring.ipynb](https://github.com/Min-Quan/QuIPI/blob/main/Demon_of_kitaev_ring.ipynb): the demonstration of solving the ground state energy of Kitaev ring with different parameter.
- [Demin_of_Ising_model_with_Trotter.ipynb](https://github.com/Min-Quan/QuIPI/blob/main/Demin_of_Ising_model_with_Trotter.ipynb): the demonstration of solving the ground state energy of tranverse field Ising model while considering the first Trotter decomposition.

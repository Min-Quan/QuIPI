# Inverse iteration quantum eigensolvers assisted with a continuous variable
This repository contains the demonstration code of the article "[Inverse iteration quantum eigensolvers assisted with a continuous variable](https://iopscience.iop.org/article/10.1088/2058-9565/ac5b30)" ([arXiv:2010.03236](https://arxiv.org/abs/2010.03236)). This simulation is based on package [*QuTip*](https://github.com/qutip).

## Introduction
Our quantum algorithm is proposed to solve the eigenstate problem of a given local Hamiltonian. The principle of this algorithm is iteratively performing the inverse Hamiltonian onto a quantum state, so that it converges to an approximate ground state of the Hamiltonian. Here, the inverse Hamiltonian is constructed by integrate of unitaries with the help of continuous-variable qumode.

## Dependecies
- Python 3.7
- Numpy
- Qutip 4.6
- Scipy


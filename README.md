# Inverse iteration quantum eigensolvers assisted with a continuous variable
This repository contains the demonstration code of the article "[Inverse iteration quantum eigensolvers assisted with a continuous variable](https://iopscience.iop.org/article/10.1088/2058-9565/ac5b30)" ([arXiv:2010.03236](https://arxiv.org/abs/2010.03236)). This simulation is based on package [*QuTip*](https://github.com/qutip).

## Introduction
Our quantum algorithm is proposed to solve the eigenstate problem of a given local Hamiltonian. The principle of this algorithm is iteratively performing the inverse Hamiltonian onto a quantum state, so that it converges to an approximate ground state of the Hamiltonian. Here, the inverse Hamiltonian is constructed by integrate of unitaries with the help of continuous-variable qumode.

### Principle
- Inverse powered iteration method: $|b \rangle^{(k+1)} = \frac{\hat{H} |b\range^{(k)}}{|\hat{H} |b\range^{(k)}|}$. The evolved state converges to the ground state of Hamiltonian $\hat{H}$ while the initial state $|b\rangle^{(0)}$ has a non-zero overlap to the ground state.
- Fourier transformation of the inverse Hamiltonian: $\hat{H}^{-1} = i \int_{0}^{\infty} e^{-i \hat{H} p} dp$
- Qumode assistance implementation of inverse Hamiltonian: $|b\rangle^{(k+1)} = -i \hat{H}^{-1} |b\rangle^{(k)} \propto \langle q=0|e^{-i \hat{H} \hat{p}}|R\rangle |b\rangle^{(k)}$ with $|R\rangle = \int_{0}^{\infty} |p\rangle dp$.

## Dependecies
- Python 3.7
- Numpy
- Qutip 4.6
- Scipy
- matplotlib

## File description

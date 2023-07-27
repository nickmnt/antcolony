# Ant Colony Optimization (ACO) in C++

## Overview
Welcome to the Ant Colony Optimization (ACO) repository! This project provides an efficient implementation of the Ant Colony Optimization algorithm in C++. The ACO algorithm is a powerful metaheuristic used to find approximate solutions for combinatorial optimization problems, inspired by the foraging behavior of ants.

This repository contains both serial and parallel implementations of the ACO algorithm. The parallel version utilizes OpenMP to take advantage of multi-core processors and accelerate the optimization process. This makes it particularly suitable for solving large-scale problems that require substantial computational resources.

## Features
- Serial and parallel (OpenMP) implementations of the Ant Colony Optimization algorithm
- Customizable parameters to fine-tune the algorithm for different problem domains
- Easily extensible for different combinatorial optimization problems
- Efficient handling of large-scale problems through parallelization

## How to run
Use the following Windows Powershell commands.
### Serial version
```console
g++ -fopenmp serial1.cpp -o serial1.exe; ./serial1.exe
```
### Parallel version
```console
g++ -fopenmp ompserial1.cpp -o ompserial1.exe; ./ompserial1.exe
```

## License

This project is licensed under the [MIT License](LICENSE), which means you can use, modify, and distribute the code freely. However, we would appreciate it if you provide attribution to this repository when using it in your projects.

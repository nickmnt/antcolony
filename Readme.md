# How to run
Use the following Windows Powershell commands.
## Serial version
```console
g++ -fopenmp serial.cpp -o serial.exe; ./serial.exe
```
## Parallel version
```console
g++ -fopenmp omp.cpp -o omp.exe; ./omp.exe
```
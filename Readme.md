# How to run
Use the following Windows Powershell commands.
## Serial version
```console
g++ -fopenmp serial1.cpp -o serial1.exe; ./serial1.exe
```
## Parallel version
```console
g++ -fopenmp ompserial1.cpp -o ompserial1.exe; ./ompserial1.exe
```
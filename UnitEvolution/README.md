# Multiprocess unitary evolution via scalapack routines

Building:
---
```shell
mkdir build
cd build
cmake ..
make
mpirun -np 4 ./unitEvolution 16 8 5
```

Command line format:
```
mpirun -np P ./unitEvolution N localN numSteps
mpirun -np P ./unitEvolution N localN numSteps
```
numSteps = number of evolution steps

It is also available to modify input binary matrix files "ro" and "H".

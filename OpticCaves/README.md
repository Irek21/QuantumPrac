# Multiprocess building of distributed Hamiltonian.

423 group Gu Yu, Kamalov Irek. Variant 1.

Building:
---
```shell
mkdir Build
cd Build
cmake ..
make
mpirun -np 4 ./opticCaves 3 3 1 3
```

Command line format:
```
mpirun -np P ./opticCaves numQBits numSteps Emin Emax
mpirun -np P ./opticCaves numQBits numSteps Emin Emax
```
numQBits = number of atoms in system which model N-qbit system

numSteps = number of evolution steps

Emin = minimum excited atoms in system

Emax = maximum excited atoms in system

Vectors a, w, phi are read from files. phi is nomralized after read. It available to modify them in Build folder.

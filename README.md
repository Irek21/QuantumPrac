# QuantumPrac

Assembly:
---
```bash
mkdir build

cd build

cmake ..

make

mpirun -np 4 ./matMul 100 50 z
```

Command line format:
```
mpirun -np P ./matMul N localN z
mpirun -np P ./matMul N localN d
```
(in case of complex values or real)

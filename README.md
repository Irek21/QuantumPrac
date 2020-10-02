# QuantumPrac

Assembly:
---
*mkdir build
*cd build
*cmake ..
*make
*mpirun -np 4 ./matMul 100 50 z

*Command line format: mpirun -np P ./matMul N localN z/d (in case of complex values or real)

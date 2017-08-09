#  MPI Fault Tolerance Library
## Feature list
+ **C/C++**
+ **MPI 3.0**
+ [**User-level checkpoint**](https://github.com/54markov/mpi_fault_tolerance/tree/master/src/ulcp_lib "link to source files")
+ [**ULFM**](http://fault-tolerance.org/category/ulfm/ "official site ULFM")

## Test Samples
+ [**head_2d**](https://github.com/54markov/mpi_fault_tolerance/tree/master/tests/heat_2d "link to source files") - [Laplace equation](https://en.wikipedia.org/wiki/Laplace%27s_equation "wiki Laplace equation") solver by [Jacobi iteration method (https://en.wikipedia.org/wiki/Jacobi_method "wiki Jacobi iteration method")
+ [**n_body**](https://github.com/54markov/mpi_fault_tolerance/tree/master/tests/nbody "link to source files") - an [n-body simulation](https://en.wikipedia.org/wiki/N-body_simulation "wiki N-body simulation") approximates the motion of particles, often specifically particles that interact with one another through some type of physical forces.

## User-level checkpoint library
+ **Rollback Recovery** - checkpoint/restart based 
+ **Failure	detection** - ULFM based
+ **Snapshot creation** - hard drive based (in place/ via NFS)
+ **Compress procedure** - zlib based
+ **Delta snapshot encoding** - XOR operation based

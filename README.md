# OU Portfolio 

A small sampling of the projects I worked on as part of my undergraduate and
masters degrees from the University of Oklahoma. These span all 5 years of my
education and vary dramatically in complexity and scope. Rather than organize
them by class or year, I decided to structure them according to what field
of CS they are in.

Note: This is very much still a work in progress. As you can imagine, 5 years
of school resulted in quite a few projects; it takes a while to dig through
them all. 

## Cryptography

There are three projects here. The first and simplest is a straight-forward
no-library implementation of the RSA algorithm (limited to 64 bit keys to avoid
the complexity of large-integer math) in C. The other two focus on the promising
applications of lattice based methods, first to the Approximate Subset Sum
problem and second to the NTRU public key encryption (and signature) scheme that
is currently being considered as part of the NIST Post-Quantum Cryptography
standardization process. These two are both implemented in SageMath, an open
source superset of python with a feature set similar to Mathematica, Maple, and
other mathematics-oriented software systems. 

## Parallel Programming

There are three projects here that relate specifically to parallel programming,
though, of course, projects in other sections also make use of parallel
constructs. The simplest project of these three is a simple application of
OpenMP to the parallelization of simple C code, including analysis of the
resulting performance. The other two projects both use MPI to distribute large
matrix computations across a super-computer grid. The first does simple matrix
multiplication, while the second is much more complex. It attempts to find
numerical solutions to the heat equations on a uniform 2 dimensional rod, by
first reformulating the problem as a simple linear recurrence, rewriting that
linear recurrence into a matrix-vector equation, with a matrix whose form is
known to be particularly amenable to parallel solutions. Two different
solution techniques (both specific to the case of tridiagonal matrices), are
applied, and performance results are compared. 


## Numerical Methods

There are three projects here, all three of which are generally smaller than
other projects in this portfolio. The first is a very simple implementation 
of Gaussian elimination in Java. The second two are IPython notebooks that
demonstrate theoretical results from basic numerical analysis: first, that the
Newton-Raphson method for root finding converges much quicker than the simpler
bisection method, especially when you have an analytic form of the derivative.
Second, that while higher-order Runge-Kutta methods converge much more quickly
than the Euler method, they do not provide the stability that a (much more
computationally complex) backwards Euler method provides. 

## Networking

There is just one project here, which demonstrates the use of Mininet to emulate
networks. Mininet is a project that uses Linux kernel namespaces to allow a
developer to emulate network traffic using the full kernel network traffic. In
particular, this allows us to define a particular physical network topology and
experiment with Software Defined Networks. In particular, this project uses two
different techniques to construct a L2 SDN on a given topology: 1) manually
create flow tables in that network using ovs and 2) Use the `pox` OpenFlow
controller to automatically construct L2 routing tables in arbitrary
network topologies. 

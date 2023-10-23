# IPM-Solver

<div align="center">

   *An Interior Point Method based on new class of kernel functions for solving linear optimization*

</div>


IPM-kF-Solver (Interior Point Method based on kernel function is a Matlab library that provides an interior point method based on a new kernel function for solving linear optimization problems of the form:

$$
\begin{align}
\min \text{  } & c^Tx \\
\text{s.t.}               & \text{  }  Ax=b\\
                            & \text{  } 0 \le x \\
\end{align}
$$

 The theoretical analysis of Interior Point Method based on new class of kernel functions for solving linear optimization (IPM-kF-Solver) was developed by Nawel BOUDJELLAL (https://github.com/nawel-boudjellal), Hayet ROUMILI, Djamel BENTERKI and Yacine Slimani (https://github.com/slimany09).

 
IPM-Solver was  implemented by Nawel BOUDJELLAL.

## How to cite IPM-Solver

### In an article
The IPM-Solver is been programmed for a publication of our paper entitiled "A new efficient kernel function with a finite double barrier term for LCCO" in TOP Journal by the authors: Nawel BOUDJELLAL, Hayet ROUMILI, Djamel BENTERKI and Yacine Slimani 

### In link

To mention IPM-Solver in link, use (https://github.com/nawel-boudjellal/IPM-kF-Solver).


## Packages

IPM-Solver contains:

Main program (IPMBoudjellal.m)

Three subprograms(PrimalDualLP.m, pointinit.m, depl.m)

Eight real problems (For other problems use  (https://sparse.tamu.edu/LPnetlib).

## Compilation

To compile the IPM-Solver code,  we perform the following instructions:

1- Execute the file PrimalDualLP.m to obtain the initial point to start our solver.

2- Execute the file IPMBoudjellal.m to obtain the optimal solution of liear problem.

## Solving a problem with IPM-kF-Solver

To solve a real problem by IPM-Solver, it is necessary to  change the nom of real problem in two files (PrimalDualLP.m and IPMBoudjellal.m)




***********************  Matlab programs ***************************

All programs in the file  Matlab_code.zip
are implemented in Matlab  and  correspond to the following chapters
in the Appendix:

A.1: Matlab Programs for Gaussian Elimination using LU
Factorization:  the main program is

Poisson2D_LU.m   (Example 8.2)

and functions which are used by this programs are in files:

DiscretePoisson2D.m      (Example 8.2)
LU_PP.m                 (Algorithm 8.1)
ForwSub.m               (Algorithm 8.2)
BackSub.m               (Algorithm 8.3)


A.2: Matlab programs for Cholesky decomposition: the main program is

Poisson2D_Chol.m          (Example 8.4.4)

and functions which are used by this programs are in files:

DiscretePoisson2D.m    (Example 8.2)
Cholesky.m             (Algorithm 8.10)
ForwSub.m              (Algorithm 8.2)
BackSub.m              (Algorithm 8.3)

A.3: Matlab Programs testing Hager’s condition estimator:
the main program is

TestHagersCondAlg.m  (Example 8.4)

and function which is used by this program is in the file: 

HagersAlg.m   (Algorithm 8.7)


A.4: Matlab Program FitFunctionNormaleq.m (Example 9.2)
to test fitting
to a polynomial using method of normal equations.

A.5: Matlab Program FitFunctionQRCGS.m to test fitting to a
polynomial using QR decomposition via classical Gram-Schmidt (CGS)
orthogonalization procedure (Algorithm 9.4).

A.6: Matlab Program CGS.m performing QR decomposition via
 classical Gram-Schmidt (CGS) orthogonalization procedure  (Algorithm 9.4).

A.7: Matlab Programs to fit a function using linear splines: the
main program is

MainHatFit.m                       (Example 9.3)

and functions which are used by this program are in files:

fihatt.m                     (Example 9.3)
LLSChol.m                    (Algorithm 8.10)
LLSQR.m                      (Algorithm 9.4)
LLSSVD.m
newtonIR.m                   (Algorithm 8.8)


A.8: Matlab Programs to fit a function using bellsplines. The
main program is MainBellspline.m (Example 9.4).


and functions which are used by this program are in files:

LLSChol.m               (Algorithm 8.10)
LLSQR.m                 (Algorithm 9.4)
LLSSVD.m
newtonIR.m              (Algorithm 8.8)

A.9 : Matlab Program PowerM.m (Example 10.1)
to test Power Method (Algorithm 10.1).

A.10: Matlab Program InverseIteration.m (Examples 10.5-10.8)
to test Inverse
Iteration Method (Algorithm 10.2).

A.11: Matlab Program MethodOrtIter.m  (Examples 10.9-10.14)  to test Method of
Orthogonal Iteration  (Algorithm 10.3)


A.12: Matlab Program MethodQR iter.m (Example 10.15) to test Method of
QR Iteration  (Algorithm 10.4).

A.13: Matlab Program MethodQR shift.m (Example 10.16) to test Method of
QR Iteration with Shifts (Algorithm 10.5).

A.14: Matlab Program MethodQR Wshift.m (Example 10.16) to test Method of
QR Iteration with shifts (Algorithm 10.5) using Wilkinson’s Shift.


A.15: Matlab Program HessenbergQR.m  (Example 10.17): first is used
Hessenberg Reduction (Algorithm 10.6) and then the Method of QR
Iteration (Algorithm 10.4).

A.16: Matlab Program  testRayleigh.m (Example 11.1)
for computation
the Rayleigh Quotient (Algorithm 11.1).

Function which is used by the main program testRayleigh 
is in the file:

RayleighQuotient.m              (Algorithm 11.1)

A.17:  Matlab Program  for computation of the
algorithm of Divide-and-Conquer: the main program is 
testDC.m  (Example 11.2)

and function which is used by this program is in the file: 

DivideandConq.m    (Algorithm 11.2)

A.18: Matlab Program Bisection.m (Example 11.3, Algorithm 11.4)
which finds all eigenvalues of the matrix A ion the input interval [a,b).

Function which is used by the main program Bisection.m is in the file:

 Negcount.m

A.19: Matlab Program testClassicalJacobi.m  (Example 11.4).

Function which is used by the main program testClassicalJacobi.m
is in the file:
RunJacobi.m             (Algorithm 11.7)

A.20: Matlab Program testSVDJacobi.m  (Example 11.5)
Function which is used by the main program  testSVDJacobi.m
is in the file:

RunSVDJacobi.m   (Algorithm 11.14)

A.21: Matlab Programs for solution of the Dirichlet problem
for the Poisson's equation in 2D on a square using iterative Jacobi method:
the main program is

Poisson2D_Jacobi.m   (Example 12.1, Algorithms 12.1, 12.2)

and function which is used by this program is in the file:

DiscretePoisson2D.m

A.22: Matlab Programs for solution of the Dirichlet problem
for the Poisson's equation in 2D on a square using iterative
Gauss-Seidel method:
the main program is

Poisson2D_Gauss_Seidel.m    (Example 12.2, Algorithms 12.3)

and function which is used by this program is in the file:

DiscretePoisson2D.m


A.23: Matlab Programs for solution of the Dirichlet problem
for the Poisson's equation in 2D on a square using iterative
Gauss-Seidel method  with red-black ordering:
the main program is

Poisson2D_Gauss_SeidelRedBlack.m  (Example 12.2, Algorithm 12.4)


and function which is used by this program is:

DiscretePoisson2D.m

A.24: Matlab Programs for solution of the Dirichlet problem
for the Poisson's equation in 2D on a square using iterative
  SOR method: the main program is

Poisson2D_SOR.m   (Example 12.3, Algorithms 12.5, 12.6)

and function which is used by this program is in the file:

DiscretePoisson2D.m


A.25: Matlab Programs for solution of the Dirichlet problem for the
Poisson's equation in 2D on a square using Conjugate Gradient method:
the main program is

Poisson2D ConjugateGrad.m    (Example 12.4, Algorithm  12.13)

and function which is used by this program is in the file:

DiscretePoisson2D.m


A.26: Matlab Programs for solution of the Dirichlet problem for the
Poisson's equation in 2D on a square using Preconditioned Conjugate
Gradient method: the main program is

Poisson2D_PrecConjugateGrad.m  (Example 12.6, Algorithm  12.14)

and function which is used by this program is in the file:

DiscretePoisson2D.m


*********************  C++ and PETSc programs ***************************

A.27: PETSc programs for the solution of the Poisson’s equation
in two dimensions  on a square using different iterative methods.

All programs in the file PETSc_code.zip are implemented in C++ and the
software package PETSc (http://www.mcs.anl.gov/petsc/).

These
programs illustrate Example 12.5: solution of the Dirichlet problem
for the Poisson's equation in 2D on a square using different iterative methods.

The different iterative methods 
are encoded by numbers 1-7 in the main program Main.cpp
in the
following order:
1 - Jacobi’s method,
2 - Gauss-Seidel method,
3 - Successive Overrelaxation method (SOR),
4 - Conjugate Gradient method,
5 - Conjugate Gradient method (Algorithm 12.13),
6 - Preconditioned Conjugate Gradient method,
7 - Preconditioned Conjugate Gradient method (Algorithm 12.14).

Methods 1-5 use inbuilt PETSc functions, and
methods 6,7 implement algorithms 12.13, 12.14, respectively. For
example,  we can run the program  Main.cpp using SOR method  as follows:

> nohup Main 3 > result.m

After running the results will
be printed in the file result.m and can be viewed in Matlab using the
command surf(result).


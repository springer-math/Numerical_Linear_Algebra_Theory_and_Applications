#ifndef _CE3_H
#define _CE3_H

#include <petsc.h>

PetscErrorCode CreateMatrix(Mat*, PetscInt, PetscInt);
PetscErrorCode CreateVector(Vec*, PetscInt);

PetscErrorCode DiscretePoisson2D(PetscInt, Mat*);
//PetscErrorCode DiscretePoisson2D_coeffs(PetscInt, PetscScalar, Mat*, Mat*, Vec*);
PetscErrorCode DiscretePoisson2D_coeffs(PetscInt, PetscScalar, Vec*);

/***********************/
/* METHODS OF SOLUTION */
/***********************/
PetscErrorCode Solve(Mat, Vec, Vec, PetscInt, bool);
PetscErrorCode Jacobi(PC);
PetscErrorCode GaussSeidel(PC);
PetscErrorCode SOR(PC);
PetscErrorCode ConjugateGradient(KSP, PC);
PetscErrorCode ConjugateGradient_full(Mat, Vec, Vec, bool);
PetscErrorCode ConjugateGradient_inner(Mat, Vec, Vec, Mat, bool);
PetscErrorCode PreconditionedConjugateGradient(KSP, PC);
PetscErrorCode PreconditionedConjugateGradient_full(Mat, Vec, Vec, bool);
PetscErrorCode PreconditionedConjugateGradient_inner(Mat, Vec, Vec, Mat, bool);

enum SolverMethod {
    METHOD_INVALID=0,
    METHOD_JACOBI=1,
    METHOD_GAUSS_SEIDEL=2,
    METHOD_SOR=3,
    METHOD_CG=4,
    METHOD_CG_FULL=5,
    METHOD_PCG=6,
    METHOD_PCG_FULL=7
};

#endif/*_CE3_H*/

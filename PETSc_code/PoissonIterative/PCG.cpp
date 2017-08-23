
/*Program for using Preconditioned Conjugate gradient method */

#include <petsc.h>
#include <petscmat.h>
#include <petscvec.h>
#include <petscksp.h>
#include <cmath>
#include "Poisson.h"


PetscErrorCode PreconditionedConjugateGradient(KSP ksp, PC preconditioner) {
    PetscErrorCode ierr;

    ierr = KSPSetType(ksp, KSPCG);

    //ierr = PCSetType(preconditioner, PCJACOBI); CHKERRQ(ierr);
    ierr = PCSetType(preconditioner, PCCHOLESKY); CHKERRQ(ierr);

    return 0;
}

/**
 * Implements the preconditioned conjugate gradient
 * method with Jacobi preconditioning.
 */
PetscErrorCode PreconditionedConjugateGradient_full(Mat A, Vec b, Vec x,
						    bool VERBOSE) {
    PetscErrorCode ierr;
    Mat Minv;
    Vec diagonal, unity;
    PetscInt n;

    ierr = MatGetSize(A, &n, NULL); CHKERRQ(ierr);
    ierr = CreateMatrix(&Minv, n, n); CHKERRQ(ierr);
    ierr = CreateVector(&diagonal, n); CHKERRQ(ierr);
    ierr = CreateVector(&unity, n); CHKERRQ(ierr);

    ierr = MatAssemblyBegin(Minv, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(Minv, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(diagonal); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(diagonal); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(unity); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(unity); CHKERRQ(ierr);

    // We use the diagonal preconditioner for simplicity
    ierr = MatGetDiagonal(A, diagonal); CHKERRQ(ierr);

    // Compute inverse of all diagonal entries
    ierr = VecSet(unity, 1.0); CHKERRQ(ierr);
    ierr = VecPointwiseDivide(diagonal, unity, diagonal);

    // Create M^{-1}
    ierr = MatDiagonalSet(Minv, diagonal, INSERT_VALUES); CHKERRQ(ierr);

    return PreconditionedConjugateGradient_inner(A, b, x, Minv, VERBOSE);
}
PetscErrorCode PreconditionedConjugateGradient_inner(Mat A, Vec b, Vec x,
					            Mat Minv, bool VERBOSE) {
    PetscErrorCode ierr;
    PetscInt k=0, n;
    PetscScalar mu, nu, yTr, pTz, rNorm, tol = 1e-12;
    Vec p, r, y, z;

    ierr = MatGetSize(A, &n, NULL); CHKERRQ(ierr);

    CreateVector(&p, n);
    CreateVector(&r, n);
    CreateVector(&y, n);
    CreateVector(&z, n);

    VecCopy(b, r);
    ierr = MatMult(Minv, b, p); CHKERRQ(ierr);
    VecCopy(p, y);

    ierr = VecAssemblyBegin(p); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(p); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(r); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(r); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(y); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(y); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(z); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(z); CHKERRQ(ierr);

    ierr = VecZeroEntries(x);

    // Pre-compute first (y^T r)
    ierr = VecDot(y, r, &yTr); CHKERRQ(ierr);

    do {
        k++;

        // z = A * p_k
        ierr = MatMult(A, p, z); CHKERRQ(ierr);

        // nu_k = y_{k-1}^T r_{k-1} / p_k^T z
        ierr = VecDot(p, z, &pTz); CHKERRQ(ierr);
        nu = yTr / pTz;

        // x_k = x_{k-1} + nu_k p_k
        ierr = VecAXPY(x, nu, p); CHKERRQ(ierr);

        // r_k = r_{k-1} - nu_k z
        ierr = VecAXPY(r, -nu, z); CHKERRQ(ierr);

        // y_k = M^{-1} r_k
        ierr = MatMult(Minv, r, y); CHKERRQ(ierr);

        // y_k^T r_k
        mu = 1 / yTr;
        ierr = VecDot(y, r, &yTr); CHKERRQ(ierr);

        // mu_{k+1}
        mu = yTr * mu;

        // p_{k+1} = r_k + mu_{k+1} p_k
        ierr = VecAYPX(p, mu, y);

        // || r_k ||_2
        ierr = VecNorm(r, NORM_2, &rNorm);
    } while (rNorm > tol);

    if (VERBOSE) {
        PetscPrintf(PETSC_COMM_WORLD, "Number of iterations: %d\n", k);
    }

    return 0;
}

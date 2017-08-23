
/*
 * Program for two versions of the Conjugate gradient method.
 */

#include <petsc.h>
#include <petscmat.h>
#include <petscvec.h>
#include <petscksp.h>
#include <cmath>
#include "Poisson.h"
/**
 * Conjugate gradient method using  inbuilt PETSc functions.
 */

PetscErrorCode ConjugateGradient(KSP ksp, PC preconditioner) {
    PetscErrorCode ierr;

    ierr = KSPSetType(ksp, KSPCG);
    ierr = PCSetType(preconditioner, PCNONE); CHKERRQ(ierr);

    return 0;
}

/**
 * An implementation of the conjugate gradient method
 * not utilizing the PETSc KSP interface, but 
 * implementing the matrix/vector operations directly.
 */
PetscErrorCode ConjugateGradient_full(Mat A, Vec b, Vec x, bool VERBOSE) {
    PetscErrorCode ierr;
    PetscInt k=0, n;
    PetscScalar mu, nu, rTr, pTz, rNorm, tol = 1e-12;
    Vec p, r, z;

    ierr = MatGetSize(A, &n, NULL); CHKERRQ(ierr);

    CreateVector(&p, n);
    CreateVector(&r, n);
    CreateVector(&z, n);

    VecCopy(b, p);
    VecCopy(b, r);

    ierr = VecAssemblyBegin(p); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(p); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(r); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(r); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(z); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(z); CHKERRQ(ierr);

    ierr = VecZeroEntries(x);

    // Pre-compute first (r^T r)
    ierr = VecDot(r, r, &rTr); CHKERRQ(ierr);

    do {
        k++;

        // z = A * p_k
        ierr = MatMult(A, p, z); CHKERRQ(ierr);

        // nu_k = r_{k-1}^T r_{k-1} / p_k^T z
        ierr = VecDot(p, z, &pTz); CHKERRQ(ierr);
        nu = rTr / pTz;

        // x_k = x_{k-1} + nu_k p_k
        ierr = VecAXPY(x, nu, p); CHKERRQ(ierr);

        // r_k = r_{k-1} - nu_k z
        ierr = VecAXPY(r, -nu, z); CHKERRQ(ierr);

        // r_k^T r_k
        mu = 1 / rTr;
        ierr = VecDot(r, r, &rTr); CHKERRQ(ierr);

        // mu_{k+1}
        mu = rTr * mu;

        // p_{k+1} = r_k + mu_{k+1} p_k
        ierr = VecAYPX(p, mu, r);

        // || r_k ||_2
        ierr = VecNorm(r, NORM_2, &rNorm);
    } while (rNorm > tol);

    if (VERBOSE) {
        PetscPrintf(PETSC_COMM_WORLD, "Number of iterations: %d\n", k);
    }

    return 0;
}


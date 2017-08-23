
/* Program for choosing different PETSc    preconditioners. */
//*****************************************************************

#include <petsc.h>
#include <petscmat.h>
#include <petscvec.h>
#include <petscksp.h>
#include <cmath>
#include "Poisson.h"

PetscErrorCode Solve(Mat S, Vec h2b, Vec u, PetscInt method, bool VERBOSE) {
    PetscErrorCode ierr;
    KSP ksp;
    KSPConvergedReason convergedReason;
    PC preconditioner;
    PetscInt number_of_iterations;

    ierr = KSPCreate(PETSC_COMM_WORLD, &ksp); CHKERRQ(ierr);
    ierr = KSPSetOperators(ksp, S, S); CHKERRQ(ierr);
    //ierr = KSPSetOperators(ksp, S, S, DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);

    ierr = KSPGetPC(ksp, &preconditioner); CHKERRQ(ierr);
    if (method == METHOD_JACOBI) {
        ierr = Jacobi(preconditioner); CHKERRQ(ierr);
    } else if (method == METHOD_GAUSS_SEIDEL) {
        ierr = GaussSeidel(preconditioner); CHKERRQ(ierr);
    } else if (method == METHOD_SOR) {
        ierr = SOR(preconditioner); CHKERRQ(ierr);
    } else if (method == METHOD_CG) {
        ierr = ConjugateGradient(ksp, preconditioner); CHKERRQ(ierr);
    } else if (method == METHOD_PCG) {
        ierr = PreconditionedConjugateGradient(ksp, preconditioner); CHKERRQ(ierr);
    }

    ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);

    ierr = KSPSolve(ksp, h2b, u); CHKERRQ(ierr);
    ierr = KSPGetIterationNumber(ksp, &number_of_iterations); CHKERRQ(ierr);

    ierr = KSPGetConvergedReason(ksp, &convergedReason); CHKERRQ(ierr);

    if (convergedReason < 0) {
        PetscPrintf(PETSC_COMM_WORLD,
       "KSP solver failed to converge! Reason: %d\n", convergedReason);
    }

    if (VERBOSE) {
        PetscPrintf(PETSC_COMM_WORLD, "Number of iterations: %d\n", number_of_iterations);
    }

    ierr = KSPDestroy(&ksp); CHKERRQ(ierr);

    return 0;
}

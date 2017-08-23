
/* Program to  create matrix  and vector in PETSc. */

#include <petsc.h>
#include <petscmat.h>
#include <petscvec.h>

PetscErrorCode CreateMatrix(Mat *A, PetscInt rows, PetscInt cols) {
    PetscErrorCode ierr;
    ierr = MatCreate(PETSC_COMM_WORLD, A); CHKERRQ(ierr);
    ierr = MatSetSizes(*A, PETSC_DECIDE, PETSC_DECIDE, rows, cols); CHKERRQ(ierr);
    ierr = MatSetFromOptions(*A); CHKERRQ(ierr);
    ierr = MatSetUp(*A); CHKERRQ(ierr);

    return 0;
}

PetscErrorCode CreateVector(Vec *v, PetscInt N) {
    PetscErrorCode ierr;

    ierr = VecCreate(PETSC_COMM_WORLD, v); CHKERRQ(ierr);
    ierr = VecSetSizes(*v, PETSC_DECIDE, N); CHKERRQ(ierr);
    ierr = VecSetFromOptions(*v); CHKERRQ(ierr);

    return 0;
}



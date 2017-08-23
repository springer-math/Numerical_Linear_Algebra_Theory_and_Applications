

/*  Jacobi's method */


#include <petsc.h>
#include <petscmat.h>
#include <petscvec.h>
#include <petscksp.h>
#include <cmath>

/**
 * Returns the preconditioner used for Jacobi's method
 */
PetscErrorCode Jacobi(PC preconditioner) {
    PetscErrorCode ierr;
    ierr = PCSetType(preconditioner, PCJACOBI); CHKERRQ(ierr);

    return 0;
}

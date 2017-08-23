

/*  Gauss-Seidel method */


#include <petsc.h>
#include <petscmat.h>
#include <petscvec.h>
#include <petscksp.h>
#include <cmath>
#include "Poisson.h"

PetscErrorCode GaussSeidel(PC preconditioner) {
    PetscErrorCode ierr;
    ierr = PCSetType(preconditioner, PCSOR); CHKERRQ(ierr);

    /**
     * To use the Gauss-Seidel method we set
     * omega = 1.
     */
    // By default, omega = 1, so the below line is not necessary
    //ierr = PCSORSetOmega(preconditioner, 1.0); CHKERRQ(ierr);

    return 0;
}


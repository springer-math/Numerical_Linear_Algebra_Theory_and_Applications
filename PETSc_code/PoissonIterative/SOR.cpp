
/* Program  for computing SOR */

#include <petsc.h>
#include <petscmat.h>
#include <petscvec.h>
#include <petscksp.h>
#include <cmath>
#include "Poisson.h"

const PetscScalar omega = 1.5;

PetscErrorCode SOR(PC preconditioner) {
    PetscErrorCode ierr;

    ierr = PCSetType(preconditioner, PCSOR); CHKERRQ(ierr);
    ierr = PCSORSetOmega(preconditioner, omega); CHKERRQ(ierr);

    return 0;
}

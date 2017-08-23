
/* Program for generatation of the discretized Laplacian */

#include <petsc.h>
#include <petscmat.h>
#include <petscvec.h>
#include <cmath>

const PetscScalar A_amplitude = 12.;
const PetscScalar f_amplitude = 1.;
const PetscScalar c_x = 1.;
const PetscScalar c_y = 1.;
const PetscScalar poisson_x0 = 0.5;
const PetscScalar poisson_y0 = 0.5;

/**
 * Compute coefficient matrices.
 *
 * n: Number of rows of matrices
 * h: Timestep length
 * C: n-by-n matrix
 * D: (n*n)-by-(n*n) matrix
 * f:
 **/
PetscErrorCode DiscretePoisson2D_coeffs(PetscInt n, PetscScalar h, Vec *h2b) {
    PetscErrorCode ierr;
    PetscInt i, j, idx2[n*n];
    PetscScalar *vecb = new PetscScalar[n*n];

    // Compute C, D and f
    PetscScalar xarg, yarg, expfunc, a, f;
    for (i = 0; i < n; i++) {
        xarg = (((i+1) * h - poisson_x0)) / c_x;

        for (j = 0; j < n; j++) {
            idx2[i*n + j] = i*n + j;

            yarg = (((j+1) * h - poisson_y0)) / c_y;
            expfunc = exp(-(xarg*xarg/2 + yarg*yarg/2));

            f = f_amplitude * expfunc;
            a = 1 + A_amplitude * expfunc;

            vecb[i*n + j] = h*h * f / a;
        }
    }

    ierr = VecSetValues(*h2b, n*n, idx2, vecb, INSERT_VALUES); CHKERRQ(ierr);

    delete [] vecb;

    return 0;
}
PetscErrorCode DiscretePoisson2D(PetscInt n, Mat *A) {
    PetscErrorCode ierr;
    PetscInt i, k, curr, next, matsize = n*n, idx[matsize];
    PetscScalar *matrep = new PetscScalar[matsize*matsize];

    // Initialize all elements to 0
    for (i = 0; i < matsize; i++) {
        // Create index vectors
        idx[i] = i;

        for (k = 0; k < matsize; k++) {
            matrep[i*matsize + k] = 0;
        }
    }

    // Set main diagonal
    for (i = 0; i < matsize; i++)
        matrep[i*matsize + i] = 4.;

    // 1st and 2nd off-diagonals
    for (k = 0; k < n; k++) {
        for (i = 0; i < n-1; i++) {
            curr = (n*k + i);
            next = (n*k + i + 1);

            matrep[curr*matsize + next] = -1;
            matrep[next*matsize + curr] = -1;
        }
    }

    // 3rd and 4th off-diagonals
    for (i = 0; i < n*(n-1); i++) {
        matrep[i*matsize + (i+n)] = -1;
        matrep[(i+n)*matsize + i] = -1;
    }

    ierr = MatSetValues(*A, matsize, idx, matsize, idx, matrep, INSERT_VALUES);
    CHKERRQ(ierr);

    delete [] matrep;

    return 0;
}

/**
 * @file main.cpp
 * @author theoneinaminion
 * @brief Oh yes, it is all coming together.
 * @version 0.1
 * @date 2023-12-28
 * 
 * @copyright Copyright (c) 2023
 * 
 */


#include "space.hpp"

static char help[] = "I came, I saw, I sighed.";

int main(int argc,char **argv){

    PetscInitialize(&argc,&argv,(char*)0,help);
    PetscErrorCode ierr;

    Mesh msh;
    msh.generate_grid();
    PetscScalar pr = 0.8; // Pressure ratio
    flux flx(msh,pr);
    ierr = flx.initialize_primitives();
    ierr = flx.assemble_conservative_vec();
    ierr = flx.assemble_flux_vec();
    //std::cout << flw.pr << std::endl;
    VecView(flx.f,PETSC_VIEWER_STDOUT_SELF);

    

    return ierr;
    PetscFinalize();

}
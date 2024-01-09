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
#include "time.hpp"

static char help[] = "I came, I saw, I sighed.";

int main(int argc,char **argv){

    PetscInitialize(&argc,&argv,(char*)0,help);
    PetscErrorCode ierr;

    //########## Structure of main.cpp begins ##########
    /*
     Mesh msh; 
     msh.generate-grid();
     PetscScalar pr = 0.8; // Pressure ratio
     flux flx(msh,pr);
     Solver slv(flx);
     slv.solve();

    */
    //########## Structure of main.cpp ends ##########

    Mesh msh;
    msh.generate_grid();
    PetscScalar pr = 0.8; // Pressure ratio
    flux flx(msh,pr);
    Solver slv(&flx);
    ierr = slv.flx->initialize_primitives();
    ierr = slv.flx->assemble_conservative_vec();
    ierr = slv.flx->assemble_flux_vec();
    ierr = slv.flx->assemble_source_vec();
    ierr = slv.flx->element_flux_jacobian(0);
    //std::cout << flw.pr << std::endl;
    // for (int i = 0; i < flx.mesh.nvars; i++){
    //     std::cout << flx.qelem[i] << std::endl;
    // }
    //VecView(flx.q,PETSC_VIEWER_STDOUT_SELF);
    MatView(slv.flx->A,PETSC_VIEWER_STDOUT_SELF);    
    

    return ierr;
    PetscFinalize();
}
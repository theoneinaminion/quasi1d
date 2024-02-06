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

    Mesh msh;
    msh.ngrid = 200; //Changes mesh size
    msh.generate_grid();
    PetscScalar pr = 0.6; // Pressure ratio
    flux flx(msh,pr);
    Solver slv(&flx);
    
    ierr = slv.solve(); CHKERRQ(ierr);
    PetscScalar pmax, rhomax;
    ierr = VecMax(slv.flx->p,NULL,&pmax); CHKERRQ(ierr);
    ierr = VecMax(slv.flx->rho,NULL,&rhomax); CHKERRQ(ierr);
    
    std::cout << "Solution written to file" << std::endl;
    ierr = slv.write_soln(); CHKERRQ(ierr); 
    return ierr;
    PetscFinalize();
}
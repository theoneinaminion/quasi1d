#include "time.hpp"

Solver::Solver(flux* flux):flx(flux)
{
    dt = 1e-6;
    CFL = 10;
    std::cout << "Setting initial time step to:" << dt<<std::endl;
    std::cout << "Setting initial CFL to:" << CFL<<std::endl;

    //Setting up the vectors
    PetscErrorCode ierr;
    ierr = VecCreateSeq(PETSC_COMM_SELF,flx->mesh.ngrid,&res);
    if(ierr)
        std::cout << "Couldn't create residual vector!\n";
    ierr = VecCreateSeq(PETSC_COMM_SELF,flx->mesh.ngrid,&dw);
    if(ierr)
        std::cout << "Couldn't create dw vector!\n";    

    std::cout << "Setting KSP Object"<<std::endl;
    
    ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);
    if(ierr)
		std::cout << "Couldn't Create KSP object!\n";
    
    res_elem = new PetscScalar[flx->mesh.nvars];    

}

Solver::~Solver()
{
    std::cout << "Destroying Solver object"<<std::endl;
    PetscErrorCode ierr;
    ierr = KSPDestroy(&ksp);

    ierr = VecDestroy(&res);
    if(ierr)
        std::cout << "Couldn't destroy residual vector!\n";
    ierr = VecDestroy(&dw);
    if(ierr)
        std::cout << "Couldn't destroy dw vector!\n";
    if (res_elem)
        delete[] res_elem;
  

}

void Solver::new_time_step()
{
    dt = prev_resnrm/resnrm;
    
}

PetscErrorCode Solver::setup_ksp(){

    PetscErrorCode ierr;
    PC pc;
    ierr = KSPSetOperators(ksp,flx->A,flx->A); CHKERRQ(ierr);
    ierr = KSPSetType(ksp,KSPPREONLY); CHKERRQ(ierr);
    ierr = KSPGetPC(ksp,&pc); CHKERRQ(ierr);
    ierr = PCSetType(pc,PCLU); CHKERRQ(ierr);
    ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);
    ierr = KSPSetUp(ksp); CHKERRQ(ierr);
    std::cout << "KSP Object setup complete"<<std::endl;
    return ierr;
}

PetscErrorCode Solver::solve()
{
    /**
     * @brief Non Linear Solver 
     * @param mesh Mesh object
     * @param pr Pressure ratio
     * 
     */
    PetscErrorCode ierr;
    ierr = 0;
    
    //################ Basic structure begins ################
    /*
        flx.initialize_primitives();
        
        ierr =setup_ksp();
        while (resnrm>tol)
        {
            prev_res = res;
            flx.assemble_conservative_vec();
            flx.assemble_flux_vec();
            flx.assemble_source_vec();
            flx.assemble_jacobian();
            flx.assemble_residual();
            KSPSolve()
            updatesoln();
            flx.conservative_to_primitive();
            VecNorm(res,NORM_2,&resnrm);

            write_solution();

        }
    */
   //################ Basic structure Ends ################

   return ierr;
    

}
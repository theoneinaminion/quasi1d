#include "time.hpp"

Solver::Solver(flux* flux):flx(flux)
{
    dt = 1e-6;
    CFL = 10;
    std::cout << "Setting initial time step to:" << dt<<std::endl;
    std::cout << "Setting initial CFL to:" << CFL<<std::endl;

    //Setting up the vectors
    PetscErrorCode ierr;
    int size = flx->mesh.ngrid*flx->mesh.nvars - 2*flx->mesh.nvars;
    ierr = VecCreateSeq(PETSC_COMM_SELF,size,&res);
    if(ierr)
        std::cout << "Couldn't create residual vector!\n";
    ierr = VecCreateSeq(PETSC_COMM_SELF,size,&dw);
    if(ierr)
        std::cout << "Couldn't create dw vector!\n";    

    std::cout << "Setting KSP Object"<<std::endl;
    
    ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);
    if(ierr)
		std::cout << "Couldn't Create KSP object!\n";
    

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

PetscErrorCode Solver::compute_residual()
{
   
   PetscErrorCode ierr;
   PetscInt idx[flx->mesh.nvars];
   PetscScalar sm, sp, vi, fi[flx->mesh.nvars], fp[flx->mesh.nvars],fm[flx->mesh.nvars], wi[flx->mesh.nvars], wp[flx->mesh.nvars], wm[flx->mesh.nvars], q[flx->mesh.nvars],res_elem[flx->mesh.nvars];
   
   ierr = flx->assemble_flux_vec(); CHKERRQ(ierr); //Build the flux vector
   ierr = flx->assemble_source_vec(); CHKERRQ(ierr); //Build the source vector

    for(PetscInt iel=1; iel < flx->mesh.ngrid-1; iel++)
    {
      // Prelimnary Data from the element
      ierr = flx->get_elem_eigen(iel); CHKERRQ(ierr);
      ierr = VecGetValues(flx->mesh.vol,1,&iel,&vi); CHKERRQ(ierr);
      ierr = VecGetValues(flx->mesh.sw,1,&iel,&sm); CHKERRQ(ierr);
      int el = iel+1;
      ierr = VecGetValues(flx->mesh.sw,1,&el,&sp); CHKERRQ(ierr);

      // Get f, w vectors at iel, iel+1, iel-1 and q at iel
      // At iel
        for (int i = 0; i < flx->mesh.nvars; i++)
        {
            idx[i] = iel*flx->mesh.nvars + i; 
        }
        ierr = VecGetValues(flx->f,flx->mesh.nvars,idx,fi); CHKERRQ(ierr);
        ierr = VecGetValues(flx->w,flx->mesh.nvars,idx,wi); CHKERRQ(ierr);
        ierr = VecGetValues(flx->q,flx->mesh.nvars,idx,q); CHKERRQ(ierr);
       
      //At iel+1
        el = iel+1;
        for (int i = 0; i < flx->mesh.nvars; i++)
        {
          idx[i] = el*flx->mesh.nvars + i; 
        }
        ierr = VecGetValues(flx->f,flx->mesh.nvars,idx,fp); CHKERRQ(ierr);
        ierr = VecGetValues(flx->w,flx->mesh.nvars,idx,wp); CHKERRQ(ierr);

      // At iel-1
        el = iel-1;
        for (int i = 0; i < flx->mesh.nvars; i++)
        {
            idx[i] = el*flx->mesh.nvars + i; 
        }
        ierr = VecGetValues(flx->f,flx->mesh.nvars,idx,fm); CHKERRQ(ierr);
        ierr = VecGetValues(flx->w,flx->mesh.nvars,idx,wm); CHKERRQ(ierr);


        // Compute the residual
        for (int i = 0; i < flx->mesh.nvars; i++)
        {
            PetscScalar fph = 0.5*(fp[i] + fi[i] - flx->epsilon*flx->lp*(wp[i] - wi[i]));
            PetscScalar fmh = 0.5*(fi[i] + fm[i] - flx->epsilon*flx->lm*(wi[i] - wm[i]));

            res_elem[i] = (fph*sp - fmh*sm -q[i])/vi;
        }

        // Assemble the residual vector (indices are corresponding to iel-1 since boundaries are ignored)
        ierr = VecSetValues(res,flx->mesh.nvars,idx,res_elem,INSERT_VALUES); CHKERRQ(ierr);
        

    }
    ierr = VecAssemblyBegin(res); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(res); CHKERRQ(ierr);
    ierr = VecNorm(res,NORM_2,&resnrm); CHKERRQ(ierr);
    std::cout << "Log Non-Linear Residual norm is: " << std::log10(resnrm) << std::endl;
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
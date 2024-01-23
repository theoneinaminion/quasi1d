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

void Solver::adapt_time_step()
{
    dt = dt*prev_resnrm/resnrm;
    
}

PetscErrorCode Solver::setup_ksp()
{
    PetscErrorCode ierr;
    PC pc;
    ierr = KSPSetOperators(ksp,flx->A,flx->A); CHKERRQ(ierr);
    ierr = KSPSetType(ksp,KSPPREONLY); CHKERRQ(ierr);
    ierr = KSPGetPC(ksp,&pc); CHKERRQ(ierr);
    ierr = PCSetType(pc,PCLU); CHKERRQ(ierr);
    //ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);
    //ierr = KSPSetUp(ksp); CHKERRQ(ierr);
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
    ierr = compute_resnrm(); CHKERRQ(ierr);
    return ierr;
}

PetscErrorCode Solver::compute_resnrm()
{
    PetscErrorCode ierr;
    PetscInt idx[flx->mesh.nvars], el;
    PetscScalar relem[flx->mesh.nvars], nrm = 0.0; 

    for (PetscInt iel = 1; iel < flx->mesh.ngrid-1; iel++)
    {
        el = iel -1;
        for (int i = 0; i < flx->mesh.nvars; i++)
        {
            idx[i] = el*flx->mesh.nvars + i; 
        }
        ierr = VecGetValues(res,flx->mesh.nvars,idx,relem); CHKERRQ(ierr);
        nrm += relem[0]*relem[0]; 
        
    }
    resnrm = std::sqrt(nrm);
    return ierr;
    
}

PetscErrorCode Solver::update_solution()
{
    PetscErrorCode ierr;
    PetscInt idx[flx->mesh.nvars], el;
    PetscScalar welem[flx->mesh.nvars], dwelem[flx->mesh.nvars];

    for (PetscInt iel = 1; iel < flx->mesh.ngrid-1; iel++)
    {
        for (int i = 0; i < flx->mesh.nvars; i++)
        {
            idx[i] = iel*flx->mesh.nvars + i;
        }
        ierr = VecGetValues(flx->w,flx->mesh.nvars,idx,welem); CHKERRQ(ierr);

        el = iel-1;
        for (int i = 0; i < flx->mesh.nvars; i++)
        {
            idx[i] = el*flx->mesh.nvars + i;
        }

        ierr = VecGetValues(dw,flx->mesh.nvars,idx,dwelem); CHKERRQ(ierr);
        for (int i = 0; i < flx->mesh.nvars; i++)
        {
            welem[i] = welem[i] + relax*dwelem[i];
        }
        ierr = VecSetValues(flx->w,flx->mesh.nvars,idx,welem,INSERT_VALUES); CHKERRQ(ierr);
    }
    ierr = VecAssemblyBegin(flx->w); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(flx->w); CHKERRQ(ierr);

    // Update Boundaries 
    ierr = flx->outlet_bc(); CHKERRQ(ierr); // Already has a VecAssembly part
    ierr = flx->inlet_bc(); CHKERRQ(ierr); // Already has a VecAssembly part
    return ierr;

}

PetscErrorCode Solver::writePetscObj(const Vec &v, std::string name)
	
{

	PetscViewer viewer;
    PetscErrorCode ierr;
	const std::string namefin = name + ".dat";   
    std::cout << "Writing the file: " << namefin << std::endl;
	ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, namefin.c_str(), &viewer); CHKERRQ(ierr);
	ierr = VecView(v, viewer); CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
	return 0;

} 

PetscErrorCode Solver::write_soln()
{
    PetscErrorCode ierr;

    std::filesystem::path dir("data");
    if (!std::filesystem::exists(dir)) {
        std::filesystem::create_directory(dir);
}
    chdir("data");
    ierr = writePetscObj(flx->M,"mach"); CHKERRQ(ierr);
    ierr = writePetscObj(flx->mesh.xc,"cell_centers"); CHKERRQ(ierr);
    ierr = writePetscObj(flx->p,"pressure"); CHKERRQ(ierr);
    ierr = writePetscObj(flx->rho,"density"); CHKERRQ(ierr);
    ierr = writePetscObj(dw,"Soln"); CHKERRQ(ierr);
    ierr = writePetscObj(res,"residual"); CHKERRQ(ierr);
    return ierr;

}

PetscErrorCode Solver::writePetscObj(const Mat &A, std::string name)
	
{

    PetscViewer viewer;
    PetscErrorCode ierr;
    chdir("data");
    const std::string namefin = name + ".m";
    std::cout << "Writing the file: " << namefin << std::endl;
    ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, namefin.c_str(), &viewer); CHKERRQ(ierr);
    ierr = PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB); CHKERRQ(ierr);
    ierr = MatView(A, viewer); CHKERRQ(ierr);

    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
    return ierr;

} 


PetscErrorCode Solver::solve()
{
 
    PetscErrorCode ierr;
    ierr = 0;
    PC pc;
    PetscInt iter = 0;

    //Initializing 
    ierr = flx->initialize_primitives();CHKERRQ(ierr);
    ierr = flx->assemble_conservative_vec();CHKERRQ(ierr);

    // Setup KSP
    ierr = setup_ksp(); CHKERRQ(ierr);

    while ((resnrm > restol) && (iter <= maxiter))
    {
        prev_resnrm = resnrm;
        iter = iter + 1;
       
        ierr = flx->assemble_jacobian(dt);CHKERRQ(ierr);
        ierr = compute_residual();CHKERRQ(ierr);

        if ((iter%10 == 0) || iter < 10)
            std::cout << "Log Non-Linear Residual norm at iteration " << iter << " is: " << std::log10(resnrm) << std::endl;
        
        ierr = KSPSolve(ksp, res, dw); CHKERRQ(ierr);
        
        ierr = update_solution(); CHKERRQ(ierr);
        writePetscObj(flx->A,"A");
        adapt_time_step();
        ierr = flx->conserved_to_primitive(); CHKERRQ(ierr); 
    }
      
   return ierr;
    

}
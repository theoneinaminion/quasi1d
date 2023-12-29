#include "space.hpp"

void Mesh::generate_grid() {

    VecCreateSeq(PETSC_COMM_SELF,ngrid+1,&xw);
    VecCreateSeq(PETSC_COMM_SELF,ngrid,&xc);
    VecCreateSeq(PETSC_COMM_SELF,ngrid+1,&sw);
    VecCreateSeq(PETSC_COMM_SELF,ngrid,&sc);
    VecCreateSeq(PETSC_COMM_SELF,ngrid,&vol);
    
    dx = (x_max - x_min)/(ngrid);
    
    PetscScalar xw_arr;
    PetscScalar xc_arr;
    PetscScalar sw_arr;
    PetscScalar sc_arr;
    PetscScalar vol_arr;


    for (int i = 0; i <= ngrid; i++) {
        xw_arr = x_min + i*dx;
        sw_arr = 1-h*pow(sin(PETSC_PI*pow(xw_arr,t1)),t2);

        //VecSetValues(xw,ngrid+1,&i,xw_arr,INSERT_VALUES);
        VecSetValues(xw,1,&i,&xw_arr,INSERT_VALUES);
        VecSetValues(sw,1,&i,&sw_arr,INSERT_VALUES);
        if (i < ngrid) {
            xc_arr = x_min + (i+0.5)*dx;
            sc_arr = 1-h*pow(sin(PETSC_PI*pow(xc_arr,t1)),t2);
            vol_arr = sc_arr*dx;

            VecSetValues(xc,1,&i,&xc_arr,INSERT_VALUES);
            VecSetValues(sc,1,&i,&sc_arr,INSERT_VALUES);
            VecSetValues(vol,1,&i,&vol_arr,INSERT_VALUES);
        }
    }
    
    
    
    
    //VecSetValues(sc,ngrid,NULL,sc_arr,INSERT_VALUES);
    

    VecAssemblyBegin(xw);
    VecAssemblyEnd(xw);
    VecAssemblyBegin(xc);
    VecAssemblyEnd(xc);
    VecAssemblyBegin(sw);
    VecAssemblyEnd(sw);
    VecAssemblyBegin(sc);
    VecAssemblyEnd(sc);
    VecAssemblyBegin(vol);
    VecAssemblyEnd(vol);

}

flow::flow(const Mesh &msh, const PetscScalar p_ratio) {

    /**
     * @brief Constructor for the fluid class.
     * 
     */

    

    pr = p_ratio;
    p_exit = pr*p_t;
    rho_t = p_t/(R*T_t); //Total density
    c_t = sqrt(gamma*R*T_t); //Total speed of sound
    E_tot = p_t/(gamma-1) + init_mach*rho_t*pow(0.5*c_t,2); //Total energy when velocity = 0.5*ct

    p_sc = p_t;
    T_sc = T_t;
    rho_sc = rho_t;

    //Performing the scaling operation
    p_t           = p_t/p_sc;
    T_t           = T_t/T_sc;
    rho_t         = rho_t/rho_sc; 

    R_sc         = p_sc/(rho_sc*T_sc);
    R            = R/R_sc;      
    cp           = gamma*R/(gamma-1);    /**< cp corresponding to scaled R*/

    u_sc        = sqrt(p_sc/rho_sc);
    c_t         = c_t/u_sc;  

    E_sc        = rho_sc*pow(u_sc,2);
    E_tot       = E_tot/E_sc;

    //Construct the vectors of flow properties at cell centres
   
    VecDuplicate(msh.xc,&p); 
    VecDuplicate(msh.xc,&T); 
    VecDuplicate(msh.xc,&rho); 
    VecDuplicate(msh.xc,&u); 
    VecDuplicate(msh.xc,&c); 
    VecDuplicate(msh.xc,&E); 
    VecDuplicate(msh.xc,&M); 

    
     
}


# if 1
flow::~flow() {

    /**
     * @brief Destroy the flow object
     * 
     */
    int ierr;
    std::cout << "Destroying flow object\n";
    ierr = VecDestroy(&p);
    if(ierr)
		std::cout << "Couldn't destroy p!\n";
    ierr = VecDestroy(&T);
    if(ierr)
		std::cout << "Couldn't destroy T!\n";
    ierr = VecDestroy(&rho);
    if(ierr)
		std::cout << "Couldn't destroy rho!\n";
    ierr = VecDestroy(&u);
    if(ierr)
		std::cout << "Couldn't destroy u!\n";
    ierr = VecDestroy(&c);
    if(ierr)
		std::cout << "Couldn't destroy c!\n";
    ierr = VecDestroy(&E);
    if(ierr)
		std::cout << "Couldn't destroy E!\n";
    ierr = VecDestroy(&M);
    if(ierr)
		std::cout << "Couldn't destroy M!\n";
} 
#endif



flux::flux(const Mesh &msh, const PetscScalar p_ratio):flow(msh,p_ratio) {

    /**
     * @brief Construct a new flux object
     * 
     */

    mesh = msh; //Assigning the mesh object in the flux class
    
    VecCreateSeq(PETSC_COMM_SELF,msh.ngrid*msh.nvars,&f);
    VecCreateSeq(PETSC_COMM_SELF,msh.ngrid*msh.nvars,&w);
    VecCreateSeq(PETSC_COMM_SELF,msh.nvars,&w_elem);
    VecCreateSeq(PETSC_COMM_SELF,msh.nvars,&f_elem);
}


#if 1
flux::~flux() {

    /**
     * @brief Destroy the flux object
     * 
     */
    std::cout << "Destroying flux object\n";
    int ierr;
    ierr = VecDestroy(&f); 
    if(ierr)
		std::cout << "Couldn't destroy f!\n";
    ierr = VecDestroy(&w); 
    if(ierr)
		std::cout << "Couldn't destroy w!\n";

    ierr = VecDestroy(&w_elem); 
    if(ierr)
		std::cout << "Couldn't destroy welem!\n";
    ierr = VecDestroy(&f_elem); 
    if(ierr)
		std::cout << "Couldn't destroy felem!\n";
} 
#endif
    

PetscErrorCode flux::initialize_primitives() {
    /**
     * @brief Initialize the primitive variables
     * 
     */


    PetscScalar temp;
    PetscErrorCode ierr; 
    ierr = VecSet(p,p_t); CHKERRQ(ierr);
    ierr = VecSet(T,T_t); CHKERRQ(ierr);
    ierr = VecSet(rho,rho_t); CHKERRQ(ierr);
    ierr = VecSet(E,E_tot); CHKERRQ(ierr);
    temp = init_mach*c_t; CHKERRQ(ierr);
    ierr = VecSet(u, temp); CHKERRQ(ierr);
    ierr = VecSet(c,c_t); CHKERRQ(ierr);
    ierr = VecSet(M,init_mach); CHKERRQ(ierr);

    return ierr;

}

PetscErrorCode flux::element_conservative(const PetscInt &elem){
    
    /**
    * @brief Generate conservative variables vector at each element.
    * 
    */
    PetscErrorCode ierr;
    PetscInt idx[mesh.nvars] = {0,1,2};
    PetscScalar temp[mesh.nvars], rho_elem, u_elem, E_elem;

    ierr = VecSet(w_elem,0.0); CHKERRQ(ierr);
    ierr = VecSet(f_elem,0.0); CHKERRQ(ierr);

    ierr = VecGetValues(rho,1, &elem,&rho_elem); CHKERRQ(ierr);
    ierr = VecGetValues(u,1, &elem,&u_elem); CHKERRQ(ierr);
    ierr = VecGetValues(E,1, &elem,&E_elem); CHKERRQ(ierr);

    temp[0] = rho_elem; 
    temp[1] = rho_elem*u_elem;
    temp[2] = E_elem;

    ierr = VecSetValues(w_elem,mesh.nvars,idx,temp,INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(w_elem); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(w_elem); CHKERRQ(ierr);

    return ierr;
    
}


PetscErrorCode flux::assemble_conservative_vec(){

    /**
     * @brief Assemble the global conservative vars vector
     * 
     */

    PetscErrorCode ierr;  
    PetscScalar temp[mesh.nvars];
     

    ierr = VecSet(w,0.0); CHKERRQ(ierr);

    for (int iel = 0; iel < mesh.ngrid; iel++) {

        ierr = element_conservative(iel);

        PetscInt idx[mesh.nvars] = {0,1,2};
        ierr = VecGetValues(w_elem,mesh.nvars,idx,temp); CHKERRQ(ierr);

        for (int j = 0; j < mesh.nvars; j++)
        {
            idx[j] = iel*mesh.nvars + j;
        }
        ierr = VecSetValues(w,mesh.nvars,idx,temp,INSERT_VALUES); CHKERRQ(ierr);
        
        
    }
    ierr = VecAssemblyBegin(w); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(w); CHKERRQ(ierr);

    return ierr;

}






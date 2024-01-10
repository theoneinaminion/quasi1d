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
    felem = new PetscScalar[mesh.nvars];
    welem = new PetscScalar[mesh.nvars];
    qelem = new PetscScalar[mesh.nvars];
    wb = new PetscScalar[mesh.nvars];

    Jel.resize(mesh.nvars,mesh.nvars);
    Qel.resize(mesh.nvars,mesh.nvars); 
    Jb.resize(mesh.nvars,mesh.nvars); 

    PetscInt ierr;
    std::cout << "Constructing System Vectors"<<std::endl;

    ierr = VecCreateSeq(PETSC_COMM_SELF,mesh.ngrid*mesh.nvars,&f);
    if(ierr)
		  std::cout << "Couldn't create f!\n";
    ierr = VecCreateSeq(PETSC_COMM_SELF,mesh.ngrid*mesh.nvars,&w);
    if(ierr)
		  std::cout << "Couldn't create w!\n";
    ierr = VecCreateSeq(PETSC_COMM_SELF,mesh.ngrid*mesh.nvars,&q);
    if(ierr)
		  std::cout << "Couldn't create w!\n";

    std::cout << "Constructing System Matrix"<<std::endl;

    ierr = MatCreate(PETSC_COMM_SELF,&A);
    if(ierr)
		  std::cout << "Couldn't create implicit matrix A!\n";
    ierr = MatSetType(A,MATSEQBAIJ);
    if(ierr)
		  std::cout << "Couldn't set type of A!\n";
    ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,mesh.nvars,mesh.nvars);
    if(ierr)
		  std::cout << "Couldn't set sizes of A!\n";
    ierr = MatSetBlockSize(A,mesh.nvars);
    if(ierr)
		  std::cout << "Couldn't set block size of A!\n";
    ierr = MatSetUp(A);
    if(ierr)
		  std::cout << "Couldn't set up internal data structures of A!\n";
    
    
}



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
    ierr = VecDestroy(&q); 
    if(ierr)
		  std::cout << "Couldn't destroy q!\n";
    ierr = MatDestroy(&A);
    if(ierr)
      std::cout << "Couldn't destroy A!\n";

    
    delete [] felem;
    delete [] welem;
    delete [] qelem;
    delete [] wb;
} 

    

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
    PetscScalar rho_elem, u_elem, E_elem;


    ierr = VecGetValues(rho,1, &elem,&rho_elem); CHKERRQ(ierr);
    ierr = VecGetValues(u,1, &elem,&u_elem); CHKERRQ(ierr);
    ierr = VecGetValues(E,1, &elem,&E_elem); CHKERRQ(ierr);

    welem[0] = rho_elem; 
    welem[1] = rho_elem*u_elem;
    welem[2] = E_elem;

    return ierr;
    
}


PetscErrorCode flux::assemble_conservative_vec(){

    /**
     * @brief Assemble the global conservative vars vector
     * 
     */

    PetscErrorCode ierr;  

    ierr = VecSet(w,0.0); CHKERRQ(ierr);

    for (int iel = 0; iel < mesh.ngrid; iel++) {

        ierr = element_conservative(iel);

        PetscInt idx[mesh.nvars];
        
        for (int j = 0; j < mesh.nvars; j++)
        {
            idx[j] = iel*mesh.nvars + j;
        }
        ierr = VecSetValues(w,mesh.nvars,idx,welem,INSERT_VALUES); CHKERRQ(ierr);
        
        
    }
    ierr = VecAssemblyBegin(w); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(w); CHKERRQ(ierr);

    return ierr;

}

PetscErrorCode flux::element_flux(const PetscInt &elem){

  /**
   * @brief Generate the flux vector for each element.
   *  
   */
    PetscErrorCode ierr;
    PetscInt idx[mesh.nvars];

    for (int i = 0; i < mesh.nvars; i++) {
        idx[i] =elem*mesh.nvars + i;
    }
    ierr = VecGetValues(w, mesh.nvars, idx, welem); CHKERRQ(ierr);

    //Components of flux vector
    felem[0] = welem[0];
    felem[1] = welem[1]/welem[0];
    felem[2] = (gamma-1)*(welem[2] - 0.5*welem[0]*pow(felem[1],2));
    return ierr;


}

PetscErrorCode flux::assemble_flux_vec(){

  /**
     * @brief Assemble the global flux vector
     * 
     */

    PetscErrorCode ierr;  
    PetscInt idx[mesh.nvars];
    
    ierr = VecSet(f,0.0); CHKERRQ(ierr);

    for (int iel = 0; iel < mesh.ngrid; iel++) {

      ierr = element_flux(iel);

     
      for (int j = 0; j < mesh.nvars; j++)
      {
          idx[j] = iel*mesh.nvars + j;
      }
        ierr = VecSetValues(f,mesh.nvars,idx,felem,INSERT_VALUES); CHKERRQ(ierr);
        
        
    }
    ierr = VecAssemblyBegin(f); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(f); CHKERRQ(ierr);

    return ierr;

}

PetscErrorCode flux::conserved_to_primitive(){
  /**
   * @brief Convert the conserved variables to primitive variables
   * 
   */

  PetscErrorCode ierr;
  PetscInt idx[mesh.nvars];
  PetscScalar relem,pelem,uelem,Eelem,celem, var;

  for (int iel = 0; iel < mesh.ngrid; iel++)
  {
    for (int i = 0; i < mesh.nvars; i++)
    {
      idx[i] = iel*mesh.nvars + i;
    }
    ierr = VecGetValues(w,mesh.nvars,idx,welem); CHKERRQ(ierr);

    relem = welem[0];
    ierr = VecSetValues(rho,1,&iel,&relem, INSERT_VALUES); CHKERRQ(ierr);
    uelem = welem[1]/welem[0];
    ierr = VecSetValues(u,1,&iel,&uelem, INSERT_VALUES); CHKERRQ(ierr);
    Eelem = welem[2];
    ierr = VecSetValues(E,1,&iel,&Eelem, INSERT_VALUES); CHKERRQ(ierr);
    pelem = (gamma-1)*(Eelem - 0.5*pow(uelem,2)*relem);
    ierr = VecSetValues(p,1,&iel,&pelem, INSERT_VALUES); CHKERRQ(ierr);
    celem = sqrt(gamma*pelem/relem);
    ierr = VecSetValues(c,1,&iel,&celem, INSERT_VALUES); CHKERRQ(ierr);
    var = pelem/(R*relem);
    ierr = VecSetValues(T,1,&iel,&var, INSERT_VALUES); CHKERRQ(ierr);
    var = uelem/celem;
    ierr = VecSetValues(M,1,&iel,&var, INSERT_VALUES); CHKERRQ(ierr);
  }
  ierr = VecAssemblyBegin(p); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(p); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(T); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(T); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(rho); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(rho); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(u); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(u); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(c); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(c); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(E); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(E); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(M); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(M); CHKERRQ(ierr);


  return ierr;
}

PetscErrorCode flux:: element_source(const PetscInt &elem){

  /**
   * @brief Generate the source vector for each element.
   *  
   */
    // We calculate the source term only for the interior cells. 
    assert (elem>0);
    assert (elem<mesh.ngrid);
    PetscErrorCode ierr;
    PetscScalar pelem, si, sm;
 

    ierr = VecGetValues(p, 1, &elem, &pelem); CHKERRQ(ierr);
    ierr = VecGetValues(mesh.sc, 1, &elem, &si); CHKERRQ(ierr);
    PetscInt el = elem -1;
    ierr = VecGetValues(mesh.sc, 1, &el, &sm); CHKERRQ(ierr);
    //Components of source vector
    qelem[0] = 0;
    qelem[1] = pelem*(si-sm);
    qelem[2] = 0;
    
    return ierr;

}

PetscErrorCode flux::assemble_source_vec(){

  PetscErrorCode ierr;
  ierr = VecSet(q,0.0); CHKERRQ(ierr);

  for (int iel = 1; iel < mesh.ngrid-1; iel++) {

    PetscInt idx[mesh.nvars];
    ierr = element_source(iel);

    for (int j = 0; j < mesh.nvars; j++)
    {
        idx[j] = iel*mesh.nvars + j;
    }
    ierr = VecSetValues(q,mesh.nvars,idx,qelem,INSERT_VALUES); CHKERRQ(ierr);
    
    
  }
  ierr = VecAssemblyBegin(q); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(q); CHKERRQ(ierr);

  return ierr;

}

PetscErrorCode flux::int_element_flux_jacobian(const PetscInt &elem)
{
 /**
  * @brief Analytical jacobian for each interior element
  * 
  */
  PetscErrorCode ierr;
  PetscScalar wel[mesh.nvars];
  PetscInt idx[mesh.nvars];

  // We calculate the jacobian only for the interior cells. Direct BCs are used for boundary cells. 
  assert (elem>0);
  assert (elem<mesh.ngrid);

  Jel = Eigen::MatrixXd::Zero(mesh.nvars,mesh.nvars);
  for (int i = 0; i < mesh.nvars; i++)
  {
    idx[i] = elem*mesh.nvars + i;
  }
  ierr = VecGetValues(w,mesh.nvars,idx,wel); CHKERRQ(ierr);
  
  // The analytical Jacobian

  //Row 1
  Jel(0,0) = 0;Jel(0,1) = 1;Jel(0,2) = 0; 

  //Row 2
  Jel(1,0) = -0.5*(3-gamma)*pow(wel[1],2)/pow(wel[0],2);
  Jel(1,1) = (3-gamma)*wel[1]/wel[0];
  Jel(1,2) = gamma-1;

  //Row 3
  Jel(2,0) = (gamma-1)*pow(wel[1],3)/pow(wel[0],3) - gamma*wel[2]*wel[1]/pow(wel[0],2); 
  Jel(2,1) = gamma*wel[2]/wel[0] - 1.5*(gamma-1)*pow(wel[1],2)/pow(wel[0],2);
  Jel(2,2) = gamma*wel[1]/wel[0];

  return ierr;
}

 PetscErrorCode flux::int_element_source_jacobian(const PetscInt &elem){

  /**
   * @brief Analytical source vector jacobian for each element
   * 
   */
  PetscErrorCode ierr;
  PetscScalar wel[mesh.nvars], si, sm, k;
  PetscInt idx[mesh.nvars], el;

  // We calculate the jacobian only for the interior cells. Direct BCs are used for boundary cells. 
  assert (elem>0);
  assert (elem<mesh.ngrid);

  Qel = Eigen::MatrixXd::Zero(mesh.nvars,mesh.nvars);
  for (int i = 0; i < mesh.nvars; i++)
  {
    idx[i] = elem*mesh.nvars + i;
  }
  ierr = VecGetValues(w,mesh.nvars,idx,wel); CHKERRQ(ierr);
  ierr = VecGetValues(mesh.sc, 1, &elem, &si); CHKERRQ(ierr);
  el = elem -1;
  ierr = VecGetValues(mesh.sc, 1, &el, &sm); CHKERRQ(ierr);

  //Row 2 (Row 1 & 3 are zeros)
  k = si-sm;
  Qel(1,0) = k*(gamma-1)*0.5*pow(wel[1],2)/pow(wel[0],2);
  Qel(2,0) = -k*(gamma-1)*wel[1]/wel[0];
  Qel(3,0) = k*(gamma-1);
 
  return ierr;

 }

PetscErrorCode flux::inlet_bc(){
  /**
     * @brief Inlet Boundary conditions update using Riemann Invariants. Overloaded function
     * 
  */

  PetscErrorCode ierr;
  PetscScalar w1[mesh.nvars], w2[mesh.nvars];
  PetscInt idx[mesh.nvars];

  // Primitives at elem = 0 and 1.   
  PetscScalar r1,u1,p1,T1,c1;
  PetscScalar r2,u2,p2,c2;

  PetscInt elem = 1;
  for (int i = 0; i < mesh.nvars; i++)
  {
    idx[i] = elem*mesh.nvars + i;
  }
  ierr = VecGetValues(w,mesh.nvars,idx,w2); CHKERRQ(ierr);
  
  elem = elem-1;
  for (int i = 0; i < mesh.nvars; i++)
  {
    idx[i] = elem*mesh.nvars + i;
  }
  ierr = VecGetValues(w,mesh.nvars,idx,w1); CHKERRQ(ierr);
  

  //Getting values at first element. (No using primitive vectors directly coz they are not updated yet)
  r1 = w1[0];
  u1 = w1[1]/w1[0];
  p1 = (gamma-1)*(w1[2] - 0.5*pow(u1,2)*r1);
  T1 = p1/(R*r1);
  c1 = sqrt(gamma*p1/r1);

  //Getting values at second element.
  r2 = w2[0];
  u2 = w2[1]/w2[0];
  p2 = (gamma-1)*(w2[2] - 0.5*pow(u2,2)*r2);
  //T2 = p2/(R*r2);
  c2 = sqrt(gamma*p2/r2);

  PetscScalar mi = (u1+u2)/(c1+c2);

  // Updating the first element only duriing subsonic flow
  if(mi<1) {

    PetscScalar Rm  = -u2-2*c2/(gamma-1);
    PetscScalar c02 = pow(c2,2) + 0.5*(gamma-1)*pow(u2,2);
    PetscScalar l   = -Rm*(gamma-1)/(gamma+1);
    PetscScalar k   = (c02/pow(Rm,2))*(gamma-1)/(gamma+1) - 0.5*(gamma-1);
    c1  = l*(1+sqrt(k));

    T1 = T_t*(pow(c1,2)/c02);
    p1 = p_t*pow((T1/T_t),(gamma/(gamma-1)));
    r1 = p1/(R*T1);
    u1 = pow((2*cp*(T_t - T1)),0.5);
    
    PetscScalar E1 = p1/(gamma-1) + 0.5*r1*pow(u1,2);

    w1[0] = r1;
    w1[1] = r1*u1;
    w1[2] = E1;
  }
  
  //Updating the first element
  ierr = VecSetValues(w,mesh.nvars,idx,w1,INSERT_VALUES); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(w); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(w); CHKERRQ(ierr);
  return ierr;
}

 

 PetscErrorCode flux::inlet_bc(PetscInt i){

  /**
     * @brief Inlet Boundary conditions function for updating jacobian.
     * 
     * 
  */

  PetscErrorCode ierr;
  PetscScalar w1[mesh.nvars];

  // Primitives at elem = 0 and 1.   
  PetscScalar r1,u1,p1,T1,c1;
  PetscScalar r2,u2,p2,c2;
  PetscInt idx[mesh.nvars]={0,1,2};

  ierr = VecGetValues(w,mesh.nvars,idx,w1); CHKERRQ(ierr);

  //Getting values at first element. (No using primitive vectors directly coz they are not updated yet)
  r1 = w1[0];
  u1 = w1[1]/w1[0];
  p1 = (gamma-1)*(w1[2] - 0.5*pow(u1,2)*r1);
  T1 = p1/(R*r1);
  c1 = sqrt(gamma*p1/r1);

  //Getting values at second element.
  r2 = wb[0];
  u2 = wb[1]/wb[0];
  p2 = (gamma-1)*(wb[2] - 0.5*pow(u2,2)*r2);
  //T2 = p2/(R*r2);
  c2 = sqrt(gamma*p2/r2);

  PetscScalar mi = (u1+u2)/(c1+c2);

  // Updating the first element only duriing subsonic flow
  if(mi<1) {

    PetscScalar Rm  = -u2-2*c2/(gamma-1);
    PetscScalar c02 = pow(c2,2) + 0.5*(gamma-1)*pow(u2,2);
    PetscScalar l   = -Rm*(gamma-1)/(gamma+1);
    PetscScalar k   = (c02/pow(Rm,2))*(gamma-1)/(gamma+1) - 0.5*(gamma-1);
    c1  = l*(1+sqrt(k));

    T1 = T_t*(pow(c1,2)/c02);
    p1 = p_t*pow((T1/T_t),(gamma/(gamma-1)));
    r1 = p1/(R*T1);
    u1 = pow((2*cp*(T_t - T1)),0.5);
    
    PetscScalar E1 = p1/(gamma-1) + 0.5*r1*pow(u1,2);

    // This is temprary update of w1 stored here. It will be used in jacobian updating process. 
    welem[0] = r1;
    welem[1] = r1*u1;
    welem[2] = E1;
  } 

  return ierr;

 }

 
  PetscErrorCode flux::outlet_bc(){
    /**
       * @brief Outlet Boundary conditions update using Riemann Invariants. Overloaded function
       * 
    */


    PetscErrorCode ierr;
    PetscScalar w1[mesh.nvars], w2[mesh.nvars]; // Conservatives at elem = ngrid-1 and ngrid respectively.
    PetscInt idx[mesh.nvars];

    // Primitives at elem = ngrid-1 and ngrid respectively.    
    PetscScalar r1,u1,p1,c1;
    PetscScalar r2,u2,p2,c2;
    PetscInt elem = mesh.ngrid-1;
    for (int i = 0; i < mesh.nvars; i++)
    {
      idx[i] = elem*mesh.nvars + i;
    }

    ierr = VecGetValues(w,mesh.nvars,idx,w1); CHKERRQ(ierr);
    elem = mesh.ngrid;
    for (int i = 0; i < mesh.nvars; i++)
    {
      idx[i] = elem*mesh.nvars + i;
    }
    ierr = VecGetValues(w,mesh.nvars,idx,w2); CHKERRQ(ierr);

    //Getting values at ngrid-1 element. (No using primitive vectors directly coz they are not updated yet)
    r1 = w1[0];
    u1 = w1[1]/w1[0];
    p1 = (gamma-1)*(w1[2] - 0.5*pow(u1,2)*r1);
    //T1 = p1/(R*r1);
    c1 = sqrt(gamma*p1/r1);

    //Getting values at ngrid element.
    r2 = w2[0];
    u2 = w2[1]/w2[0];
    p2 = (gamma-1)*(w2[2] - 0.5*pow(u2,2)*r2);
    //T2 = p2/(R*r2);
    c2 = sqrt(gamma*p2/r2);

    PetscScalar me = (u1+u2)/(c1+c2);

    if(me<1)
    {
      PetscScalar c0 = sqrt(pow(c1,2) + 0.5*(gamma-1)*pow(u1,2));  // Total speed of sound
      p2 = p_exit;
      r2 = r1 + (p2 - p1)/(pow(c0,2));
      u2 = u1 + (p1 - p2)/(r1*c0);
      PetscScalar E2 = p2/(gamma-1) + 0.5*r2*pow(u2,2);

      w2[0] = r2;
      w2[1] = r2*u2;
      w2[2] = E2;
    }

    //Updating the last element
    ierr = VecSetValues(w,mesh.nvars,idx,w2,INSERT_VALUES); CHKERRQ(ierr); //idx is already at ngrid
    ierr = VecAssemblyBegin(w); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(w); CHKERRQ(ierr);
    return ierr;
  }

  PetscErrorCode flux::outlet_bc(PetscInt i){
    /**
       * @brief Outlet Boundary conditions update using Riemann Invariants. Overloaded function
       * 
    */


    PetscErrorCode ierr;
    PetscScalar w2[mesh.nvars]; // Conservatives at elem = ngrid-1 and ngrid respectively.
    PetscInt idx[mesh.nvars];

    // Primitives at elem = ngrid-1 and ngrid respectively.    
    PetscScalar r1,u1,p1,c1;
    PetscScalar r2,u2,p2,c2;
    PetscInt elem = mesh.ngrid;
    
    for (int i = 0; i < mesh.nvars; i++)
    {
      idx[i] = elem*mesh.nvars + i;
    }
    ierr = VecGetValues(w,mesh.nvars,idx,w2); CHKERRQ(ierr);

    //Getting values at ngrid-1 element. (No using primitive vectors directly coz they are not updated yet)
    r1 = wb[0];
    u1 = wb[1]/wb[0];
    p1 = (gamma-1)*(wb[2] - 0.5*pow(u1,2)*r1);
    //T1 = p1/(R*r1);
    c1 = sqrt(gamma*p1/r1);

    //Getting values at ngrid element.
    r2 = w2[0];
    u2 = w2[1]/w2[0];
    p2 = (gamma-1)*(w2[2] - 0.5*pow(u2,2)*r2);
    //T2 = p2/(R*r2);
    c2 = sqrt(gamma*p2/r2);

    PetscScalar me = (u1+u2)/(c1+c2);

    if(me<1)
    {
      PetscScalar c0 = sqrt(pow(c1,2) + 0.5*(gamma-1)*pow(u1,2));  // Total speed of sound
      p2 = p_exit;
      r2 = r1 + (p2 - p1)/(pow(c0,2));
      u2 = u1 + (p1 - p2)/(r1*c0);
      PetscScalar E2 = p2/(gamma-1) + 0.5*r2*pow(u2,2);

      // Temporarily storing it in welem to be used in Jacobian update
      welem[0] = r2;
      welem[1] = r2*u2;
      welem[2] = E2;
    }


    return ierr;
  }


 





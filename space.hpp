#pragma once

#include <iostream>
#include <vector>
#include <cmath>
#include "petscksp.h"
#include "petscmath.h"


class Mesh {

 public:    

 /**
  * @brief Generate the grid at the constructor level.
  * 
  */
    

    PetscScalar x_max = 1.; /**< domain begins*/
    PetscScalar x_min = 0.; /**< domain ends*/
    PetscInt    ngrid = 50; /**< total grid points*/
    PetscInt    nvars = 3; /**< total conservatives variables for Quasi 1D*/
    PetscScalar dx; /**< grid spacing*/


    /**< Nozzle Geometry*/
    PetscScalar h = 0.15;
    PetscScalar t1 = 0.8;
    PetscScalar t2 = 3;

    Vec xw; /**< wall points*/
    Vec xc; /**< cell centers*/
    Vec sw; /**< Area at wall points*/
    Vec sc; /**<Area at  cell centers*/
    Vec vol; /**<Volume*/

  /**
   * @brief Generate the grid and initialize the vectors.
   * 
   */
    void generate_grid();


};

class flow {

  /**
   * @brief Construct a new flow object and initialize the primitives after scaling.s
   * 
   * @param msh Mesh object
   * @param pr Pressure ratio
   */


  public:

    flow(const Mesh &msh, const PetscScalar p_ratio); // Constructor
    virtual ~flow(); // Destructor
    Vec p; /**< Pressure*/
    Vec T; /**< Temperature*/
    Vec rho; /**< Density*/
    Vec u; /**< Velocity*/
    Vec c; /**< speed of sound*/
    Vec E; /**< energy*/
    Vec M; /**< Mach number*/
    

    // Exit conditions
    PetscScalar p_exit; /**< Exit pressure*/ 
    PetscScalar T_exit; /**< Exit temperature*/
    PetscScalar rho_exit; /**< Exit density*/
    PetscScalar u_exit;   /**< Exit velocity*/
    PetscScalar c_exit;   /**< Exit speed of sound*/ 
    PetscScalar E_exit;  /**< Exit energy*/
    PetscScalar M_exit; /**< Exit Mach number*/

    //Constants
    PetscScalar gamma    = 1.5;                /**< Specific heat ratio*/
    PetscScalar p_t      = 2117;               /**< Total pressure*/
    PetscScalar T_t      = 531.2;              /**< Total temperature*/
    PetscScalar c_t;                            /**< Total speed of sound*/
    PetscScalar E_tot;                           /**< Energy correspondng to total conditions*/         
    PetscScalar pr;                           /**< Pressure ratio*/
    PetscScalar rho_t;                        /**< Total Density*/
    PetscScalar R = 1716;                     /**< Gas constant*/
    PetscScalar cp;                           /**< Total pressure*/
    PetscScalar init_mach = 0.1;              /**< Initial Mach number just to initialize things*/


   
    //Scales 
    PetscScalar p_sc;                     /**< Pressure scale*/
    PetscScalar T_sc;                     /**< Temperature scale*/
    PetscScalar rho_sc;                   /**< Density scale*/
    PetscScalar u_sc;                     /**< Velocity scale*/
    PetscScalar E_sc;                     /**< Energy scale*/
    PetscScalar L_sc = 1.;                /**< Length scale*/
    PetscScalar R_sc;                     /**< Gas constant scale*/ 


};

class flux:public flow {
  /**
   * @brief Construct a new flux object
   * 
   * @param msh Mesh object
   * @param pr Pressure ratio
   */

  public:
    flux(const Mesh &msh, const PetscScalar p_ratio);
    virtual ~flux();
    Vec f; /**< Flux vector*/
    Vec w; /**< Conservative variables vector*/
    Mesh mesh; /**< Mesh object*/

    Vec f_elem; /**< Flux vector at each element*/
    Vec w_elem; /**< Conservative variables vector at each element*/

    PetscScalar epsilon = 0.08; /**< Scalar Dissipation*/

    /**
     * @brief Initialize the primitive variables. 
     * 
     */
    PetscErrorCode initialize_primitives();

    /**
     * @brief Generate the flux vector for each element.
     * 
     * @param elem  Element number
     */
    PetscErrorCode element_flux(const PetscInt &elem);

    /**
     * @brief Generate the conservative vars vector for each element.
     * 
     * @param elem  Element number
     */
    PetscErrorCode element_conservative(const PetscInt &elem);

    /**
     * @brief Assemble flux vector in the entire domain
     * 
     */
    PetscErrorCode assemble_flux_vec();

    /**
     * @brief Conservative variables vector for the entire domain.
     * 
     */
    PetscErrorCode assemble_conservative_vec();

    /**
     * @brief Conserved to Primitive variables.
     * 
     */

    PetscErrorCode conserved_to_primitive();

   


};

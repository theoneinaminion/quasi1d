#pragma once
#include "space.hpp"
#include "petscmat.h" 
#include <filesystem>
#include <unistd.h>

#define DEBUG 0

class Solver{
    /**
     * @brief Solver class.
     * 
     */

    public: 
    Solver(flux* flux);
    virtual ~Solver();

    flux* flx; /**< Flux object*/
    
    Vec res; /**< Residual vector*/
    Vec dw; /**< Increment conservative vector*/
    KSP ksp; /**< KSP Object*/

    PetscScalar resnrm = 10; /**< Non-Linear Residual norm initialized to 10*/
    PetscScalar prev_resnrm; /**< Non-Linear Residual norm at previous time step*/

    PetscScalar dt; /**< Time step*/
    PetscScalar CFL; /**< CFL number*/
    PetscScalar relax = 0.01; /**< Under Relaxation factor*/
    PetscScalar restol = 1e-6; /**< Tolerance for Non-Linear solver*/
    PetscInt    maxiter ; /**< Maximum number of iterations for Non-Linear solver*/

    /**
     * @brief Setup KSP object
     * 
     */
    PetscErrorCode setup_ksp();
    
    /**
     * @brief New time step after for next iteration 
     * dt = prev_resnrm/resnrm
     * 
     */
    void adapt_time_step();

 
    /**
     * @brief Assemble the global non-lin residual vector
     * 
     */
    PetscErrorCode compute_residual();

    /**
     * @brief Compute the non-linear residual norm of mass eqn
     * 
     */
    PetscErrorCode compute_resnrm();

    /**
     * @brief Solve the non-linear system
     */
    PetscErrorCode solve();

    /**
     * @brief Update the solution vector flx.w = flx.w + dw (including updating boundaries)
     * 
     */

    PetscErrorCode update_solution();

    /**
     * @brief Write a Petsc Object as a dat file
     * @param v Petsc Object
     * @param name Name of the file excluding extension
     */
    PetscErrorCode writePetscObj(const Vec &v, std::string name);

    /**
     * @brief Write a Petsc Object as a dat file (to be used while denugging)
     * @param A Petsc Object
     * @param name Name of the file excluding extension
     */
    PetscErrorCode writePetscObj(const Mat &A, std::string name);


    /**
     * @brief Writes Mach, xw, p as dat files. 
     * 
     * @return PetscErrorCode 
     */
    PetscErrorCode write_soln();

};
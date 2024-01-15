#include "space.hpp"

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

    PetscScalar resnrm; /**< Non-Linear Residual norm*/
    PetscScalar prev_resnrm; /**< Non-Linear Residual norm at previous time step*/

    PetscScalar dt; /**< Time step*/
    PetscScalar CFL; /**< CFL number*/

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
    void new_time_step();

 
    /**
     * @brief Assemble the global non-lin residual vector
     * 
     */
    PetscErrorCode compute_residual();

    /**
     * @brief Solve the non-linear system
     */
    PetscErrorCode solve();

    /**
     * @brief Update the solution vector flx.w = flx.w + dw
     * 
     */

    PetscErrorCode update_solution();

};
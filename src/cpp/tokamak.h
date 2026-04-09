/**
 * @file tokamak.h
 * 
 * @brief Functions for solving the Grad-Shafranov equation
 * 
 * @date 2026/4/8
 */

#ifndef TOKAMAK_H
#define TOKAMAK_H

#include "mole.h"

/**
 * @brief Path to file containing x-coordinates of tokamak boundary
 */
#define XLIM_PATH ""

/**
 * @brief Path to file containing y-coordinates of tokamak boundary
 */
#define YLIM_PATH ""

/**
 * @brief Path to file containing coil information
 */
#define COIL_PATH ""

/**
 * @brief Constructs the laplacian matrix for u_xx / (alpha^2) + u_yy / (beta^2)
 * 
 * @param k Order of accuracy
 * @param m Number of cells in x-direction
 * @param dx Step size in x-direction
 * @param n Number of cells in y-direction
 * @param dy Step size in y-direction
 * @param alpha parameter to scale x derivative
 * @param beta parameter to scale y derivative
 */
sp_mat weighted_laplacian(u16 k, u32 m, Real dx, u32 n, Real dy, Real alpha, Real beta);

/**
 * @brief Constructs the x or y coordinates of the plasma domain using a reference grid
 * 
 * @param k Order of accuracy
 * @param left coordinates of the left boundary, ordered bottom to top
 * @param right coordinates of the right boundary, ordered bottom to top
 * @param bottom coordinates of the bottom boundary, ordered left to right
 * @param top coordinates of the top boundary, ordered left to right
 */
vec plasma_grid(u16 k, vec left, vec right, vec bottom, vec top);

/**
 * @brief Constructs to x- or y- coordinates of the vacuum domain using a reference grid
 * 
 * @param k Order of accuracy
 * @param bottom coordinates of the plasma boundary, ordered left to right
 * @param top coordinates of the tokamak boundary, ordered left to right
 */
vec vacuum_grid(u16 k, u32 n, vec bottom, vec top);

/**
 * 
 */
mat read_coils(const char* coil_path);

#endif //TOKAMAK_H
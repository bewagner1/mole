/**
 * @file jacobian.h
 * 
 * @brief Jacobian Metrics for Curvilinear Operators
 * 
 * @date 2026/4/9
 */

#ifndef JACOBIAN_H
#define JACOBIAN_H

#include "gradient.h"
#include "interpolFtoC.h"
#include "utils.h"
#include <cassert>

/**
 * @brief Returns the 2-D Jacobian metrics (J, xe, xn, ye, yn) at the centers
 * 
 * @param X x-coordinates of Omega at the centers
 * @param Y y-coordinates of Omega at the centers
 * @param k Order of accuracy
 * @param m Number of cells in the xi-direction
 * @param dx Step size in the xi-direction
 * @param n Number of cells in the eta-direction
 * @param dy Step size in the eta-direction
 * @param dc Dirichlet coefficients of left, right, bottom, and top boundaries
 * @param nc Neumann coefficients of left, right, bottom, and top boundaries
 */
mat jacobian(const vec& X, const vec& Y,
             u16 k, u32 m, Real dx, u32 n, Real dy,
             const vec& dc, const vec& nc);

/**
 * @brief Returns the 3-D Jacobian metrics (J, xe, xn, xk, ye, yn, yk, ze, zn, zk) at the centers
 * 
 * @param X x-coordinates of Omega at the centers
 * @param Y y-coordinates of Omega at the centers
 * @param Z z-coordinates of Omega at the centers
 * @param k Order of accuracy
 * @param m Number of cells in the xi-direction
 * @param dx Step size in the xi-direction
 * @param n Number of cells in the eta-direction
 * @param dy Step size in the eta-direction
 * @param o Number of cells in the kappa-direction
 * @param dz Step size in the kappa-direction
 * @param dc Dirichlet coefficients of left, right, bottom, top, front, and back boundaries
 * @param nc Neumann coefficients of left, right, bottom, top, front, and back boundaries
 */
mat jacobian(const vec& X, const vec& Y, const vec& Z,
             u16 k, u32 m, Real dx, u32 n, Real dy, u32 o, Real dz,
             const vec& dc, const vec& nc);

#endif //JACOBIAN_H
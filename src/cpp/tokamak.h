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
#include <fstream>

/**
 * @brief Default path to file containing x-coordinates of tokamak boundary
 */
#define XLIM_PATH ""

/**
 * @brief Default path to file containing y-coordinates of tokamak boundary
 */
#define YLIM_PATH ""

/**
 * @brief Default path to file containing coil information
 */
#define COIL_PATH ""

/**
 * @brief Alias for pi
 */
#define PI arma::datum::pi

/**
 * @brief Alias for mu0 (vacuum magnetic permeability)
 */
#define MU0 arma::datum::mu_0

/**
 * @brief Default number of points to sample along the plasma boundary
 */
#define NUM_PLASMA_BDRY 200

/**
 * 
 */
struct Point2D {
    Real x, y;
};

using Segment = std::pair<Point2D, Point2D>;

using Polyline = std::vector<Point2D>;

/**
 * 
 */
Point2D interp(Point2D a, Point2D b, Real va, Real vb, Real level);

/**
 * 
 */
std::vector<Segment> marchingSquares(const mat& psi, const mat& R, const mat& Z, Real level);

/**
 * 
 */
Real ptDist(Point2D a, Point2D b);

/**
 * 
 */
std::vector<Polyline> chainSegments(const std::vector<Segment>& segs, Real tol = 1e-9);

/**
 * 
 */
bool allClosed(const std::vector<Polyline>& chains, Real tol = 1e-6);

/**
 * 
 */
Real lastClosedContour(const mat& psi, const mat& R, const mat& Z, int iterations = 50);

/**
 * 
 */
Polyline sampleClosedContour(const Polyline& chain, int n);

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
 * 
 * @returns A weighted laplacian operator
 */
sp_mat weighted_laplacian(u16 k, u32 m, Real dx, u32 n, Real dy, Real alpha, Real beta);

/**
 * @brief Constructs the x or y coordinates of the plasma domain using a reference grid
 * 
 * @param k Order of accuracy
 * @param left coordinates of left boundary nodes, ordered bottom to top
 * @param right coordinates of right boundary nodes, ordered bottom to top
 * @param bottom coordinates of bottom boundary nodes, ordered left to right
 * @param top coordinates of top boundary nodes, ordered left to right
 * 
 * @returns A vector of coordiantes of the domain inside of left, right, bottom, and top
 */
vec plasma_grid(u16 k, const vec& left, const vec& right, const vec& bottom, const vec& top);

/**
 * @brief Constructs to x- or y- coordinates of the vacuum domain using a reference grid
 * 
 * @param k Order of accuracy
 * @param bottom coordinates of plasma boundary centers, ordered left to right
 * @param top coordinates of tokamak boundary centers, ordered left to right
 * 
 * @returns A vector of coordinates of the domain between bottom and top
 */
vec vacuum_grid(u16 k, u32 n, const vec& bottom, const vec& top);

/**
 * @brief Reads coil data from file
 * 
 * @param coil_path Path to coil file
 * 
 * @returns A matrix consisting of (R, Z, I) of a coil per row
 */
mat read_coils(const char* coil_path);

/**
 * @brief Reads boundary information from file
 * 
 * @param bdry_path Path to boundary file
 * 
 * @returns Vector containing coordinates of tokamak boundary
 */
vec read_boundary(const char* bdry_path);

/**
 * @brief Segments the boundary into 4 pieces for making a mesh
 * 
 * @param r_bdry R-coordinates of boundary
 * @param z_bdry Z-coordinates of boundary
 * @param left_r R-coordinates of left boundary piece
 * @param right_r R-coordinates of right boundary piece
 * @param bottom_r R-coordinates of bottom boundary piece
 * @param top_r R-coordinates of top boundary piece
 * @param left_z Z-coordinates of left boundary piece
 * @param right_z Z-coordinates of right boundary piece
 * @param bottom_z Z-coordinates of bottom boundary piece
 * @param top_z Z-coordinates of top boundary piece
 * @param m Number of cells in xi-direction
 * @param n Number of cells in eta-direction
 */
void segment_boundary(const vec& r_bdry, const vec& z_bdry, vec& left_r, vec& right_r, vec& bottom_r, vec& top_r, vec& left_z, vec& right_z, vec& bottom_z, vec& top_z, const u32 m, const u32 n);

/**
 * @brief 4th Order Lagrangian interpolation
 * 
 * @param m Number of cells in x-direction
 * @param n Number of cells in y-direction
 * 
 * @returns A 4th order Lagrangian interpolation operator from nodes to centers
 */
sp_mat interpolNodesToCentersCurv(u32 m, u32 n);

/**
 * @brief Finds the last closed flux surface of psi
 * 
 * @param R R-coordinates of domain
 * @param Z Z-coordinates of domain
 * @param psi Poloidal magnetic flux per radian in domain
 * 
 * @returns A matrix of points along the last closed flux surface
 */
mat get_last_closed_flux_surface(const mat& R, const mat& Z, const mat& psi);

/**
 * @brief Finds the points on or inside of the last closed flux surface
 * 
 * @param R R-coordinates of domain
 * @param Z Z-coordinates of domain
 * @param psi Polodial magnetic flux per radian of domain
 * 
 * @returns A vector of the indices that correspond to points on or inside of the last closed flux surface
 */
uvec get_plasma_indices(const vec& R, const vec& Z, const vec& psi, const u32 m, const u32 n);

#endif //TOKAMAK_H
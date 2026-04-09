/**
 * @file tokamak.cpp
 * 
 * @brief Implementation of functions for solving the Grad-Shafranov equation
 */

#include "tokamak.h"

// Used for reference grid
sp_mat weighted_laplacian(u16 k, u32 m, Real dx, u32 n, Real dy, Real alpha, Real beta)
{
    Divergence De(k, m, dx);
    Gradient Ge(k, m, dx);

    Divergence Dn(k, n, dy);
    Gradient Gn(k, n, dy);

    sp_mat Im = speye(m + 2, m + 2);
    Im.shed_col(0);
    Im.shed_col(m);

    sp_mat In = speye(n + 2, n + 2);
    In.shed_col(0);
    In.shed_col(m);

    sp_mat Dx = Utils::spkron(In, De);
    sp_mat Gx = Utils::spkron(In.t(), Ge);
    sp_mat Dy = Utils::spkron(Dn, Im);
    sp_mat Gy = Utils::spkron(Gn, Im.t());

    sp_mat A = (1.0 / (alpha * alpha)) * speye(Dx.n_cols, Gx.n_rows);
    sp_mat B = (1.0 / (beta * beta)) * speye(Dy.n_cols, Gy.n_rows);

    sp_mat L = Dx * A * Gx + Dy * B * Gy;
    return L;
}

// Generates x- or y- coordinates of plasma grid
vec plasma_grid(u16 k, vec left, vec right, vec bottom, vec top)
{
    u32 m = bottom.n_elem - 1;
    u32 n = top.n_elem - 1;

    Real dx = 1.0 / m;
    Real dy = 1.0 / n;

    Real alpha = 0.0;
    Real beta = 0.0;

    ivec dc = {1, 1, 1, 1};
    ivec nc = {0, 0, 0, 0};

    InterpolNtoC INC(k, m, n, dc, nc);
    InterpolCtoN ICN(k, m, n, dc, nc);

    sp_mat L = INC * weighted_laplacian(k, m, dx, n, dy, alpha, beta) * ICN;

    AddScalarBC::BC2D bc;
    bc.dc = {1.0, 1.0, 1.0, 1.0};
    bc.nc = {0.0, 0.0, 0.0, 0.0};
    bc.v[0] = left(span(1, n));
    bc.v[1] = right(span(1, n));
    bc.v[2] = bottom;
    bc.v[3] = top;

    vec b = zeros(L.n_rows);

    AddScalarBC::addScalarBC(L, b, k, m, dx, n, dy, bc);

#ifdef EIGEN
    vec sol = Utils::spsolve_eigen(L, b); 
#else
    vec sol = spsolve(L, b);
#endif

    return sol;
}

// Generates x- or y- coordinates of vacuum grid
vec vacuum_grid(u16 k, u32 n, vec bottom, vec top)
{
    Real alpha = 0.0;
    Real beta = 0.0;

    u32 m = top.n_elem;

    Real dx = 1.0 / m;
    Real dy = 1.0 / n;

    ivec dc = {1, 1};
    ivec nc = {0, 0};

    Divergence De(k, m, dx, nc, nc); // Periodic in xi
    Gradient Ge(k, m, dx, nc, nc);   // Periodic in xi

    Divergence Dn(k, n, dy, dc, nc); // Nonperiodic in eta
    Gradient Gn(k, n, dy, dc, nc);   // Nonperiodic in eta

    sp_mat Im = speye(m, m);
    sp_mat In = speye(n + 2, n + 2);
    In.shed_col(0);
    In.shed_col(n);

    sp_mat Dx = Utils::spkron(In, De);
    sp_mat Gx = Utils::spkron(In.t(), Ge);

    sp_mat Dy = Utils::spkron(Dn, Im);
    sp_mat Gy = Utils::spkron(Gn, Im); // Im.t() == Im

    sp_mat A = (1.0 / (alpha * alpha)) * speye(Dx.n_cols, Gx.n_rows);
    sp_mat B = (1.0 / (beta * beta)) * speye(Dy.n_cols, Gy.n_cols);

    sp_mat L = Dx * A * Gx + Dy * B * Gy;

    vec b = zeros(L.n_rows);

    AddScalarBC::BC2D bc;
    bc.dc = {0.0, 0.0, 1.0, 1.0};
    bc.nc = {0.0, 0.0, 0.0, 0.0};
    bc.v[0] = 0.0;
    bc.v[1] = 0.0;
    bc.v[2] = bottom;
    bc.v[3] = top;

    AddScalarBC::addScalarBC(L, b, k, m, dx, n, dy, bc);

#ifdef EIGEN
    vec sol = Utils::spsolve_eigen(L, b);
#else
    vec sol = spsolve(L, b);
#endif

    return sol;
}

// Reads coil information
mat read_coils(const char* coil_path)
{

}

// 
sp_mat interpolNodesToCentersCurv(u32 m, u32 n)
{

}

// 
vec get_boundary(vec domain, u32 m, u32 n)
{
    return get_boundary((mat)reshape(domain, n, m));
}

// 
vec get_boundary(mat domain)
{

}

// 
vec get_separatrix(vec plasma_r, vec plasma_z, vec plasma_p, vec vacuum_r, vec vacuum_z, vec vacuum_p, const u32 num_plasma_bdry)
{
    
}
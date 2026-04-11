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
vec plasma_grid(u16 k, const vec& left, const vec& right, const vec& bottom, const vec& top)
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
vec vacuum_grid(u16 k, u32 n, const vec& bottom, const vec& top)
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

// Reads boundary information
vec read_boundary(const char* bdry_path)
{

}

// Segments boundary into 4 pieces, left, right, bottom, and top
void segment_boundary(const vec& r_bdry, const vec& z_bdry, vec& left_r, vec& right_r, vec& bottom_r, vec& top_r, vec& left_z, vec& right_z, vec& bottom_z, vec& top_z, const u32 m, const u32 n)
{
    uword idx = index_min(r_bdry + z_bdry);

    vec r_bdry_shifted = circshift(r_bdry, idx);
    vec z_bdry_shifted = circshift(z_bdry, idx);

    left_r = r_bdry_shifted(span(0, n));
    left_z = z_bdry_shifted(span(0, n));

    top_r = r_bdry_shifted(span(n, n + m + 1));
    top_z = z_bdry_shifted(span(n, n + m + 1));

    right_r = r_bdry_shifted(span(n + m + 1, n + m + n + 2));
    right_z = z_bdry_shifted(span(n + m + 1, n + m + n + 2));

    bottom_r = r_bdry_shifted(span(n + m + n + 2, n + m + n + m + 3));
    bottom_z = z_bdry_shifted(span(n + m + n + 2, n + m + n + m + 3));

    right_r = flipud(right_r); // Order needs to be bottom to top
    right_z = flipud(right_z);

    bottom_r = flipud(bottom_r); // Order needs to be left to right
    bottom_z = flipud(bottom_z);
}

// 
sp_mat interpolNodesToCentersCurv(u32 m, u32 n)
{
    sp_mat Ix(m + 2, m + 1);
    sp_mat Iy(n + 2, n + 1);

    Ix(0, 0) = 1.0;
    Ix(1, 0) = 5.0 / 16.0;
    Ix(1, 1) = 15.0 / 16.0;
    Ix(1, 2) = -5.0 / 16.0;
    Ix(1, 3) = 1.0 / 16.0;
    for (int i = 2; i < m; ++i) {
        Ix(i, i - 2) = -1.0 / 16.0;
        Ix(i, i - 1) = 9.0 / 16.0;
        Ix(i, i) = 9.0 / 16.0;
        Ix(i, i + 1) = -1.0 / 16.0;
    }
    Ix(m, m - 3) = 1.0 / 16.0;
    Ix(m, m - 2) = -5.0 / 16.0;
    Ix(m, m - 1) = 15.0 / 16.0;
    Ix(m, m) = 5.0 / 16.0;
    Ix(m + 1, m) = 1.0;

    Iy(0, 0) = 1.0;
    Iy(1, 0) = 5.0 / 16.0;
    Iy(1, 1) = 15.0 / 16.0;
    Iy(1, 2) = -5.0 / 16.0;
    Iy(1, 3) = 1.0 / 16.0;
    for (int i = 2; i < n; ++i) {
        Iy(i, i - 2) = -1.0 / 16.0;
        Iy(i, i - 1) = 9.0 / 16.0;
        Iy(i, i) = 9.0 / 16.0;
        Iy(i, i + 1) = -1.0 / 16.0;
    }
    Iy(n, n - 3) = 1.0 / 16.0;
    Iy(n, n - 2) = -5.0 / 16.0;
    Iy(n, n - 1) = 15.0 / 16.0;
    Iy(n, n) = 5.0 / 16.0;
    Iy(n + 1, n) = 1.0;

    return Utils::spkron(Iy, Ix);
}

// 
mat get_last_closed_flux_surface(const vec& R, const vec& Z, const vec& psi)
{
    /* TODO: Find last closed flux surface */
    const u32 num_plasma_bdry = 200;
}

// 
uvec get_plasma_indices(const vec& R, const vec& Z, const vec& psi)
{
    mat plasma_bdry = get_last_closed_flux_surface(R, Z, psi);

    /* TODO: Find indices of R and Z that fall within plasma_bdry */
}
/**
 * @file grad_shaf.cpp
 * 
 * @brief Solves the Grad-Shafranov equation
 * 
 * @date 2026/04/10
 */

#include "mole.h"
#include "tokamak.h"

/**
 * @brief Provides an initial guess for psi
 * 
 * @param R R-coordinates of domain
 * @param Z Z-coordinates of domain
 * 
 * @returns Initial guess of psi
 */
vec initial_guess(const vec& R, const vec& Z);

/**
 * @brief Solves the Grad-Shafranov equation
 * 
 * @param R R-coordinates of domain
 * @param Z Z-coordinates of domain
 * @param k Order of accuracy
 * @param m Number of cells in xi-direction
 * @param n Number of cells in eta-direction
 * 
 * @returns Solution to the Grad-Shafranov equation
 */
vec solve(const vec& R, const vec& Z, const u16 k, const u32 m, const u32 n);

/**
 * @brief Calculates the source term for the next iteration
 * 
 * @param R R-coordinates of domain
 * @param Z Z-coordinates of domain
 * @param psi Solution at last iteration
 * @param RHS Right hand side source term
 * @param m Number of cells in xi-direction
 * @param n Number of cells in eta-direction
 */
void apply_rhs(const vec& R, const vec& Z, const vec& psi, vec& RHS, const u32 m, const u32 n);

/**
 * @brief Adds the poloidal coil flux to psi
 * 
 * @param R R-coordinates of domain
 * @param Z Z-coordinates of domain
 * @param psi Poloidal magnetic flux per radian of domain
 * @param coils Matrix of coil information. Each row is r, z, I for a coil
 */
void coil_contribution(const vec& R, const vec& Z, vec& psi, const mat& coils);

int main(int argc, char* argv[])
{
    if (/* check proper number of arguments */) {

    }

    // Parameters and things
    constexpr u16 k = 8;                   // Order of accuracy
    constexpr u32 num_bdry = 200;          // Should this be an argument?
    constexpr u32 m = 35 * num_bdry / 100; // Number of cells in xi-direction
    constexpr u32 n = 15 * num_bdry / 100; // Number of cells in eta-direction

    /**********************
        Read Input Data
    **********************/

    vec r_bdry, z_bdry;
    mat coils;
    if (/*something to do with argc*/) {

        // Tokamak boundary points
        r_bdry = read_boundary(argv[0]);
        z_bdry = read_boundary(argv[1]);

        // Coil data
        coils = read_coils(argv[2]);

        // efit.input things

    } else {

        // Tokamak boundary points
        r_bdry.load(XLIM_PATH);
        z_bdry.load(YLIM_PATH);
        
        // Coil data
        coils.load(COIL_PATH);

        // efit.input things

    }

    /* Interpolate boundary ? */

    vec left_r(n + 1), right_r(n + 1), bottom_r(m + 1), top_r(m + 1);
    vec left_z(n + 1), right_z(n + 1), bottom_z(m + 1), top_z(m + 1);

    segment_boundary(r_bdry, z_bdry, left_r, right_r, bottom_r, top_r, left_z, right_z, bottom_z, top_z, m, n);

    vec R = plasma_grid(k, left_r, right_r, bottom_r, top_r);
    vec Z = plasma_grid(k, left_z, right_z, bottom_z, top_z);

    sp_mat INC = interpolNodesToCentersCurv(m, n);
    R = INC * R;
    Z = INC * Z;

    vec psi = solve(R, Z, k, m, n);
    coil_contribution(R, Z, psi, coils);

    /* Time dependency ? */

    // Write to solution to a file (for plotting elsewhere)
    if (/* something to do with argc */) {

    } else {
        mat rzp = join_rows(R, Z, psi);
        char buff[256];
        int _ = sprintf(buff, "rzp%um%un.csv", m, n);
        rzp.save(buff, arma::csv_ascii);
    }

    return 0;
}

// Gaussian
vec initial_guess(const vec& R, const vec& Z)
{
    Real rmag = (R.max() - R.min()) / 2.0;
    Real zmag = 0.0;
    Real rsig = 1.0;
    Real zsig = 1.0;
    vec psi = exp(-((R - rmag) % (R - rmag) / (2.0 * rsig * rsig) + (Z - zmag) % (Z - zmag) / (2.0 * zsig * zsig)));

    return psi;
}

// Picard Iterations
vec solve(const vec& R, const vec& Z, const u16 k, const u32 m, const u32 n)
{
    Real dx = 1.0 / m;
    Real dy = 1.0 / n;
    vec dc = {1.0, 1.0, 1.0, 1.0};
    vec nc = {0.0, 0.0, 0.0, 0.0};

    GradCurv G(R, Z, k, m, dx, n, dy, dc, nc);
    DivCurv D(R, Z, k, m, dx, n, dy, dc, nc);

    InterpolCtoF Ix(k, m, dc.subvec(0, 2), nc.subvec(0, 2));
    sp_mat In = speye(n + 2, n + 2);

    // R^2 div( grad(psi) / R^2 ) = psi_RR + psi_ZZ - (2 / R) psi_R
    sp_mat L = D * G - (2 / R) % D * Utils::spkron(In, Ix);

    vec RHS(L.n_rows);

    AddScalarBC::BC2D bc;
    bc.dc = dc;                        // Tokamak Dirichlet coefficients
    bc.nc = nc;                        // Tokamak Neumann coefficients
    bc.v[0] = vec(n, fill::zeros);     // Left tokamak boundary
    bc.v[1] = vec(n, fill::zeros);     // Right tokamak boundary
    bc.v[2] = vec(m + 2, fill::zeros); // Bottom tokamak boundary
    bc.v[3] = vec(m + 2, fill::zeros); // Top tokamak boundary

    vec psi = initial_guess(R, Z);
    vec psi0 = psi;

    Real e = 1.0;
    Real tol = 1e-2;
    u32 i = 0;
    while (e > tol) {

        apply_rhs(R, Z, psi, RHS, m, n);

        AddScalarBC::addScalarBC(L, RHS, k, m, dx, n, dy, bc);

#ifdef EIGEN
        psi = Utils::spsolve_eigen(L, RHS);
#else
        spsolve(psi, L, RHS);
#endif

        e = norm(psi - psi0);
        psi0 = psi;
        if (++i > 1000000) {
            std::cerr << "No solution reached after 1,000,000 iterations" << std::endl;
            exit(1);
        }
    }

    std::cout << "Solution reached after " << i << " iterations" << std::endl;
    return psi;
}

// Source term
void apply_rhs(const vec& R, const vec& Z, const vec& psi, vec& RHS, const u32 m, const u32 n)
{
    RHS = RHS.zeros();
    uvec pidx = get_plasma_indices(R, Z, psi, m, n);
    Real simag = min(psi(pidx));
    Real sibry = max(psi(pidx)); // Is this the best ?
    vec upsi = (psi(pidx) - simag) / (sibry - simag);
    RHS(pidx) = -MU0 * R(pidx) % (upsi) - (upsi); // TODO: Figure out parameterization of dp/dpsi and dF^2 / dpsi
}

// Coils
void coil_contribution(const vec& R, const vec& Z, vec& psi, const mat& coils)
{
    vec cc = zeros(R.n_elem);
    int n = 15;
    int r, z, I;
    for (int i = 0; i < coils.n_rows; ++i) {
        r = coils(i, 0);
        z = coils(i, 2);
        I = coils(i, 3);
        for (int k = 0; k < n; ++k) {
            cc += cos(2 * PI * k / n) / sqrt(pow(R - r * cos(2 * PI * k / n), 2) + pow(- r * sin(2 * PI * k / n), 2) + (Z - z) % (Z - z));
        }
        psi += (2 * PI * I * 1e-7 / n) * r * R % cc; // mu0 / 4 pi = 1e-7
    }
}
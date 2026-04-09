/**
 * 
 */

#include "mole.h"
#include "tokamak.h"

/**
 * 
 */
sp_mat interpolNodesToCentersCurv(u32 m, u32 n);

/**
 * 
 */
vec get_boundary(vec domain);

/**
 * 
 */
vec solve();

int main(int argc, char* argv[])
{
    if (/* chech proper number of arguments */) {

    }

    // Parameters and things
    constexpr u16 k = 8;
    constexpr u32 mp = 70; // For plasma domain
    constexpr u32 np = 30; // For plasma domain
    constexpr u32 nv = 25; // For vacuum domain (mv is 2 * (mp + np) + 4)

    ivec pdc = {1, 1, 1, 1};
    ivec pnc = {0, 0, 0, 0};
    ivec vdc = {0, 0, 1, 1};
    ivec vnc = {0, 0, 0, 0};

    /**********************
        Read Input Data
    **********************/

    vec tok_r_bdry, tok_z_bdry;
    mat coils;
    if (/*something to do with argc*/) {

        // Tokamak boundary points

        // Coil data

        // efit.input things

    } else {

        // Tokamak boundary points
        tok_r_bdry.load(XLIM_PATH);
        tok_z_bdry.load(YLIM_PATH);
        
        // Coil data
        coils.load(COIL_PATH);

        // efit.input things

    }

    /**********************
        Nonlinear Solve
    **********************/

    vec plasma_left_r(np), plasma_right_r(np), plasma_bottom_r(mp), plama_top_r(mp);
    vec plasma_left_z(np), plasma_right_z(np), plasma_bottom_z(mp), plama_top_z(mp);

    /* initial guess of plama boundary */

    vec plasma_r_nodes = plasma_grid(k, plasma_left_r, plasma_right_r, plasma_bottom_r, plama_top_r);
    vec plasma_z_nodes = plasma_grid(k, plasma_left_z, plasma_right_z, plasma_bottom_z, plama_top_z);

    sp_mat INC = interpolNodesToCentersCurv(mp, np);
    vec plasma_r = INC * plasma_r_nodes;
    vec plasma_z = INC * plasma_z_nodes;

    vec plasma_r_bdry = get_boundary(plasma_r);
    vec plasma_z_bdry = get_boundary(plasma_z);
    vec vacuum_r = vacuum_grid(k, nv, plasma_r_bdry, tok_r_bdry);
    vec vacuum_z = vacuum_grid(k, nv, plasma_z_bdry, tok_z_bdry);

    // Picard Iterations

    // Coil Contributions

    // Find separatrix


    // Time dependency ?

    return 0;
}

// Function implementations

sp_mat interpolNodesToCentersCurv(u32 m, u32 n)
{

}

vec get_boundary(vec domain)
{

}

vec solve()
{
    
}
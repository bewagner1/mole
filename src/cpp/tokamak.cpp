/**
 * @file tokamak.cpp
 * 
 * @brief Implementation of functions for solving the Grad-Shafranov equation
 */

#include "tokamak.h"

// 
Point2D interp(const Point2D& a, const Point2D& b, Real va, Real vb, Real level)
{
    Real t = (level - va) / (vb - va);
    return {a.x + t * (b.x - a.x), a.y + t * (b.y - a.y)};
}

// 
std::vector<Segment> marchingSquares(const mat& psi, const mat& R, const mat& Z, Real level)
{
    std::vector<Segment> segs;
    uword nr = psi.n_rows;
    uword nc = psi.n_cols;
    int edges[4][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};

    for (uword i = 0; i < nr; ++i) {
        for (uword j = 0; j < nc; ++j) {
            Real v[4] = { psi(i, j),
                          psi(i, j + 1),
                          psi(i + 1, j + 1),
                          psi(i + 1, j)
                        };

            Point2D p[4] = { { R(i, j), Z(i, j) },
                             { R(i, j + 1), Z(i, j + 1) },
                             { R(i + 1, j + 1), Z(i + 1, j + 1) },
                             { R(i + 1, j), Z(i + 1, j) }
                           };

            int idx = 0;
            for (int k = 0; k < 4; ++k) {
                if (v[k] >= level) idx |= (1 << k);
            }
            if (idx == 0 || idx == 15) continue;

            std::vector<Point2D> crossings;
            for (auto& e : edges) {
                int a = e[0];
                int b = e[1];
                bool aAbove = (v[a] >= level);
                bool bAbove = (v[b] >= level);
                if (aAbove != bAbove) crossings.push_back(interp(p[a], p[b], v[a], v[b], level));
            }

            if (crossings.size() >= 2) segs.push_back({crossings[0], crossings[1]});
            if (crossings.size() == 4) segs.push_back({crossings[2], crossings[3]});
        }
    }

    return segs;
}

// 
Real ptDist(Point2D a, Point2D b)
{
    return std::hypot(a.x - b.x, a.y - b.y);
}

// 
std::vector<Polyline> chainSegments(const std::vector<Segment>& segs, Real tol = 1e-9)
{
    std::vector<bool> used(segs.size(), false);
    std::vector<Polyline> chains;

    for (size_t i = 0; i < segs.size(); ++i) {
        if (used[i]) continue;
        Polyline chain = {segs[i].first, segs[i].second};
        used[i] = true;
        bool grew = true;
        while (grew) {
            grew = false;
            for (size_t j = 0; j < segs.size(); ++j) {
                if (used[j]) continue;
                Point2D& tail = chain.back();
                if (ptDist(tail, segs[j].first) < tol) {
                    chain.push_back(segs[j].second);
                    used[j] = true;
                    grew = true;
                } else if (ptDist(tail, segs[j].second) < tol) {
                    chain.push_back(segs[j].first);
                    used[j] = true;
                    grew = true;
                }
            }
        }
        chains.push_back(chain);
    }

    return chains;
}

// 
bool allClosed(const std::vector<Polyline>& chains, Real tol = 1e-6)
{
    if (chains.empty()) return false;
    for (const auto& c : chains) {
        if (c.size() < 3 || ptDist(c.front(), c.back()) > tol) return false;
    }
    return true;
}

// 
Real lastClosedContour(const mat& psi, const mat& R, const mat& Z, int iterations = 50)
{
    Real lo = psi.min() + 1e-9;
    Real hi = psi.max() - 1e-9;

    for (int iter = 0; iter < iterations; ++iter) {
        Real mid = 0.5 * (lo + hi);
        auto segs = marchingSquares(psi, R, Z, mid);
        auto chains = chainSegments(segs);
        if (allClosed(chains)) {
            lo = mid;
        } else {
            hi = mid;
        }
    }

    return lo;
}

// 
Polyline sampleClosedContour(const Polyline& chain, int n)
{
    int m = chain.size();
    std::vector<Real> arcLen(m, 0.0);
    for (int i = 1; i < m; ++i) {
        arcLen[i] = arcLen[i - 1] + ptDist(chain[i - 1], chain[i]);
    }

    Real totalLen = arcLen.back();

    Polyline result;
    result.reserve(n);

    for (int k = 0; k < n; ++k) {
        Real target = totalLen * k / n;

        auto it = std::lower_bound(arcLen.begin(), arcLen.end(), target);
        int idx = std::max(0, (int)std::distance(arcLen.begin(), it) - 1);

        Real segLen = arcLen[idx + 1] - arcLen[idx];
        Real t = (segLen > 1e-12) ? (target - arcLen[idx]) / segLen : 0.0;

        result.push_back({chain[idx].x + t * (chain[idx + 1].x - chain[idx].x),
                          chain[idx].y + t * (chain[idx + 1].y - chain[idx].y)});
    }

    return result;
}

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
    bc.v[0] = left(arma::span(1, n));
    bc.v[1] = right(arma::span(1, n));
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

    left_r = r_bdry_shifted(arma::span(0, n)); // n + 1 points
    left_z = z_bdry_shifted(arma::span(0, n));

    top_r = r_bdry_shifted(arma::span(n, n + m)); // m + 1 points
    top_z = z_bdry_shifted(arma::span(n, n + m));

    right_r = r_bdry_shifted(arma::span(n + m, n + m + n)); // n + 1 points
    right_z = z_bdry_shifted(arma::span(n + m, n + m + n));

    bottom_r(arma::span(0, m - 1)) = r_bdry_shifted(arma::span(n + m + n, n + m + n + m - 1)); // m points
    bottom_z(arma::span(0, m - 1)) = z_bdry_shifted(arma::span(n + m + n, n + m + n + m - 1));

    bottom_r(m) = r_bdry_shifted(0); // m + 1 point
    bottom_z(m) = z_bdry_shifted(0);

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
mat get_last_closed_flux_surface(const mat& R, const mat& Z, const mat& psi)
{
    Real sibry = lastClosedContour(psi, R, Z);

    auto segs = marchingSquares(psi, R, Z, sibry);
    auto chains = chainSegments(segs);

    Polyline sampled = sampleClosedContour(chains[0], NUM_PLASMA_BDRY);

    mat surf(NUM_PLASMA_BDRY, 2);
    for (int i = 0; i < NUM_PLASMA_BDRY; ++i) {
        surf(i, 0) = sampled[i].x;
        surf(i, 1) = sampled[i].y;
    }

    return surf;
}

// 
uvec get_plasma_indices(const vec& R, const vec& Z, const vec& psi)
{
    mat plasma_bdry = get_last_closed_flux_surface(R, Z, psi);
    vec plasma_r = plasma_bdry.col(0);
    vec plasma_z = plasma_bdry.col(1);

    int n = plasma_bdry.n_rows;
    int j = n - 1;

    Col<bool> inside = Col<bool>(R.n_elem, fill::zeros);
    bool intersect;
    Real ri, zi, rj, zj;
    for (int i = 0; i < n; ++i) {
        ri = plasma_r(i);
        zi = plasma_z(i);
        rj = plasma_r(j);
        zj = plasma_z(j);

        for (int k = 0; k < R.n_elem; ++k) {
            intersect = ((zi > Z(k)) != (zj > Z(k))) && (R(k) < (rj - ri) * (Z(k) - zi) / (zj - zi) + ri);
            inside[k] = inside[k] ^ intersect;
        }
        
        j = i;
    }

    return find(inside);
}
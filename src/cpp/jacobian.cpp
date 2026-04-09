/**
 * @file jacobian.cpp
 * 
 * @brief Jacobian Metrics for Curvilinear Operators
 * 
 * @date 2026/4/9
 */

#include "jacobian.h"

// 2-D Jacobian
mat jacobian(const vec& X, const vec& Y, u16 k, u32 m, Real dx, u32 n, Real dy, const vec& dc, const vec& nc)
{
    assert(dc.n_elem == 4);
    assert(nc.n_elem == 4);

    Gradient Ge(k, m, dx, dc.subvec(0, 2), nc.subvec(0, 2));
    Gradient Gn(k, n, dy, dc.subvec(2, 2), dc.subvec(2, 2));

    InterpolFtoC IFCx(k, m, dc.subvec(0, 2), nc.subvec(0, 2));
    InterpolFtoC IFCy(k, n, dc.subvec(2, 2), nc.subvec(2, 2));

    sp_mat Im = speye(m + 2, m + 2);
    sp_mat In = speye(n + 2, n + 2);

    Ge = Utils::spkron(In, IFCx) * Utils::spkron(In, Ge);
    Gn = Utils::spkron(IFCy, Im) * Utils::spkron(Gn, Im);

    vec xe = Ge * X;
    vec xn = Gn * X;
    vec ye = Ge * Y;
    vec yn = Gn * Y;
    vec J = xe % yn - xn % ye;

    mat Jmat = join_rows(J, xe, xn, ye); // Can only do 4 at a time
    Jmat = join_rows(Jmat, yn);

    return Jmat;
}

// 3-D Jacobian
mat jacobian(const vec& X, const vec& Y, const vec& Z, u16 k, u32 m, Real dx, u32 n, Real dy, u32 o, Real dz, const vec& dc, const vec& nc)
{
    assert(dc.n_elem == 6);
    assert(nc.n_elem == 6);

    Gradient Ge(k, m, dx, dc.subvec(0, 2), nc.subvec(0, 2));
    Gradient Gn(k, n, dy, dc.subvec(2, 2), nc.subvec(2, 2));
    Gradient Gk(k, o, dz, dc.subvec(4, 2), nc.subvec(4, 2));

    InterpolFtoC IFCx(k, m, dc.subvec(0, 2), nc.subvec(0, 2));
    InterpolFtoC IFCy(k, n, dc.subvec(2, 2), nc.subvec(2, 2));
    InterpolFtoC IFCz(k, o, dc.subvec(4, 2), nc.subvec(4, 2));

    sp_mat Im = speye(m + 2, m + 2);
    sp_mat In = speye(n + 2, n + 2);
    sp_mat Io = speye(o + 2, o + 2);

    Ge = Utils::spkron(Utils::spkron()) * Utils::spkron(Utils::spkron());
    Gn = Utils::spkron(Utils::spkron()) * Utils::spkron(Utils::spkron());
    Gk = Utils::spkron(Utils::spkron()) * Utils::spkron(Utils::spkron());

    vec xe = Ge * X;
    vec xn = Gn * X;
    vec xk = Gk * X;
    vec ye = Ge * Y;
    vec yn = Gn * Y;
    vec yk = Gk * Y;
    vec ze = Ge * Z;
    vec zn = Gn * Z;
    vec zk = Gk * Z;
    vec J = xe % (yn % zk - yk % zn) - xn % (ye % zk - yk % ze) + xk % (ye % zn - yn % ze);

    mat Jmat = join_rows(J, xe, xn, xk);
    Jmat = join_rows(Jmat, ye, yn, yk);
    Jmat = join_rows(Jmat, ze, zn, zk);

    return Jmat;
}
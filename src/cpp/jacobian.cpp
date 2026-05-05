/*
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (c) 2008-2024 San Diego State University Research Foundation
 * (SDSURF).
 * See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html for details.
 */

/**
 * @file jacobian.cpp
 * 
 * @brief Jacobian Metrics
 * 
 * @date 2026/05/05
 */

#include "jacobian.h"

// 2-D Jacobian
arma::mat jacobian(u16 k, const arma::vec& X, const arma::vec& Y,
                   u32 m, Real dx, u32 n, Real dy,
                   const arma::ivec& dc, const arma::ivec& nc)
{
    assert(dc.n_elem == 4);
    assert(nc.n_elem == 4);
    assert(X.n_elem == Y.n_elem);

    Gradient Ge(k, m, dx, dc.subvec(0, 1), nc.subvec(0, 1));
    Gradient Gn(k, n, dy, dc.subvec(2, 3), nc.subvec(2, 3));

    InterpolFtoC IFCx(k, m, dc.subvec(0, 1), nc.subvec(0, 1));
    InterpolFtoC IFCy(k, n, dc.subvec(2, 3), nc.subvec(2, 3));

    arma::sp_mat Im;
    if (dc(0) == 0 && dc(1) == 0 && nc(0) == 0 && nc(1) == 0) {
        Im = arma::speye(m, m);
    } else {
        Im = arma::speye(m + 2, m + 2);
    }

    arma::sp_mat In;
    if (dc(2) == 0 && dc(3) == 0 && nc(2) == 0 && nc(3) == 0) {
        In = arma::speye(n, n);
    } else {
        In = arma::speye(n + 2, n + 2);
    }

    u32 num_centers = X.n_elem;

    arma::sp_mat G1 = Utils::spkron(In, IFCx) * Utils::spkron(In, Ge);
    arma::sp_mat G2 = Utils::spkron(IFCy, Im) * Utils::spkron(Gn, Im);

    arma::sp_mat G = Utils::spjoin_cols(G1, G2);

    arma::vec x_metrics = G * X;
    arma::vec y_metrics = G * Y;

    arma::vec xe = x_metrics.subvec(0, num_centers - 1);
    arma::vec xn = x_metrics.subvec(num_centers, 2 * num_centers - 1);
    arma::vec ye = y_metrics.subvec(0, num_centers - 1);
    arma::vec yn = y_metrics.subvec(num_centers, 2 * num_centers - 1);
    arma::vec J = xe % yn - xn % ye;

    arma::mat Jmat = arma::join_rows(J, xe, xn, ye); // Can only do 4 at a time
    return arma::join_rows(Jmat, yn);
}

// 3-D Jacobian
arma::mat jacobian(u16 k, const arma::vec& X, const arma::vec& Y, const arma::vec& Z,
                   u32 m, Real dx, u32 n, Real dy, u32 o, Real dz,
                   const arma::ivec& dc, const arma::ivec& nc)
{
    assert(dc.n_elem == 6);
    assert(nc.n_elem == 6);
    assert(X.n_elem == Y.n_elem && X.n_elem == Z.n_elem);

    Gradient Ge(k, m, dx, dc.subvec(0, 1), nc.subvec(0, 1));
    Gradient Gn(k, n, dy, dc.subvec(2, 3), nc.subvec(2, 3));
    Gradient Gk(k, o, dz, dc.subvec(4, 5), nc.subvec(4, 5));

    InterpolFtoC IFCx(k, m, dc.subvec(0, 1), nc.subvec(0, 1));
    InterpolFtoC IFCy(k, n, dc.subvec(2, 3), nc.subvec(2, 3));
    InterpolFtoC IFCz(k, o, dc.subvec(4, 5), nc.subvec(4, 5));

    arma::sp_mat Im;
    if (dc(0) == 0 && dc(1) == 0 && nc(0) == 0 && nc(1) == 0) {
        Im = arma::speye(m, m);
    } else {
        Im = arma::speye(m + 2, m + 2);
    }

    arma::sp_mat In;
    if (dc(2) == 0 && dc(3) == 0 && nc(2) == 0 && nc(3) == 0) {
        In = arma::speye(n, n);
    } else {
        In = arma::speye(n + 2, n + 2);
    }

    arma::sp_mat Io;
    if (dc(4) == 0 && dc(5) == 0 && nc(4) == 0 && nc(5) == 0) {
        In = arma::speye(o, o);
    } else {
        In = arma::speye(o + 2, o + 2);
    }

    u32 num_centers = X.n_elem;

    arma::sp_mat G1 = Utils::spkron(Utils::spkron(Io, In), IFCx)
                    * Utils::spkron(Utils::spkron(Io, In), Ge);
    arma::sp_mat G2 = Utils::spkron(Utils::spkron(Io, IFCy), Im)
                    * Utils::spkron(Utils::spkron(Io, Gn), Im);
    arma::sp_mat G3 = Utils::spkron(Utils::spkron(IFCz, In), Im)
                    * Utils::spkron(Utils::spkron(Gk, In), Im);

    arma::sp_mat G = Utils::spjoin_cols(Utils::spjoin_cols(G1, G2), G3);

    arma::vec x_metrics = G * X;
    arma::vec y_metrics = G * Y;
    arma::vec z_metrics = G * Z;

    arma::vec xe = x_metrics.subvec(0, num_centers - 1);
    arma::vec xn = x_metrics.subvec(num_centers, 2 * num_centers - 1);
    arma::vec xk = x_metrics.subvec(2 * num_centers, 3 * num_centers - 1);

    arma::vec ye = y_metrics.subvec(0, num_centers - 1);
    arma::vec yn = y_metrics.subvec(num_centers, 2 * num_centers - 1);
    arma::vec yk = y_metrics.subvec(2 * num_centers, 3 * num_centers - 1);

    arma::vec ze = z_metrics.subvec(0, num_centers - 1);
    arma::vec zn = z_metrics.subvec(num_centers, 2 * num_centers - 1);
    arma::vec zk = z_metrics.subvec(2 * num_centers, 3 * num_centers - 1);

    arma::vec J = xe % (yn % zk - yk % zn)
                - xn % (ye % zk - yk % ze)
                + xk % (ye % zn - yn % ze);

    arma::mat Jmat = arma::join_rows(J, xe, xn, xk);
    Jmat = arma::join_rows(Jmat, ye, yn, yk);
    return arma::join_rows(Jmat, ze, zn, zk);

}
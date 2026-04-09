/**
 * 
 */

#include "gradcurv.h"

// 2-D Implementation
GradCurv::GradCurv(const vec& X, const vec& Y, u16 k, u32 m, Real dx, u32 n, Real dy, const vec& dc, const vec& nc)
{
    assert(dc.n_elem == 4);
    assert(nc.n_elem == 4);

    Gradient Ge(k, m, dx, dc.subvec(0, 2), nc.subvec(0, 2));
    Gradient Ge(k, n, dy, dc.subvec(2, 2), nc.subvec(2, 2));

    InterpolFtoC IFCx(k, m, dc.subvec(0, 2), nc.subvec(0, 2));
    InterpolFtoC IFCy(k, n, dc.subvec(2, 2), nc.subvec(2, 2));

    InterpolCtoF ICFx(k, m, dc.subvec(0, 2), nc.subvec(0, 2));
    InterpolCtoF ICFy(k, n, dc.subvec(2, 2), nc.subvec(2, 2));

    sp_mat Im = speye(m + 2, m + 2);
    sp_mat In = speye(n + 2, n + 2);

    Ge = Utils::spkron(In, IFCx) * Utils::spkron(In, Ge);
    Gn = Utils::spkron(IFCy, Im) * Utils::spkron(Gn, Im);

    mat Jmat = jacobian(X, Y, k, m, dx, n, dy, dc, nc);
    vec J = Jmat.col(0);

    sp_mat Gx = (Jmat.col(4) / J) % Ge - (Jmat.col(3) / J) % Gn;
    sp_mat Gy = (Jmat.col(1) / J) % Gn - (Jmat.col(2) / J) % Ge;

    Gx = Utils::spkron(In, ICFx) * Gx;
    Gy = Utils::spkron(ICFy, Im) * Gy;

    if (m != n) {
        *this = Utils::spjoin_cols(Gx, Gy);
    } else {
        sp_mat A1(2, 1);
        sp_mat A2(2, 1);
        A1(0, 0) = A2(1, 0) = 1.0;
        *this = Utils::spkron(A1, Gx) + Utils::spkron(A2, Gy);
    }
}

// 3-D Implementation
GradCurv::GradCurv(const vec& X, const vec& Y, const vec& Z, u16 k, u32 m, Real dx, u32 n, Real dy, u32 o, Real dz, const vec& dc, const vec& nc)
{
    
}
/**
 * 
 */

#include "jacobian.h"

// 2-D Jacobian Metrics
mat jacobian(vec X, vec Y, u16 k, u32 m, Real dx, u32 n, Real dy, const ivec& dc, const ivec& nc)
{

    Gradient Ge(k, m, dx, dc, nc);
    Gradient Gn(k, n, dy, dc, nc);
    InterpolFtoC IFCx(k, m, dc, nc);
    InterpolFtoC IFCy(k, n, dc, nc);
    sp_mat Gx, Gy, Im, In;

    if (dc[0] == 0 && dc[1] == 0 && nc[0] == 0 && nc[1] == 0)
    {
        Im = speye(m, m);
    } else {
        Im = speye(m + 2, m + 2);
    }

    if (dc[3] == 0 && dc[4] == 0 && nc[3] == 0 && nc[4] == 0)
    {
        Im = speye(n, n);
    } else {
        Im = speye(n + 2, n + 2);
    }

    Gx = Utils::spkron(In, IFCx) * Utils::spkron(In, Ge);
    Gy = Utils::spkron(IFCy, Im) * Utils::spkron(Gn, Im);

    vec Xe = Gx * X;
    vec Xn = Gy * X;
    vec Ye = Gx * Y;
    vec Yn = Gy * Y;
    vec J = Xe % Yn - Xn % Ye;

    // Can only do 4 at a time
    mat Jmat = join_cols(J, Xe, Xn, Ye);
    Jmat = join_cols(Jmat, Yn);

    return Jmat;
}

// 3-D Jacobian Metrics
mat jacobian(vec X, vec Y, u16 k, u32 m, Real dx, u32 n, Real dy, u32 o, Real dz, const ivec& dc, const ivec& nc)
{
    std::cerr << "3-D Jacobian is not currently implemented" << std::endl;
    exit(1);
}
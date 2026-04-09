/**
 * 
 */

#ifndef GRADCURV_H
#define GRADCURV_H

#include "gradient.h"
#include "jacobian.h"
#include "interpolFtoC.h"
#include "interpolCtoF.h"
#include "utils.h"

class GradCurv : public sp_mat {

public:
    using sp_mat::operator=;

    /**
     * 
     */
    GradCurv(const vec& X, const vec& Y,
             u16 k, u32 m, Real dx, u32 n, Real dy,
             const vec& dc, const vec& nc);

    /**
     * 
     */
    GradCurv(const vec& X, const vec& Y, const vec& Z,
             u16 k, u32 m, Real dx, u32 n, Real dy, u32 o, Real dz,
             const vec& dc, const vec& nc);
};

#endif //GRADCURV_H
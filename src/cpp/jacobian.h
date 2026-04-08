/**
 * 
 */

#ifndef JACOBIAN_H
#define JACOBIAN_H

#include "gradient.h"
#include "interpolFtoC.h"

/**
 * 
 */
mat jacobian(vec X, vec Y, u16 k, u32 m, Real dx, u32 n, Real dy, const ivec& dc, const ivec& nc);

/**
 * 
 */
mat jacobian(vec X, vec Y, u16 k, u32 m, Real dx, u32 n, Real dy, u32 o, Real dz, const ivec& dc, const ivec& nc);

#endif //JACOBIAN_H
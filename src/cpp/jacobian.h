/*
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (c) 2008-2024 San Diego State University Research Foundation
 * (SDSURF).
 * See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html for details.
 */

/**
 * @file jacobian.h
 * 
 * @brief Jacobian Metrics
 * 
 * @date 2026/05/05
 */

#ifndef JACOBIAN_H
#define JACOBIAN_H

#include "gradient.h"
#include "interpolCtoF.h"
#include "interpolFtoC.h"
#include "utils.h"
#include <cassert>

/**
 * @brief 2-D Jacobian Metrics
 * 
 * @param k Order of accuracy
 * @param X X-coordinates of mesh centers
 * @param Y Y-coordinates of mesh centers
 * @param m Number of cells in xi-direction
 * @param dx Step size in xi-direction
 * @param n Number of cells in eta-direction
 * @param dy Step size in eta-direction
 * @param dc Robin coefficient a0; 4-element integer vector
 *           [left, right, bottom, top]. All-zero -> periodic BC.
 * @param nc Robin coefficient b0; 4-element integer vector 
 *           [left, right, bottom, top]. All-zero -> periodic BC.
 */
arma::mat jacobian(u16 k, const arma::vec& X, const arma::vec& Y,
                   u32 m, Real dx, u32 n, Real dy,
                   const arma::ivec& dc, const arma::ivec& nc);

/**
 * @brief 2-D Jacobian Metrics
 * 
 * @param k Order of accuracy
 * @param X X-coordinates of mesh centers
 * @param Y Y-coordinates of mesh centers
 * @param Z Z-coordinates of mesh centers
 * @param m Number of cells in xi-direction
 * @param dx Step size in xi-direction
 * @param n Number of cells in eta-direction
 * @param dy Step size in eta-direction
 * @param o Number of cells in kappa-direction
 * @param dz Step size in kappa-direction
 * @param dc Robin coefficient a0; 6-element integer vector 
 *           [left, right, bottom, top, front, back]. All-zero -> periodic BC.
 * @param nc Robin coefficient b0; 6-element integer vector 
 *           [left, right, bottom, top, front, back]. All-zero -> periodic BC.
 */
arma::mat jacobian(u16 k, const arma::vec& X, const arma::vec& Y, const arma::vec& Z,
                   u32 m, Real dx, u32 n, Real dy, u32 o, Real dz,
                   const arma::ivec& dc, const arma::ivec& nc);

#endif // JACOBIAN_H
using LinearAlgebra
using SparseArrays
using Plots

import MOLE: Operators, BCs

# Domain limits
west = 0.0
east = 1.0

k  = 6 # operator order of accuracy
m  = 2k + 1 # minimum number of cells to attain desired accuracy
dx = (east - west) / m # step length

path = joinpath(@__DIR__, "imgs") # Output path to store generated plots
mkpath(path)

# 1D staggered grid
grid = [west; (west + dx/2):dx:(east - dx/2); east]

# Mimetic Laplacian operator
L = Operators.lap(k, m, dx)

# Robin boundary conditions:  dc*u + nc*(du/dn) = v
dc = (1.0, 1.0)
nc = (1.0, 1.0)
v  = (0.0, 2 * exp(1.0))

bc = BCs.ScalarBC1D(dc, nc, v)

# RHS
rhs = exp.(grid)

# Apply boundary conditions and construct new system matrix/array
L0, rhs0 = BCs.addScalarBC!(sparse(L), rhs, bc, k, m, dx)

# Solve linear system
u = L0 \ rhs0

# Plot
p = scatter(grid, u; label="Approximated", marker=:circle, show = false)
plot!(p, grid, exp.(grid); label="Analytical", show = false)
plot!(
    p, 
    title="Poisson's equation with Robin BC",
    xlabel="x",
    ylabel="u(x)",
    legend=:topleft,
    show = false
)
Plots.png(p, joinpath(path, "elliptic1DaddScalarBC_output.png"))
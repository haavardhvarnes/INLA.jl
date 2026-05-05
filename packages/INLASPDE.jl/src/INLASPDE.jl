"""
    INLASPDE

SPDE–Matérn machinery for [`LatentGaussianModels`](@ref). Provides
FEM-assembled sparse precision matrices `Q(τ, κ)` on triangulated
meshes and exposes `SPDE2 <: AbstractLatentComponent` for use in
`LatentGaussianModel` stacks.

See `plans/plan.md` for the package roadmap and
`/plans/decisions.md` ADR-001, ADR-007, ADR-013 for load-bearing
design decisions.

### Module layout (populated incrementally per milestone)

- `assembly/` — per-element FEM matrices and `Q(α, τ, κ)` (M1).
- `components/` — `SPDE2 <: AbstractLatentComponent` (M2).
- `priors/` — PC-Matérn prior on `(range, σ)` (M2).
- `mesh/` — `inla_mesh_2d` and refinement helpers (M3).
- `projector.jl` — `MeshProjector` for observation-to-vertex maps (M4).
"""
module INLASPDE

using LinearAlgebra
using SparseArrays

using LatentGaussianModels
using LatentGaussianModels: AbstractLatentComponent, AbstractHyperPrior

using Meshes: Meshes
using DelaunayTriangulation: DelaunayTriangulation
using CoordRefSystems: CoordRefSystems
using SciMLOperators: SciMLOperators

# --- milestone includes (uncomment as each milestone lands) -----------
# M1 — FEM assembly
include("assembly/fem.jl")
include("assembly/lumping.jl")
include("assembly/precision.jl")
include("assembly/fem_1d.jl")

export assemble_fem_matrices, assemble_fem_matrices_1d
export lumped_mass, stiffness_squared
export FEMMatrices, spde_precision

# M2 — SPDE2 component + PC-Matérn prior
using GMRFs: GMRFs
include("priors/pc_matern.jl")
include("priors/gaussian_basis.jl")
include("components/spde2.jl")

export PCMatern, pc_matern_log_density
export GaussianBasisPrior, log_basis_prior_density
export SPDE2, spde_user_scale, spde_internal_scale

# M3 — Mesh generation
include("coercion.jl")
include("mesh/boundary.jl")
include("mesh/inla_mesh.jl")
include("mesh/inla_mesh_1d.jl")

export convex_hull_polygon, expand_polygon, cutoff_dedup
export INLAMesh, inla_mesh_2d, num_vertices, num_triangles
export INLAMesh1D, inla_mesh_1d, num_segments

# Phase M PR-2 — 1D SPDE component (depends on INLAMesh1D from M3)
include("components/spde1d.jl")
export SPDE1D

# Phase M PR-3 — non-stationary 2D SPDE component (depends on INLAMesh from M3)
include("components/spde2_nonstationary.jl")
export SPDE2NonStationary, spde_precision_nonstationary

# M4 — Projector
include("projector.jl")
include("projector_1d.jl")

export MeshProjector, MeshProjector1D, scimloperator

end # module

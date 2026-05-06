# Meuse zinc SPDE oracle — *prediction-side* fixture (Phase O PR-4).
#
# Companion to `meuse_spde.R`. The fit-side fixture there ships
# `summary_fixed`, `summary_hyperpar`, `mlik`, the fmesher mesh, and
# the observation projector. This script extends the same fit with
# R-INLA's posterior-mean *raster* over the Meuse hull, generated via
# `inla.mesh.project(mesh, loc = grid_centres)` — the canonical
# reference for `predict_raster`.
#
# Both implementations (R-INLA via fmesher, Julia via `MeshProjector`)
# are exact P1 barycentric interpolation on the same triangulation, so
# the Julia-side oracle test gates pixel-wise agreement at 1e-10 — the
# tightest "byte-for-byte vs R-INLA" gate in the repo.
#
# Output: fixtures/rasters/meuse_spde_predict.json
# (→ packages/INLASPDERasters.jl/test/oracle/fixtures/meuse_spde_predict.jld2)

here <- tryCatch(
    normalizePath(dirname(sys.frame(1)$ofile), mustWork = FALSE),
    error = function(e) getwd()
)
if (!nzchar(here)) here <- getwd()
source(file.path(here, "..", "_helpers.R"))

suppressPackageStartupMessages({
    library(INLA)
    library(fmesher)
    library(sp)
})

set.seed(20260424)

data(meuse)
coords <- as.matrix(meuse[, c("x", "y")]) / 1000.0
y <- log(meuse$zinc)
dist_cov <- meuse$dist

mesh <- fmesher::fm_mesh_2d_inla(
    loc       = coords,
    max.edge  = c(0.2, 0.5),
    offset    = c(0.3, 1.0),
    cutoff    = 0.05,
    min.angle = 25
)

spde <- inla.spde2.pcmatern(
    mesh = mesh, alpha = 2,
    prior.range = c(0.5, 0.5),
    prior.sigma = c(1, 0.5)
)

A_obs <- inla.spde.make.A(mesh = mesh, loc = coords)

stk <- inla.stack(
    data = list(y = y),
    A = list(A_obs, 1),
    effects = list(
        field = seq_len(spde$n.spde),
        list(intercept = rep(1, nrow(coords)), dist = dist_cov)
    ),
    tag = "est"
)

form <- y ~ 0 + intercept + dist + f(field, model = spde)

fit <- INLA::inla(
    form, family = "gaussian",
    data = inla.stack.data(stk),
    control.predictor = list(A = inla.stack.A(stk), compute = FALSE),
    control.fixed = list(prec.intercept = 1e-3, prec = 1e-3),
    control.family = list(hyper = list(prec = list(
        prior = "pc.prec", param = c(1, 0.01)
    ))),
    control.compute = list(return.marginals = FALSE)
)

# ----- Vertex-level posterior summary on the SPDE field ------------
# `summary.random$field` is a frame with one row per mesh vertex:
# (ID, mean, sd, 0.025quant, 0.5quant, 0.975quant, mode, kld).
field_summary <- fit$summary.random$field
stopifnot(nrow(field_summary) == mesh$n)
vertex_mean <- as.numeric(field_summary$mean)
vertex_sd   <- as.numeric(field_summary$sd)

# ----- Regular prediction grid over the Meuse hull -----------------
# 0.05 km = 50 m resolution matches the vignette template; covers
# the mesh-interior region that the user ultimately renders. Cells
# whose centre falls outside the mesh get NA in the projected output.
grid_x_centres <- seq(178.55, 181.45, by = 0.05)
grid_y_centres <- seq(329.55, 333.65, by = 0.05)
nx <- length(grid_x_centres)
ny <- length(grid_y_centres)

grid_centres <- cbind(
    rep(grid_x_centres, times = ny),
    rep(grid_y_centres, each  = nx)
)

# ----- R-INLA mesh→pixel projector ---------------------------------
# `inla.mesh.project` was deprecated in INLA 23.06.07 in favour of
# `fmesher::fm_evaluator`; the underlying barycentric formula is the
# same. The returned `proj$A` is the sparse projector with one row
# per grid point; `proj$ok` flags points falling inside the mesh.
ev <- fmesher::fm_evaluator(mesh, loc = grid_centres)
A_predict <- ev$proj$A
ok        <- as.logical(ev$proj$ok)
stopifnot(nrow(A_predict) == nrow(grid_centres))
stopifnot(length(ok)      == nrow(grid_centres))

pixel_mean_vec <- as.numeric(A_predict %*% vertex_mean)
pixel_sd_vec   <- as.numeric(A_predict %*% vertex_sd)
pixel_mean_vec[!ok] <- NA_real_
pixel_sd_vec[!ok]   <- NA_real_

# Reshape into nx × ny matrices indexed by (ix, iy) — column-major
# matches the (length(xs), length(ys)) layout that `Rasters.Raster`
# and Julia store internally.
pixel_mean <- matrix(pixel_mean_vec, nrow = nx, ncol = ny)
pixel_sd   <- matrix(pixel_sd_vec,   nrow = nx, ncol = ny)
ok_mat     <- matrix(ok,             nrow = nx, ncol = ny)

# ----- Emit fixture -------------------------------------------------
out_path <- file.path(here, "..", "fixtures", "rasters",
    "meuse_spde_predict.json")

# Reuse `write_inla_fixture` for the fit-side scaffolding so the
# downstream oracle can also exercise the (model, res) overload.
write_inla_fixture(
    fit = fit,
    path = out_path,
    name = "meuse_spde_predict",
    component_names = character(0),
    include_marginals = FALSE,
    meta = list(
        dataset = "sp::meuse",
        n = nrow(coords),
        family = "gaussian",
        formula = "log(zinc) ~ intercept + dist + f(field, model=spde)",
        coordinate_unit = "km",
        prior_range = c(0.5, 0.5),
        prior_sigma = c(1, 0.5),
        prec_prior = "pc.prec(1, 0.01)",
        fixed_prec = 1e-3,
        grid_resolution_km = 0.05,
        grid_dim = c(nx, ny)
    )
)

# Splice in the input data, the fmesher mesh, the observation
# projector, and the new prediction-grid block. JSON's `null` round-
# trips to Julia `nothing`, which the JLD2 converter then collapses
# to `NaN` on the matrix path.
fixture <- jsonlite::fromJSON(out_path, simplifyVector = FALSE)
fixture$input <- list(
    y         = as.numeric(y),
    dist      = as.numeric(dist_cov),
    locations = lapply(seq_len(nrow(coords)), function(i) as.numeric(coords[i, ]))
)
fixture$mesh <- mesh_to_list(mesh)
fixture$A_field <- sparse_to_triplet(A_obs)
fixture$predict <- list(
    grid_x      = as.numeric(grid_x_centres),
    grid_y      = as.numeric(grid_y_centres),
    nx          = as.integer(nx),
    ny          = as.integer(ny),
    A           = sparse_to_triplet(A_predict),
    ok          = lapply(seq_len(nrow(ok_mat)), function(i) as.logical(ok_mat[i, ])),
    vertex_mean = as.numeric(vertex_mean),
    vertex_sd   = as.numeric(vertex_sd),
    pixel_mean  = lapply(seq_len(nrow(pixel_mean)), function(i) as.numeric(pixel_mean[i, ])),
    pixel_sd    = lapply(seq_len(nrow(pixel_sd)),   function(i) as.numeric(pixel_sd[i, ]))
)

jsonlite::write_json(
    fixture, out_path,
    auto_unbox = TRUE, digits = 16, pretty = FALSE, na = "null"
)

cat("wrote", out_path, "\n")

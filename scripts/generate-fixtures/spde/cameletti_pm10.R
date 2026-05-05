# Synthetic Cameletti-style separable space-time SPDE — R-INLA reference
# fit for `KroneckerComponent(SPDE2, AR1{fix_τ=true})`. Third of three
# Phase M oracle gates.
#
# Cameletti et al. (2013) study daily PM₁₀ on a regional monitoring
# network in Piemonte, Italy, with the separable model
#   y_{s,t} = β + u(s, t) + ε_{s,t},
#   u(s, t) ~ SPDE-Matérn (space) ⊗ AR(1) (time),
#   f(field, model = spde, group = field.group,
#     control.group = list(model = "ar1")).
# Internal hyperparameter layout (R-INLA SPDE2-PCmatern):
#   θ = (τ_obs, range_spde, σ_spde, ρ_AR1).
# The AR1 group has *implicit* prec = 1 for identifiability against
# the spatial precision — Julia matches via `AR1(...; fix_τ=true)`.
#
# This fixture uses a small synthetic dataset (8 stations × 12 time
# points) on a unit square, so the R-INLA fit and the Julia oracle
# both terminate inside CI budgets while still exercising the full
# SPDE2 ⊗ AR1 surface.
#
# Output: fixtures/spde/cameletti_pm10.json
# (→ packages/INLASPDE.jl/test/oracle/fixtures/cameletti_pm10.jld2)

here <- tryCatch(
    normalizePath(dirname(sys.frame(1)$ofile), mustWork = FALSE),
    error = function(e) getwd()
)
if (!nzchar(here)) here <- getwd()
source(file.path(here, "..", "_helpers.R"))

suppressPackageStartupMessages({
    library(INLA)
    library(fmesher)
    library(Matrix)
})

set.seed(20260506)

# --- domain + station/time grid --------------------------------------
n_s   <- 8                                       # stations
n_t   <- 12                                      # time points
loc   <- cbind(runif(n_s, 0, 1), runif(n_s, 0, 1))
times <- seq_len(n_t)

# Long-format: every (station, time) pair generates one observation.
station_id <- rep(seq_len(n_s), times = n_t)
time_id    <- rep(seq_len(n_t), each = n_s)
loc_long   <- loc[station_id, , drop = FALSE]
n_obs      <- n_s * n_t

# --- mesh -------------------------------------------------------------
mesh <- fmesher::fm_mesh_2d_inla(
    loc      = loc,
    max.edge = c(0.15, 0.4),
    offset   = c(0.1, 0.3),
    cutoff   = 0.05
)
n_v <- mesh$n

# --- SPDE2-PCmatern ---------------------------------------------------
# PC priors mirror `synthetic_spde_1d.R` / `meuse_spde.R`:
#   range ~ PC(0.5, 0.5),  σ ~ PC(1, 0.5).
spde <- inla.spde2.pcmatern(
    mesh        = mesh, alpha = 2,
    prior.range = c(0.5, 0.5),
    prior.sigma = c(1.0, 0.5)
)

# --- ground-truth simulation -----------------------------------------
# Internal-scale "truth" used only to *simulate* y. Recovery target is
# the R-INLA posterior, which is what the oracle test compares against.
range_true <- 0.5
sigma_true <- 1.0
rho_true   <- 0.7
beta_true  <- 1.0
tau_obs_true <- 4.0    # noise precision (σ_y = 0.5)

# Simulate u(v, t) at vertex × time, by drawing each AR1 path
# separately at every vertex and then projecting per-time-slice.
# inla.spde2.precision returns Q at user-scale (range, σ).
Q_spde <- inla.spde2.precision(spde, theta = c(log(range_true), log(sigma_true)))
# Cholesky factor for sampling: L L' = Q⁻¹.
# Use chol(Q) and then solve(L', z) for one independent draw per slice.
L_spde <- chol(Q_spde)            # upper-triangular: Q = Lᵀ L
indep <- matrix(rnorm(n_v * n_t), n_v, n_t)
u_indep <- backsolve(L_spde, indep)   # solves Lᵀ x = z ⇒ x ~ N(0, Q⁻¹)

# Apply AR1 in the time direction: u_t = ρ u_{t-1} + sqrt(1 - ρ²) ν_t.
u_field <- matrix(0.0, n_v, n_t)
u_field[, 1] <- u_indep[, 1]
for (t in 2:n_t) {
    u_field[, t] <- rho_true * u_field[, t - 1] +
                    sqrt(1 - rho_true^2) * u_indep[, t]
}

# Observation projector: per-time-slice spatial A. The full observation
# matrix maps `c(u_field)` (length n_v · n_t, time-major) to long-format
# (s, t)-indexed observations. We store the time-major A_field below.
A_space <- inla.spde.make.A(mesh = mesh, loc = loc)        # n_s × n_v
# Build A_field: row k = (s_k, t_k); column index = (t_k - 1) · n_v + v.
# Time-major flattening: c(u_field) lays out u[, 1], u[, 2], ..., u[, n_t].
A_full_rows <- vector("list", n_obs)
for (k in seq_len(n_obs)) {
    s <- station_id[k]
    tt <- time_id[k]
    row_in_space <- A_space[s, , drop = FALSE]
    A_full_rows[[k]] <- cbind(
        i = rep(k, length(row_in_space@x)),
        j = (tt - 1) * n_v + (row_in_space@j + 1L),
        v = row_in_space@x
    )
}
trip <- do.call(rbind, A_full_rows)
A_field <- sparseMatrix(i = trip[, 1], j = trip[, 2], x = trip[, 3],
                        dims = c(n_obs, n_v * n_t))

# Mean signal at each observation = (A_field) %*% c(u_field).
u_vec <- as.numeric(u_field)         # time-major (column-stacked)
mu    <- as.numeric(A_field %*% u_vec)
y     <- beta_true + mu + rnorm(n_obs, sd = 1 / sqrt(tau_obs_true))

# --- fit (separable space-time via control.group=list(model="ar1")) --
# inla.spde.make.A with `group` builds the same time-blocked projector
# we just constructed by hand. Feed it directly.
A_inla <- inla.spde.make.A(
    mesh  = mesh,
    loc   = loc_long,
    group = time_id,
    n.group = n_t
)
# Sanity check: A_inla and A_field should agree (column convention is
# the same: time-major (vertex, time) → c(u_field) layout).
stopifnot(all(dim(A_inla) == dim(A_field)))

field_index <- inla.spde.make.index(
    name      = "field",
    n.spde    = spde$n.spde,
    n.group   = n_t
)

stk <- inla.stack(
    data = list(y = y),
    A = list(A_inla, 1),
    effects = list(
        c(field_index,
          list(intercept = rep(1, n_obs))),
        list()
    ),
    tag = "est"
)
# (intercept is part of `effects` block 1 above so the data layout is
# consistent with `field_index`. The formula references `intercept`
# explicitly.)

# Need a *separate* intercept block, so re-build the stack with two
# effects blocks the standard Cameletti way:
stk <- inla.stack(
    data = list(y = y),
    A = list(A_inla, 1),
    effects = list(
        c(field_index),
        list(intercept = rep(1, n_obs))
    ),
    tag = "est"
)

form <- y ~ 0 + intercept +
    f(field, model = spde,
      group = field.group, control.group = list(model = "ar1"))

fit <- INLA::inla(
    form, family = "gaussian",
    data = inla.stack.data(stk),
    control.predictor = list(A = inla.stack.A(stk), compute = FALSE),
    control.fixed = list(prec.intercept = 1.0e-3, prec = 1.0e-3),
    control.family = list(hyper = list(prec = list(
        prior = "pc.prec", param = c(1, 0.01)
    ))),
    control.compute = list(return.marginals = FALSE)
)

# --- write fixture ---------------------------------------------------
out_path <- file.path(here, "..", "fixtures", "spde", "cameletti_pm10.json")

write_inla_fixture(
    fit = fit,
    path = out_path,
    name = "cameletti_pm10",
    component_names = character(0),
    include_marginals = FALSE,
    meta = list(
        dataset = "synthetic_cameletti_pm10",
        n_obs = n_obs,
        n_stations = n_s,
        n_time = n_t,
        n_vertices = n_v,
        family = "gaussian",
        formula = "y ~ intercept + f(field, model=spde, group, control.group=ar1)",
        true = list(
            range = range_true, sigma = sigma_true, rho = rho_true,
            beta = beta_true, tau_obs = tau_obs_true
        ),
        prior_range = c(0.5, 0.5),
        prior_sigma = c(1.0, 0.5),
        prec_prior = "pc.prec(1, 0.01)",
        fixed_prec = 1.0e-3
    )
)

# Append the data + projector + station/time index for the Julia rebuild.
fixture <- jsonlite::fromJSON(out_path, simplifyVector = FALSE)
fixture$input <- list(
    y          = as.numeric(y),
    locations  = lapply(seq_len(nrow(loc)), function(i) as.numeric(loc[i, ])),
    station_id = as.integer(station_id),
    time_id    = as.integer(time_id),
    times      = as.integer(times)
)
fixture$mesh    <- mesh_to_list(mesh)
# Ship the *spatial* A (n_obs × n_v): one row per observation, column
# = vertex. Julia builds `KroneckerMapping(A_space, A_time)` with
# `A_time` = identity-row-selector at `time_id`. We export both for
# transparency.
fixture$A_field_full <- sparse_to_triplet(A_inla)
A_per_obs_space <- A_space[station_id, , drop = FALSE]
fixture$A_space <- sparse_to_triplet(A_per_obs_space)

jsonlite::write_json(
    fixture, out_path,
    auto_unbox = TRUE, digits = 16, pretty = FALSE, na = "null"
)

cat("wrote", out_path, "\n")

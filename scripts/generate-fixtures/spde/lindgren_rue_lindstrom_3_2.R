# Reference R-INLA non-stationary SPDE fit for Lindgren-Rue-Lindström
# (2011, JRSSB §3.2). Second of three Phase M oracle gates.
#
# Domain: unit square [0, 1]^2.
# Truth:
#   y_i  = β + u(s_i) + ε_i,   i = 1,…,n_obs,  s_i ∈ [0, 1]^2
#   u    ~ non-stationary SPDE-Matérn 2D, α = 2 (ν = 1)
#   log τ_v = θ_τ                                     (constant)
#   log κ_v = θ_κ_1 · 1[x_v < 0.5] + θ_κ_2 · 1[x_v ≥ 0.5]   (two regions)
#   ε_i ~ N(0, 1/τ_obs)
#   β ~ N(0, 1000)
#
# Truth on internal scale: θ_τ = 0, θ_κ_1 = 1, θ_κ_2 = 2.
#   ⇒ τ = 1 everywhere, κ = e ≈ 2.72 in region 1, κ = e² ≈ 7.39 in region 2.
#
# R-INLA B-matrix convention: column 1 is a fixed offset (×1), columns
# 2..(p+1) are scaled by θ. Julia's `SPDE2NonStationary` scales every
# column by its own θ, so the Julia oracle test drops R-INLA's fixed-
# offset column on read-back.
#
# The fixture ships the mesh, projector A, B_τ / B_κ matrices and the
# observations so the Julia side rebuilds the LGM directly from the
# fmesher mesh — isolates SPDE2NonStationary + INLA from any 2D mesh-
# parity gap.
#
# Output: fixtures/spde/lindgren_rue_lindstrom_3_2.json
# (→ packages/INLASPDE.jl/test/oracle/fixtures/lindgren_rue_lindstrom_3_2.jld2)

here <- tryCatch(
    normalizePath(dirname(sys.frame(1)$ofile), mustWork = FALSE),
    error = function(e) getwd()
)
if (!nzchar(here)) here <- getwd()
source(file.path(here, "..", "_helpers.R"))

suppressPackageStartupMessages({
    library(INLA)
    library(fmesher)
})

set.seed(20260505)

# --- domain + observations --------------------------------------------
n_obs <- 200
loc <- cbind(runif(n_obs, 0, 1), runif(n_obs, 0, 1))

# --- mesh -------------------------------------------------------------
mesh <- fmesher::fm_mesh_2d_inla(
    loc      = loc,
    max.edge = c(0.05, 0.2),
    offset   = c(0.1, 0.3),
    cutoff   = 0.02
)
n_v <- mesh$n

# --- non-stationary basis matrices ------------------------------------
# Region indicator at vertex level.
region1 <- as.numeric(mesh$loc[, 1] < 0.5)
region2 <- 1 - region1

# R-INLA convention: B-matrix has shape (n_v, 1 + n_theta) where n_theta
# is shared across τ and κ. Column 1 is a fixed offset (multiplied by 1
# always); columns 2..(n_theta+1) are scaled by θ_1, …, θ_{n_theta}.
# We pick the global ordering θ = (θ_τ_1, θ_κ_1, θ_κ_2):
#   log τ_v = 0·1 + 1·θ_τ_1 + 0·θ_κ_1 + 0·θ_κ_2 = θ_τ_1
#   log κ_v = 0·1 + 0·θ_τ_1 + region1·θ_κ_1 + region2·θ_κ_2
B.tau <- cbind(rep(0, n_v), rep(1, n_v), rep(0, n_v), rep(0, n_v))
B.kappa <- cbind(rep(0, n_v), rep(0, n_v), region1, region2)

# --- ground-truth simulation -----------------------------------------
true_theta <- c(0.0, 1.0, 2.0)         # θ_τ_1, θ_κ_1, θ_κ_2

spde_for_sim <- inla.spde2.matern(
    mesh, alpha = 2,
    B.tau = B.tau, B.kappa = B.kappa,
    theta.prior.mean = c(0.0, 0.0, 0.0),
    theta.prior.prec = c(1.0, 1.0, 1.0)
)
Q_true <- inla.spde2.precision(spde_for_sim, theta = true_theta)
u_true <- as.numeric(inla.qsample(n = 1, Q = Q_true, seed = 20260505))

A <- inla.spde.make.A(mesh = mesh, loc = loc)
mu <- as.numeric(A %*% u_true)
beta_true <- 1.5
sigma_y <- 0.5
y <- beta_true + mu + rnorm(n_obs, sd = sigma_y)

# --- fit --------------------------------------------------------------
spde <- inla.spde2.matern(
    mesh, alpha = 2,
    B.tau = B.tau, B.kappa = B.kappa,
    theta.prior.mean = c(0.0, 0.0, 0.0),
    theta.prior.prec = c(1.0, 1.0, 1.0)
)

stk <- inla.stack(
    data = list(y = y),
    A = list(A, 1),
    effects = list(
        field = seq_len(spde$n.spde),
        list(intercept = rep(1, n_obs))
    ),
    tag = "est"
)

form <- y ~ 0 + intercept + f(field, model = spde)

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

# --- write fixture ---------------------------------------------------
out_path <- file.path(here, "..", "fixtures", "spde",
                     "lindgren_rue_lindstrom_3_2.json")

write_inla_fixture(
    fit = fit,
    path = out_path,
    name = "lindgren_rue_lindstrom_3_2",
    component_names = character(0),
    include_marginals = FALSE,
    meta = list(
        dataset = "synthetic_lrl_3_2",
        n_obs = n_obs,
        n_vertices = n_v,
        family = "gaussian",
        formula = "y ~ intercept + f(field, model = spde2.matern{B.tau,B.kappa})",
        true_theta = true_theta,
        beta_true = beta_true,
        sigma_y_true = sigma_y,
        prior_theta_mean = c(0.0, 0.0, 0.0),
        prior_theta_prec = c(1.0, 1.0, 1.0),
        prec_prior = "pc.prec(1, 0.01)",
        fixed_prec = 1e-3
    )
)

# Append input + mesh + A + non-stationary basis matrices for the
# Julia-side oracle test. Julia drops R-INLA's fixed-offset column 1
# on read-back, so we ship the trimmed (θ-scaled-only) basis here too.
fixture <- jsonlite::fromJSON(out_path, simplifyVector = FALSE)
fixture$input <- list(
    y         = as.numeric(y),
    locations = lapply(seq_len(nrow(loc)), function(i) as.numeric(loc[i, ]))
)
fixture$mesh <- mesh_to_list(mesh)
fixture$A_field <- sparse_to_triplet(A)
# Julia-friendly basis: Julia uses separate B_τ (n_v × p_τ) and
# B_κ (n_v × p_κ). Drop R-INLA's fixed-offset column 1, then for each
# matrix retain only the columns whose corresponding θ is "owned" by
# that side: τ owns θ_1; κ owns θ_2, θ_3.
# Wrap each row in I() (AsIs) so jsonlite's auto_unbox=TRUE does not
# collapse length-1 numeric vectors (B_tau row) to scalars.
fixture$B_tau <- lapply(seq_len(n_v),
    function(i) I(as.numeric(B.tau[i, 2, drop = FALSE])))
fixture$B_kappa <- lapply(seq_len(n_v),
    function(i) I(as.numeric(B.kappa[i, 3:4, drop = FALSE])))

jsonlite::write_json(
    fixture, out_path,
    auto_unbox = TRUE, digits = 16, pretty = FALSE, na = "null"
)

cat("wrote", out_path, "\n")

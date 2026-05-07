## R-INLA benchmark harness for the same three oracle fixtures as the
## Julia side. Reads the JLD2 fixtures via the python-style HDF5 backend
## (rhdf5) is overkill — we re-derive the inputs from the original raw
## datasets that R-INLA ships with, so both sides start from identical
## inputs without an HDF5 round-trip.
##
## Three fits, timed end-to-end with `system.time()` median over `n_runs`
## (default 5). First run discarded as warm-up. Memory peak via
## `gc(reset=TRUE)` + `gc()$used` delta in MiB.
##
## Pinned: R-INLA "stable" channel; version reported in the results
## header by `run.jl`.

suppressPackageStartupMessages({
    library(INLA)
    library(SpatialEpi)   # Pennsylvania, Scotland lip cancer
    library(spdep)
    library(sp)
})

inla.setOption("num.threads", "1:1")  # 1 outer x 1 inner thread
inla.setOption(fmesher.evolution.warn = FALSE)

INLA_VERSION <- as.character(packageVersion("INLA"))
R_VERSION    <- paste(R.version$major, R.version$minor, sep = ".")

## --- Scotland BYM2 ---------------------------------------------------------

scotland_data <- function() {
    data(scotland)
    Y <- scotland$data$cases
    E <- scotland$data$expected
    x <- scotland$data$AFF
    nb <- poly2nb(scotland$spatial.polygon)
    nb2INLA("/tmp/inla_bench_scotland.adj", nb)
    list(
        df    = data.frame(idx = seq_along(Y), y = Y, E = E, x = x),
        graph = "/tmp/inla_bench_scotland.adj"
    )
}

scotland_fit <- function(d) {
    formula <- y ~ x + f(idx, model = "bym2", graph = d$graph,
                         scale.model = TRUE,
                         hyper = list(
                             prec = list(prior = "pc.prec",
                                         param = c(1, 0.01)),
                             phi  = list(prior = "pc",
                                         param = c(0.5, 0.5))
                         ))
    inla(formula, family = "poisson", data = d$df, E = d$df$E,
         control.predictor = list(compute = FALSE),
         control.compute   = list(dic = FALSE, waic = FALSE, cpo = FALSE),
         num.threads       = "1:1")
}

## --- Pennsylvania BYM2 -----------------------------------------------------

pennsylvania_data <- function() {
    data(pennLC)
    df_county <- aggregate(cbind(cases, population) ~ county,
                           data = pennLC$data, FUN = sum)
    df_county$cases      <- as.numeric(df_county$cases)
    df_county$population <- as.numeric(df_county$population)
    smk <- pennLC$smoking
    df  <- merge(df_county, smk, by = "county")
    df$expected <- df$population * sum(df$cases) / sum(df$population)
    poly <- pennLC$spatial.polygon
    poly_order <- sapply(poly@polygons, function(p) p@ID)
    df  <- df[match(poly_order, df$county), ]
    nb <- poly2nb(poly)
    nb2INLA("/tmp/inla_bench_pa.adj", nb)
    df$idx <- seq_len(nrow(df))
    list(df = df, graph = "/tmp/inla_bench_pa.adj")
}

pennsylvania_fit <- function(d) {
    formula <- cases ~ smoking +
        f(idx, model = "bym2", graph = d$graph,
          scale.model = TRUE,
          hyper = list(
              prec = list(prior = "pc.prec", param = c(1, 0.01)),
              phi  = list(prior = "pc",       param = c(0.5, 0.5))
          ))
    inla(formula, family = "poisson", data = d$df, E = d$df$expected,
         control.predictor = list(compute = FALSE),
         control.compute   = list(dic = FALSE, waic = FALSE, cpo = FALSE),
         num.threads       = "1:1")
}

## --- Meuse SPDE ------------------------------------------------------------

meuse_data <- function() {
    data(meuse)
    coordinates(meuse) <- ~x + y
    loc  <- coordinates(meuse) / 1000  # km
    mesh <- inla.mesh.2d(loc, max.edge = c(0.3, 1.5),
                         cutoff = 0.05, offset = c(0.4, 1.5))
    spde <- inla.spde2.pcmatern(
        mesh = mesh, alpha = 2,
        prior.range = c(0.5, 0.5),
        prior.sigma = c(1.0, 0.5)
    )
    A <- inla.spde.make.A(mesh = mesh, loc = loc)
    stk <- inla.stack(
        data    = list(y = log(meuse$zinc)),
        A       = list(A, 1),
        effects = list(field = seq_len(spde$n.spde),
                       data.frame(intercept = 1,
                                  dist      = meuse$dist)),
        tag = "est"
    )
    list(stk = stk, spde = spde)
}

meuse_fit <- function(d) {
    formula <- y ~ -1 + intercept + dist + f(field, model = d$spde)
    inla(formula, family = "gaussian",
         data              = inla.stack.data(d$stk),
         control.predictor = list(A = inla.stack.A(d$stk), compute = FALSE),
         control.compute   = list(dic = FALSE, waic = FALSE, cpo = FALSE),
         control.family    = list(
             hyper = list(prec = list(prior = "pc.prec",
                                      param = c(1, 0.01)))
         ),
         num.threads = "1:1")
}

## --- Driver ----------------------------------------------------------------

time_one <- function(name, data_fn, fit_fn, n_runs = 5) {
    d <- data_fn()
    invisible(fit_fn(d))  # warm-up

    times <- numeric(n_runs)
    mlik  <- NA_real_
    mem_used_start <- sum(gc(reset = TRUE)[, 2])
    for (i in seq_len(n_runs)) {
        gc()
        t0  <- Sys.time()
        fit <- fit_fn(d)
        t1  <- Sys.time()
        times[i] <- as.numeric(difftime(t1, t0, units = "secs"))
        if (i == n_runs) mlik <- fit$mlik[1, 1]
    }
    mem_used_end <- sum(gc()[, 2])

    list(
        name           = name,
        median_seconds = median(times),
        min_seconds    = min(times),
        max_seconds    = max(times),
        iqr_seconds    = unname(diff(quantile(times, c(0.25, 0.75)))),
        mem_peak_MiB   = mem_used_end - mem_used_start,
        n_samples      = n_runs,
        mlik           = mlik
    )
}

run_all_R <- function(n_runs = 5) {
    list(
        scotland_bym2     = time_one("scotland_bym2",
                                     scotland_data,     scotland_fit,     n_runs),
        pennsylvania_bym2 = time_one("pennsylvania_bym2",
                                     pennsylvania_data, pennsylvania_fit, n_runs),
        meuse_spde        = time_one("meuse_spde",
                                     meuse_data,        meuse_fit,        n_runs)
    )
}

using Test

# Tier-3 triangulation gate (Phase P PR-1 / ADR-044).
#
# Scotland and Pennsylvania BYM2 NUTS chains run on every CI pass —
# each is ~30-60 s on a workstation. The Meuse SPDE chain runs only
# when explicitly opted in via `--triangulation` (or
# `LGMTURING_TRIANGULATION=1`), because each NUTS leapfrog step costs
# a full Laplace fit of a 355-dim latent and the chain takes ~10 min.
# The nightly CI workflow flips this on; PR builds skip it.
#
# `INLASPDE` is intentionally *not* in `[extras]` because the personal
# registry caps it at 0.1.x while the in-tree source is 0.3.0; Pkg.test's
# sandbox would re-resolve to the stale registry version. The nightly
# workflow runs `Pkg.develop` on `INLASPDE` before invoking
# `Pkg.test(test_args=["--triangulation"])`, so the dev-link is
# discoverable via `Base.find_package` at that point.
const RUN_FULL_TRIANGULATION =
    ("--triangulation" in ARGS) ||
    get(ENV, "LGMTURING_TRIANGULATION", "0") == "1"

@testset "LGMTuring.jl" begin
    @testset "regression" begin
        include("regression/test_logdensity.jl")
        include("regression/test_nuts_sample.jl")
        include("regression/test_compare.jl")
    end
    @testset "triangulation" begin
        include("triangulation/test_scotland_bym2.jl")
        include("triangulation/test_pennsylvania_bym2.jl")
        if RUN_FULL_TRIANGULATION
            if Base.find_package("INLASPDE") === nothing
                @warn "INLASPDE not installed in the test environment — skipping " *
                      "Meuse SPDE triangulation. Run `Pkg.develop(path=\"" *
                      "packages/INLASPDE.jl\")` in the test env first."
            else
                include("triangulation/test_meuse_spde.jl")
            end
        else
            @info "Skipping Meuse SPDE triangulation (pass `--triangulation` " *
                  "or set LGMTURING_TRIANGULATION=1 to enable; ~10 min)."
        end
    end
end

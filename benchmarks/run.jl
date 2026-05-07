#!/usr/bin/env julia

"""
Performance benchmark entry point. Run as:

    julia --project=benchmarks benchmarks/run.jl

Pins single-threaded BLAS on both sides, fits Scotland BYM2, Pennsylvania
BYM2, and Meuse SPDE on each engine, writes the comparison to
`benchmarks/results/<date>_<arch>.md`.

Methodology, version pins, hardware spec all land in the results file
header. The `[sources]` table in `benchmarks/Project.toml` points at
the in-tree packages; `Pkg.instantiate()` resolves them on first run.
"""

using Pkg
Pkg.activate(@__DIR__)

using LinearAlgebra
BLAS.set_num_threads(1)

using Dates
using Printf
using RCall

include(joinpath(@__DIR__, "harness.jl"))
include(joinpath(@__DIR__, "compare.jl"))

const N_SAMPLES = 5
const SECONDS = 120.0  # generous per-dataset budget; only Meuse SPDE approaches it

function _hardware_spec()
    if Sys.isapple()
        cpu = strip(read(`sysctl -n machdep.cpu.brand_string`, String))
        cores_phys = strip(read(`sysctl -n hw.physicalcpu`, String))
        cores_log = strip(read(`sysctl -n hw.logicalcpu`, String))
        ram_bytes = parse(Int, strip(read(`sysctl -n hw.memsize`, String)))
        ram_gib = round(ram_bytes / 2^30; digits=1)
        return "$(cpu) — $(cores_phys) physical / $(cores_log) logical cores, $(ram_gib) GiB RAM"
    else
        return "$(Sys.cpu_info()[1].model) — $(Sys.CPU_THREADS) logical threads"
    end
end

function _os_string()
    if Sys.isapple()
        prod = strip(read(`sw_vers -productName`, String))
        ver = strip(read(`sw_vers -productVersion`, String))
        return "$(prod) $(ver) ($(Sys.MACHINE))"
    else
        return "$(Sys.KERNEL) $(Sys.MACHINE)"
    end
end

function _r_versions()
    R"version_str <- paste(R.version$major, R.version$minor, sep = '.')"
    R"inla_version <- as.character(packageVersion('INLA'))"
    return (
        r_version=String(rcopy(R"version_str")),
        inla_version=String(rcopy(R"inla_version"))
    )
end

function _julia_pkg_versions()
    deps = Pkg.dependencies()
    pkgs = (:GMRFs, :LatentGaussianModels, :INLASPDE)
    out = NamedTuple()
    for p in pkgs
        for (_, info) in deps
            if info.name === String(p) && info.version !== nothing
                out = merge(out, NamedTuple{(p,)}((string(info.version),)))
                break
            end
        end
    end
    return out
end

function main()
    println("== INLA.jl performance benchmark ==")
    println("BLAS threads: ", BLAS.get_num_threads())
    println()

    println("[1/3] Pinning R-INLA to 1 thread + sourcing benchmarks/r_inla.R …")
    r_path = replace(joinpath(@__DIR__, "r_inla.R"), "\\" => "/")
    R"source($r_path)"

    println("[2/3] Running R-INLA timed fits ($N_SAMPLES samples per dataset, warm-up discarded) …")
    R"r_results <- run_all_R(n_runs = $N_SAMPLES)"
    r_results = rcopy(R"r_results")  # Dict{Symbol,Dict{Symbol,Any}}
    for k in (:scotland_bym2, :pennsylvania_bym2, :meuse_spde)
        @printf("   %-18s median %.2fs (n=%d)\n", String(k),
            r_results[k][:median_seconds], r_results[k][:n_samples])
    end

    println("[3/3] Running INLA.jl timed fits …")
    julia_rows = run_all_julia(samples=N_SAMPLES, seconds=SECONDS)
    for jr in julia_rows
        @printf("   %-18s median %.2fs (n=%d)\n", String(jr.name),
            jr.median_seconds, jr.n_samples)
    end

    r_rows = (r_results[jr.name] for jr in julia_rows) |> collect

    rv = _r_versions()
    header_meta = (
        date=Dates.format(Dates.today(), "yyyy-mm-dd"),
        hardware=_hardware_spec(),
        os=_os_string(),
        blas_threads=BLAS.get_num_threads(),
        julia_version=string(VERSION),
        r_version=rv.r_version,
        inla_version=rv.inla_version,
        julia_pkgs=_julia_pkg_versions(),
        n_samples=N_SAMPLES
    )

    arch = Sys.isapple() ? "apple_$(Sys.ARCH)" : String(Sys.ARCH)
    results_dir = joinpath(@__DIR__, "results")
    isdir(results_dir) || mkpath(results_dir)
    outpath = joinpath(results_dir, "$(header_meta.date)_$(arch).md")

    write_results(julia_rows, r_rows, header_meta; outpath=outpath)
    println()
    println("Wrote ", relpath(outpath, dirname(@__DIR__)))
    return outpath
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

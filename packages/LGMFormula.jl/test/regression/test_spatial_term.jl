# Phase N PR-7b / PR-7c — tuple-coordinate `f(...)` syntax for SPDE
# spatial and separable space-time terms (ADR-037, ADR-038, ADR-039).
# The parser accepts a tuple of column symbols as the first argument of
# `f(...)`; length 2 (spatial) and length 3 (space-time) are accepted,
# others rejected. At runtime the term dispatches through
# `_build_spatial_block`, overloaded for `SPDE2` (2-tuple) and
# `KroneckerComponent` (3-tuple) in `LGMFormulaINLASPDEExt` (loaded
# automatically when both `LGMFormula` and `INLASPDE` are imported).
#
# Tests cover: parse-time errors (tuple arity, non-Symbol entries,
# `replicate`/`group` rejection), macroexpand structure (no
# `_wrap_term` wrap on tuple-coord terms), the missing-extension error
# path, the SPDE2-via-mesh roundtrip, the no-mesh error from the
# `(points, triangles)` `SPDE2` constructor (ADR-036), the
# Cameletti-shape gridded SPDE2 ⊗ AR1 roundtrip, and the per-obs
# (non-gridded) Khatri-Rao construction (PR-7c).

using INLASPDE: SPDE2, MeshProjector, INLAMesh, inla_mesh_2d,
                PCMatern, num_vertices
using LinearAlgebra: I
using SparseArrays: findnz

@testset "PR-7b parse-time errors for tuple-coord f-term" begin
    df = (y=Float64[], east=Float64[], north=Float64[],
        time=Float64[], idx=Int[], id=Int[], grp=Int[])

    @testset "tuple of length 1 rejected" begin
        ex = try
            @eval @lgm y~1 + f((east,), IID(5)) data=$df family=GaussianLikelihood()
            nothing
        catch e
            e
        end
        @test ex isa LoadError || ex isa ErrorException
        msg = sprint(showerror, ex)
        @test occursin("2 entries", msg) || occursin("3 entries", msg)
    end

    @testset "tuple of length 4 rejected" begin
        ex = try
            @eval @lgm y~1 + f((east, north, time, idx), IID(5)) data=$df family=GaussianLikelihood()
            nothing
        catch e
            e
        end
        @test ex isa LoadError || ex isa ErrorException
        msg = sprint(showerror, ex)
        @test occursin("2 entries", msg) || occursin("3 entries", msg) ||
              occursin("length 4", msg)
    end

    @testset "non-Symbol tuple entry rejected" begin
        ex = try
            @eval @lgm y~1 + f((east, 5), IID(5)) data=$df family=GaussianLikelihood()
            nothing
        catch e
            e
        end
        @test ex isa LoadError || ex isa ErrorException
        msg = sprint(showerror, ex)
        @test occursin("Symbol", msg) || occursin("column name", msg)
    end

    @testset "replicate kwarg on tuple-coord f-term rejected" begin
        ex = try
            @eval @lgm y~1 + f((east, north), IID(5); replicate=id) data=$df family=GaussianLikelihood()
            nothing
        catch e
            e
        end
        @test ex isa LoadError || ex isa ErrorException
        msg = sprint(showerror, ex)
        @test occursin("replicate", msg) || occursin("group", msg) ||
              occursin("tuple", msg) || occursin("spatial", msg)
    end

    @testset "group kwarg on tuple-coord f-term rejected" begin
        ex = try
            @eval @lgm y~1 + f((east, north), IID; group=grp) data=$df family=GaussianLikelihood()
            nothing
        catch e
            e
        end
        @test ex isa LoadError || ex isa ErrorException
        msg = sprint(showerror, ex)
        @test occursin("replicate", msg) || occursin("group", msg) ||
              occursin("tuple", msg) || occursin("spatial", msg)
    end
end

@testset "PR-7b @macroexpand structure for tuple-coord term" begin
    df_dummy = (y=Float64[], east=Float64[], north=Float64[])

    # We use a placeholder symbol `spde_dummy` for the component slot —
    # we only care that the AST shape is right. The runtime is not
    # invoked here, so the symbol need not resolve.
    @testset "tuple-coord f-term: comp_expr appears literally (no _wrap_term)" begin
        ex = @macroexpand @lgm y~1 + f((east, north), spde_dummy) data=df_dummy family=GaussianLikelihood()
        call = _find_lgm_call(ex)
        @test call !== nothing
        components = call.args[3]
        @test components.head === :tuple
        @test length(components.args) == 2
        @test _is_call_to(components.args[1], :Intercept)
        # Second component is the literal `spde_dummy`, NOT a `_wrap_term(...)` call.
        @test components.args[2] === :spde_dummy
        @test !_is_call_to(components.args[2], :_wrap_term)
    end

    @testset "tuple-coord f-term: mapping is _build_design_matrix(...)" begin
        ex = @macroexpand @lgm y~1 + f((east, north), spde_dummy) data=df_dummy family=GaussianLikelihood()
        call = _find_lgm_call(ex)
        mapping = call.args[4]
        @test _is_call_to(mapping, :_build_design_matrix)
    end
end

@testset "PR-7b missing-extension error (non-SPDE component)" begin
    rng = Random.Xoshiro(20260506)
    n = 5
    df = (y=randn(rng, n),
        east=rand(rng, n),
        north=rand(rng, n))

    # `IID(5)` is *not* an SPDE component, but the parser accepts the
    # tuple-coord syntax syntactically. At runtime, dispatch falls
    # through to the default `_build_spatial_block`, which raises a
    # user-readable error pointing at the INLASPDE extension.
    @test_throws ArgumentError begin
        @lgm y~1 + f((east, north), IID(5)) data=df family=GaussianLikelihood()
    end
    try
        @lgm y~1 + f((east, north), IID(5)) data=df family=GaussianLikelihood()
    catch e
        msg = sprint(showerror, e)
        @test occursin("INLASPDE", msg) || occursin("SPDE", msg) ||
              occursin("extension", msg)
    end
end

@testset "PR-7b SPDE2 spatial roundtrip via LGMFormulaINLASPDEExt" begin
    rng = Random.Xoshiro(20260506)
    sq = [0.0 0.0; 1.0 0.0; 1.0 1.0; 0.0 1.0]
    mesh = inla_mesh_2d(; boundary=sq, max_edge=0.5, min_angle=25.0)

    # Random observations strictly inside the mesh domain.
    n = 12
    df = (
        logzinc=randn(rng, n),
        dist=randn(rng, n),
        east=0.1 .+ 0.8 .* rand(rng, n),
        north=0.1 .+ 0.8 .* rand(rng, n)
    )
    pc = PCMatern(range_U=0.5, range_α=0.05, sigma_U=1.0, sigma_α=0.01)
    spde = SPDE2(mesh; pc=pc)

    @testset "Meuse-shape: 1 + dist + f((east, north), spde)" begin
        model_macro = @lgm logzinc~1 + dist + f((east, north), spde) data=df family=GaussianLikelihood()

        # Hand-written equivalent: same intercept + covariate blocks,
        # spatial block built via `MeshProjector(mesh, [east north])`.
        A_int = sparse(ones(n, 1))
        A_dist = sparse(reshape(Float64.(df.dist), n, 1))
        P = MeshProjector(mesh, [df.east df.north])
        A_spde = SparseMatrixCSC{Float64, Int}(P.A)
        A_hand = hcat(A_int, A_dist, A_spde)
        model_hand = LatentGaussianModel(
            GaussianLikelihood(),
            (Intercept(), FixedEffects(1), spde),
            A_hand
        )

        @test _struct_isequal(model_macro.likelihoods, model_hand.likelihoods)
        @test _struct_isequal(model_macro.components, model_hand.components)

        # `MeshProjector` is invoked twice independently (once inside
        # the macro expansion, once in the hand form). Per ADR-036,
        # the projector matrices agree to fp roundoff but not exactly:
        # `find_triangle` returns the enclosing triangle in different
        # cyclic vertex orders between independent calls, and
        # `_barycentric` is not symmetric under that rotation. Compare
        # via `A * x ≈ A_hand * x` on a random vector instead of `==`.
        n_v = num_vertices(mesh)
        n_cols = 1 + 1 + n_v
        @test size(model_macro.mapping.A) == (n, n_cols)
        @test size(model_hand.mapping.A) == (n, n_cols)
        x = randn(rng, n_cols)
        @test (model_macro.mapping.A * x)≈(model_hand.mapping.A * x) rtol=1.0e-12
    end

    @testset "no-intercept: f((east, north), spde) only" begin
        model_macro = @lgm logzinc~0 + f((east, north), spde) data=df family=GaussianLikelihood()
        @test length(model_macro.components) == 1
        @test _struct_isequal(model_macro.components[1], spde)
        @test size(model_macro.mapping.A) == (n, num_vertices(mesh))
    end
end

@testset "PR-7b SPDE2 without retained mesh raises a user-readable error" begin
    # The (points, triangles) constructor leaves `spde.mesh = nothing`
    # per ADR-036. The macro path needs the mesh to build a barycentric
    # projector, so the extension raises an error pointing at the
    # mesh-aware constructor.
    rng = Random.Xoshiro(20260506)
    sq = [0.0 0.0; 1.0 0.0; 1.0 1.0; 0.0 1.0]
    mesh = inla_mesh_2d(; boundary=sq, max_edge=0.5, min_angle=25.0)
    spde_no_mesh = SPDE2(mesh.points, mesh.triangles)

    n = 4
    df = (y=randn(rng, n),
        east=0.1 .+ 0.8 .* rand(rng, n),
        north=0.1 .+ 0.8 .* rand(rng, n))

    @test_throws ArgumentError begin
        @lgm y~1 + f((east, north), spde_no_mesh) data=df family=GaussianLikelihood()
    end
    try
        @lgm y~1 + f((east, north), spde_no_mesh) data=df family=GaussianLikelihood()
    catch e
        msg = sprint(showerror, e)
        @test occursin("mesh", msg)
        @test occursin("SPDE2(mesh", msg) || occursin("inla_mesh_2d", msg)
    end
end

@testset "PR-7b lgmformula function form" begin
    rng = Random.Xoshiro(20260506)
    sq = [0.0 0.0; 1.0 0.0; 1.0 1.0; 0.0 1.0]
    mesh = inla_mesh_2d(; boundary=sq, max_edge=0.5, min_angle=25.0)
    n = 8
    df = (
        y=randn(rng, n),
        east=0.1 .+ 0.8 .* rand(rng, n),
        north=0.1 .+ 0.8 .* rand(rng, n)
    )
    spde = SPDE2(mesh)

    # Function form: 4-tuple `randoms` entry with tuple-of-Symbols col.
    model_func = lgmformula(df;
        lhs=:y,
        intercept=true,
        randoms=[((:east, :north), spde, nothing, nothing)],
        family=GaussianLikelihood())
    model_macro = @lgm y~1 + f((east, north), spde) data=df family=GaussianLikelihood()

    @test _struct_isequal(model_func.likelihoods, model_macro.likelihoods)
    @test _struct_isequal(model_func.components, model_macro.components)
    n_cols = 1 + num_vertices(mesh)
    @test size(model_func.mapping.A) == (n, n_cols)
    x = randn(rng, n_cols)
    @test (model_func.mapping.A * x)≈(model_macro.mapping.A * x) rtol=1.0e-12
end

# Phase N PR-7c — `f((east, north, time), KroneckerComponent(spde,
# time_comp))` separable space-time syntax (ADR-038). The 3-tuple
# coordinate path lowers to a sparse Khatri-Rao (row-product) design
# matrix with column layout `(s - 1) · n_t + t`, matching
# `KroneckerComponent`'s `vec(X)` (`X` of shape `(n_t × n_s)`)
# convention. The Cameletti et al. (2013) PM₁₀ space-time SPDE is the
# canonical use case.

@testset "PR-7c @macroexpand structure for 3-tuple Kronecker term" begin
    df_dummy = (y=Float64[], east=Float64[], north=Float64[],
        time=Int[])

    @testset "3-tuple f-term: comp_expr appears literally (no _wrap_term)" begin
        ex = @macroexpand @lgm y~1 + f((east, north, time), kron_dummy) data=df_dummy family=GaussianLikelihood()
        call = _find_lgm_call(ex)
        @test call !== nothing
        components = call.args[3]
        @test components.head === :tuple
        @test length(components.args) == 2
        @test _is_call_to(components.args[1], :Intercept)
        @test components.args[2] === :kron_dummy
        @test !_is_call_to(components.args[2], :_wrap_term)
    end

    @testset "3-tuple f-term: mapping is _build_design_matrix(...)" begin
        ex = @macroexpand @lgm y~1 + f((east, north, time), kron_dummy) data=df_dummy family=GaussianLikelihood()
        call = _find_lgm_call(ex)
        mapping = call.args[4]
        @test _is_call_to(mapping, :_build_design_matrix)
    end
end

@testset "PR-7c mis-pair errors" begin
    rng = Random.Xoshiro(20260506)
    sq = [0.0 0.0; 1.0 0.0; 1.0 1.0; 0.0 1.0]
    mesh = inla_mesh_2d(; boundary=sq, max_edge=0.5, min_angle=25.0)
    n = 6
    n_t = 3
    df = (y=randn(rng, n),
        east=0.1 .+ 0.8 .* rand(rng, n),
        north=0.1 .+ 0.8 .* rand(rng, n),
        time=rand(rng, 1:n_t, n))
    spde = SPDE2(mesh)
    ar1 = AR1(n_t)
    spt = KroneckerComponent(spde, ar1)

    @testset "SPDE2 + 3-tuple → user-readable KroneckerComponent error" begin
        @test_throws ArgumentError begin
            @lgm y~1 + f((east, north, time), spde) data=df family=GaussianLikelihood()
        end
        try
            @lgm y~1 + f((east, north, time), spde) data=df family=GaussianLikelihood()
        catch e
            msg = sprint(showerror, e)
            @test occursin("KroneckerComponent", msg)
            @test occursin("3-tuple", msg) || occursin("(east, north, time)", msg)
        end
    end

    @testset "KroneckerComponent + 2-tuple → user-readable 3-tuple error" begin
        @test_throws ArgumentError begin
            @lgm y~1 + f((east, north), spt) data=df family=GaussianLikelihood()
        end
        try
            @lgm y~1 + f((east, north), spt) data=df family=GaussianLikelihood()
        catch e
            msg = sprint(showerror, e)
            @test occursin("3-tuple", msg) || occursin("(east, north, time)", msg)
            @test occursin("KroneckerComponent", msg)
        end
    end
end

@testset "PR-7c KroneckerComponent error: non-SPDE2 spatial child" begin
    rng = Random.Xoshiro(20260506)
    n = 4
    n_s = 5
    n_t = 3
    df = (y=randn(rng, n),
        east=rand(rng, n),
        north=rand(rng, n),
        time=rand(rng, 1:n_t, n))
    # Spatial child is `IID`, not `SPDE2` — the extension's overload
    # accepts the type signature but raises a user-readable error.
    spt_bad = KroneckerComponent(IID(n_s), AR1(n_t))

    @test_throws ArgumentError begin
        @lgm y~1 + f((east, north, time), spt_bad) data=df family=GaussianLikelihood()
    end
    try
        @lgm y~1 + f((east, north, time), spt_bad) data=df family=GaussianLikelihood()
    catch e
        msg = sprint(showerror, e)
        @test occursin("SPDE2", msg)
        @test occursin("spatial", msg) || occursin("space", msg)
    end
end

@testset "PR-7c KroneckerComponent error: SPDE2 without retained mesh" begin
    rng = Random.Xoshiro(20260506)
    sq = [0.0 0.0; 1.0 0.0; 1.0 1.0; 0.0 1.0]
    mesh = inla_mesh_2d(; boundary=sq, max_edge=0.5, min_angle=25.0)
    spde_no_mesh = SPDE2(mesh.points, mesh.triangles)
    n_t = 3
    spt_no_mesh = KroneckerComponent(spde_no_mesh, AR1(n_t))

    n = 4
    df = (y=randn(rng, n),
        east=0.1 .+ 0.8 .* rand(rng, n),
        north=0.1 .+ 0.8 .* rand(rng, n),
        time=rand(rng, 1:n_t, n))

    @test_throws ArgumentError begin
        @lgm y~1 + f((east, north, time), spt_no_mesh) data=df family=GaussianLikelihood()
    end
    try
        @lgm y~1 + f((east, north, time), spt_no_mesh) data=df family=GaussianLikelihood()
    catch e
        msg = sprint(showerror, e)
        @test occursin("mesh", msg)
        @test occursin("SPDE2(mesh", msg) || occursin("inla_mesh_2d", msg)
    end
end

@testset "PR-7c time-coord column validation" begin
    rng = Random.Xoshiro(20260506)
    sq = [0.0 0.0; 1.0 0.0; 1.0 1.0; 0.0 1.0]
    mesh = inla_mesh_2d(; boundary=sq, max_edge=0.5, min_angle=25.0)
    spde = SPDE2(mesh)
    n = 4
    n_t = 3
    spt = KroneckerComponent(spde, AR1(n_t))

    @testset "out-of-range time index" begin
        df = (y=randn(rng, n),
            east=0.1 .+ 0.8 .* rand(rng, n),
            north=0.1 .+ 0.8 .* rand(rng, n),
            time=[1, 2, n_t + 1, 1])
        @test_throws ArgumentError begin
            @lgm y~1 + f((east, north, time), spt) data=df family=GaussianLikelihood()
        end
        try
            @lgm y~1 + f((east, north, time), spt) data=df family=GaussianLikelihood()
        catch e
            msg = sprint(showerror, e)
            @test occursin("time", msg) && occursin("outside", msg)
        end
    end

    @testset "non-integer time entry" begin
        df = (y=randn(rng, n),
            east=0.1 .+ 0.8 .* rand(rng, n),
            north=0.1 .+ 0.8 .* rand(rng, n),
            time=[1.0, 2.0, 3.0, 1.0])
        @test_throws ArgumentError begin
            @lgm y~1 + f((east, north, time), spt) data=df family=GaussianLikelihood()
        end
        try
            @lgm y~1 + f((east, north, time), spt) data=df family=GaussianLikelihood()
        catch e
            msg = sprint(showerror, e)
            @test occursin("time", msg) && occursin("integer", msg)
        end
    end
end

@testset "PR-7c Cameletti-shape (gridded) SPDE2 ⊗ AR1 roundtrip" begin
    rng = Random.Xoshiro(20260506)
    sq = [0.0 0.0; 1.0 0.0; 1.0 1.0; 0.0 1.0]
    mesh = inla_mesh_2d(; boundary=sq, max_edge=0.5, min_angle=25.0)
    n_v = num_vertices(mesh)

    n_s = 4
    n_t = 5
    n = n_s * n_t
    station_east = 0.1 .+ 0.8 .* rand(rng, n_s)
    station_north = 0.1 .+ 0.8 .* rand(rng, n_s)
    east = Vector{Float64}(undef, n)
    north = Vector{Float64}(undef, n)
    time = Vector{Int}(undef, n)
    # Time-fast-within-space ordering — matches `KroneckerComponent`'s
    # `vec(X)` flattening (`x[(s - 1) · n_t + t]`).
    @inbounds for s in 1:n_s
        for t in 1:n_t
            i = (s - 1) * n_t + t
            east[i] = station_east[s]
            north[i] = station_north[s]
            time[i] = t
        end
    end
    df = (y=randn(rng, n), east=east, north=north, time=time)
    pc = PCMatern(range_U=0.5, range_α=0.05, sigma_U=1.0, sigma_α=0.01)
    spde = SPDE2(mesh; pc=pc)
    ar1 = AR1(n_t)
    spt = KroneckerComponent(spde, ar1)

    model_macro = @lgm y~1 + f((east, north, time), spt) data=df family=GaussianLikelihood()

    # Hand-form: per-station spatial projector × I_{n_t} (Cameletti
    # canonical form). The first n_t rows of the per-obs spatial
    # projector all share station 1, so we can grab one row per station.
    A_int = sparse(ones(n, 1))
    P_per_obs = MeshProjector(mesh, [df.east df.north])
    A_per_obs = SparseMatrixCSC{Float64, Int}(P_per_obs.A)
    A_space_j = A_per_obs[1:n_t:n, :]   # one row per station
    A_kron = kron(A_space_j, sparse(1.0I, n_t, n_t))
    A_hand = hcat(A_int, A_kron)
    model_hand = LatentGaussianModel(
        GaussianLikelihood(),
        (Intercept(), spt),
        A_hand
    )

    @test _struct_isequal(model_macro.likelihoods, model_hand.likelihoods)
    @test _struct_isequal(model_macro.components, model_hand.components)

    n_cols = 1 + n_v * n_t
    @test size(model_macro.mapping.A) == (n, n_cols)
    @test size(model_hand.mapping.A) == (n, n_cols)
    # `MeshProjector` is invoked twice independently (once by the macro,
    # once in the hand form). Per ADR-036 the two projector matrices
    # agree to fp roundoff but not bitwise — compare via `A * x`.
    x = randn(rng, n_cols)
    @test (model_macro.mapping.A * x)≈(model_hand.mapping.A * x) rtol=1.0e-12
end

@testset "PR-7c per-obs (non-gridded) SPDE2 ⊗ AR1 roundtrip" begin
    rng = Random.Xoshiro(20260506)
    sq = [0.0 0.0; 1.0 0.0; 1.0 1.0; 0.0 1.0]
    mesh = inla_mesh_2d(; boundary=sq, max_edge=0.5, min_angle=25.0)
    n_v = num_vertices(mesh)

    n = 14
    n_t = 4
    df = (
        y=randn(rng, n),
        east=0.1 .+ 0.8 .* rand(rng, n),
        north=0.1 .+ 0.8 .* rand(rng, n),
        time=rand(rng, 1:n_t, n)
    )
    spde = SPDE2(mesh)
    spt = KroneckerComponent(spde, AR1(n_t))

    model_macro = @lgm y~1 + f((east, north, time), spt) data=df family=GaussianLikelihood()

    # Hand-built Khatri-Rao reference: same column convention as the
    # extension implementation, built independently via `findnz`.
    A_int = sparse(ones(n, 1))
    P = MeshProjector(mesh, [df.east df.north])
    A_space = SparseMatrixCSC{Float64, Int}(P.A)
    rows_s, cols_s, vals_s = findnz(A_space)
    new_cols = [(cols_s[k] - 1) * n_t + df.time[rows_s[k]] for k in eachindex(rows_s)]
    A_st = sparse(rows_s, new_cols, vals_s, n, n_v * n_t)
    A_hand = hcat(A_int, A_st)
    model_hand = LatentGaussianModel(
        GaussianLikelihood(),
        (Intercept(), spt),
        A_hand
    )

    @test _struct_isequal(model_macro.components, model_hand.components)
    n_cols = 1 + n_v * n_t
    @test size(model_macro.mapping.A) == (n, n_cols)
    x = randn(rng, n_cols)
    @test (model_macro.mapping.A * x)≈(model_hand.mapping.A * x) rtol=1.0e-12
end

@testset "PR-7c no-intercept: f((east, north, time), spt) only" begin
    rng = Random.Xoshiro(20260506)
    sq = [0.0 0.0; 1.0 0.0; 1.0 1.0; 0.0 1.0]
    mesh = inla_mesh_2d(; boundary=sq, max_edge=0.5, min_angle=25.0)
    n_v = num_vertices(mesh)
    n = 6
    n_t = 3
    df = (
        y=randn(rng, n),
        east=0.1 .+ 0.8 .* rand(rng, n),
        north=0.1 .+ 0.8 .* rand(rng, n),
        time=rand(rng, 1:n_t, n)
    )
    spde = SPDE2(mesh)
    spt = KroneckerComponent(spde, AR1(n_t))

    model_macro = @lgm y~0 + f((east, north, time), spt) data=df family=GaussianLikelihood()
    @test length(model_macro.components) == 1
    @test _struct_isequal(model_macro.components[1], spt)
    @test size(model_macro.mapping.A) == (n, n_v * n_t)
end

@testset "PR-7c lgmformula function form (3-tuple Kronecker)" begin
    rng = Random.Xoshiro(20260506)
    sq = [0.0 0.0; 1.0 0.0; 1.0 1.0; 0.0 1.0]
    mesh = inla_mesh_2d(; boundary=sq, max_edge=0.5, min_angle=25.0)
    n_v = num_vertices(mesh)
    n = 10
    n_t = 4
    df = (
        y=randn(rng, n),
        east=0.1 .+ 0.8 .* rand(rng, n),
        north=0.1 .+ 0.8 .* rand(rng, n),
        time=rand(rng, 1:n_t, n)
    )
    spde = SPDE2(mesh)
    spt = KroneckerComponent(spde, AR1(n_t))

    model_func = lgmformula(df;
        lhs=:y,
        intercept=true,
        randoms=[((:east, :north, :time), spt, nothing, nothing)],
        family=GaussianLikelihood())
    model_macro = @lgm y~1 + f((east, north, time), spt) data=df family=GaussianLikelihood()

    @test _struct_isequal(model_func.likelihoods, model_macro.likelihoods)
    @test _struct_isequal(model_func.components, model_macro.components)
    n_cols = 1 + n_v * n_t
    @test size(model_func.mapping.A) == (n, n_cols)
    x = randn(rng, n_cols)
    @test (model_func.mapping.A * x)≈(model_macro.mapping.A * x) rtol=1.0e-12
end

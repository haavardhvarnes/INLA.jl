# Phase N PR-7b — `f((east, north), SPDE2(mesh))` tuple-coordinate
# syntax (ADR-037 / ADR-039). The parser accepts a tuple of column
# symbols as the first argument of `f(...)`; lengths 2 (spatial) and 3
# (space-time, reserved for PR-7c) are accepted, others rejected.
# At runtime, the term dispatches through `_build_spatial_block`,
# overloaded for `SPDE2` in `LGMFormulaINLASPDEExt` (loaded
# automatically when both `LGMFormula` and `INLASPDE` are imported).
#
# Tests cover: parse-time errors (tuple arity, non-Symbol entries,
# `replicate`/`group` rejection), macroexpand structure (no
# `_wrap_term` wrap on tuple-coord terms), the missing-extension error
# path (non-SPDE component with tuple syntax raises a user-readable
# error), the SPDE2-via-mesh roundtrip against an explicit
# `LatentGaussianModel(...)` constructor, and the no-mesh error from
# the `(points, triangles)` `SPDE2` constructor (ADR-036).

using INLASPDE: SPDE2, MeshProjector, INLAMesh, inla_mesh_2d,
    PCMatern, num_vertices

@testset "PR-7b parse-time errors for tuple-coord f-term" begin
    df = (y = Float64[], east = Float64[], north = Float64[],
        time = Float64[], idx = Int[], id = Int[], grp = Int[])

    @testset "tuple of length 1 rejected" begin
        ex = try
            @eval @lgm y ~ 1 + f((east,), IID(5)) data=$df family=GaussianLikelihood()
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
            @eval @lgm y ~ 1 + f((east, north, time, idx), IID(5)) data=$df family=GaussianLikelihood()
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
            @eval @lgm y ~ 1 + f((east, 5), IID(5)) data=$df family=GaussianLikelihood()
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
            @eval @lgm y ~ 1 + f((east, north), IID(5); replicate = id) data=$df family=GaussianLikelihood()
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
            @eval @lgm y ~ 1 + f((east, north), IID; group = grp) data=$df family=GaussianLikelihood()
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
    df_dummy = (y = Float64[], east = Float64[], north = Float64[])

    # We use a placeholder symbol `spde_dummy` for the component slot —
    # we only care that the AST shape is right. The runtime is not
    # invoked here, so the symbol need not resolve.
    @testset "tuple-coord f-term: comp_expr appears literally (no _wrap_term)" begin
        ex = @macroexpand @lgm y ~ 1 + f((east, north), spde_dummy) data=df_dummy family=GaussianLikelihood()
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
        ex = @macroexpand @lgm y ~ 1 + f((east, north), spde_dummy) data=df_dummy family=GaussianLikelihood()
        call = _find_lgm_call(ex)
        mapping = call.args[4]
        @test _is_call_to(mapping, :_build_design_matrix)
    end
end

@testset "PR-7b missing-extension error (non-SPDE component)" begin
    rng = Random.Xoshiro(20260506)
    n = 5
    df = (y = randn(rng, n),
        east = rand(rng, n),
        north = rand(rng, n))

    # `IID(5)` is *not* an SPDE component, but the parser accepts the
    # tuple-coord syntax syntactically. At runtime, dispatch falls
    # through to the default `_build_spatial_block`, which raises a
    # user-readable error pointing at the INLASPDE extension.
    @test_throws ArgumentError begin
        @lgm y ~ 1 + f((east, north), IID(5)) data=df family=GaussianLikelihood()
    end
    try
        @lgm y ~ 1 + f((east, north), IID(5)) data=df family=GaussianLikelihood()
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
        logzinc = randn(rng, n),
        dist = randn(rng, n),
        east = 0.1 .+ 0.8 .* rand(rng, n),
        north = 0.1 .+ 0.8 .* rand(rng, n),
    )
    pc = PCMatern(range_U=0.5, range_α=0.05, sigma_U=1.0, sigma_α=0.01)
    spde = SPDE2(mesh; pc=pc)

    @testset "Meuse-shape: 1 + dist + f((east, north), spde)" begin
        model_macro = @lgm logzinc ~ 1 + dist + f((east, north), spde) data=df family=GaussianLikelihood()

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
            A_hand,
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
        model_macro = @lgm logzinc ~ 0 + f((east, north), spde) data=df family=GaussianLikelihood()
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
    df = (y = randn(rng, n),
        east = 0.1 .+ 0.8 .* rand(rng, n),
        north = 0.1 .+ 0.8 .* rand(rng, n))

    @test_throws ArgumentError begin
        @lgm y ~ 1 + f((east, north), spde_no_mesh) data=df family=GaussianLikelihood()
    end
    try
        @lgm y ~ 1 + f((east, north), spde_no_mesh) data=df family=GaussianLikelihood()
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
        y = randn(rng, n),
        east = 0.1 .+ 0.8 .* rand(rng, n),
        north = 0.1 .+ 0.8 .* rand(rng, n),
    )
    spde = SPDE2(mesh)

    # Function form: 4-tuple `randoms` entry with tuple-of-Symbols col.
    model_func = lgmformula(df;
        lhs = :y,
        intercept = true,
        randoms = [((:east, :north), spde, nothing, nothing)],
        family = GaussianLikelihood())
    model_macro = @lgm y ~ 1 + f((east, north), spde) data=df family=GaussianLikelihood()

    @test _struct_isequal(model_func.likelihoods, model_macro.likelihoods)
    @test _struct_isequal(model_func.components, model_macro.components)
    n_cols = 1 + num_vertices(mesh)
    @test size(model_func.mapping.A) == (n, n_cols)
    x = randn(rng, n_cols)
    @test (model_func.mapping.A * x)≈(model_macro.mapping.A * x) rtol=1.0e-12
end

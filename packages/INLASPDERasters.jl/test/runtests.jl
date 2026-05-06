using Test
using INLASPDERasters

@testset "INLASPDERasters.jl" begin
    @testset "smoke — module loads" begin
        @test isdefined(INLASPDERasters, :INLASPDERasters)
        @test isdefined(INLASPDERasters, :extract_at_mesh)
    end

    @testset "M1 — extraction" begin
        include("regression/test_extract_synthetic.jl")
    end

    @testset "M2 — prediction" begin
        include("regression/test_predict_synthetic.jl")
    end

    @testset "M3 — uncertainty surfaces" begin
        include("regression/test_quantile_rasters.jl")
    end

    @testset "Phase O PR-1 — predict_raster(model, res, template)" begin
        include("regression/test_predict_model.jl")
    end

    @testset "Phase O PR-1 — extract_at_mesh CRS handling" begin
        include("regression/test_extract_crs.jl")
    end

    @testset "Phase O PR-2 — sample-based predict_raster + Exceedance" begin
        include("regression/test_predict_sample.jl")
    end

    @testset "Phase O PR-4 — Meuse R-INLA predict oracle (1e-10 gate)" begin
        include("oracle/test_meuse_predict.jl")
    end

    @testset "Phase O PR-5 — SPDE2NonStationary + KroneckerComponent dispatch" begin
        include("regression/test_predict_phase_o_pr5.jl")
    end

    @testset "Quality" begin
        include("quality/test_aqua.jl")
        include("quality/test_jet.jl")
    end
end

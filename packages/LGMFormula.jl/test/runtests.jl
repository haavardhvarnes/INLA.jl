using LGMFormula
using LatentGaussianModels
using Random
using SparseArrays
using Tables
using Test

include("test_utils.jl")

@testset "LGMFormula" begin
    include("regression/test_macroexpand.jl")
    include("regression/test_roundtrip.jl")
    include("regression/test_error_messages.jl")
    include("regression/test_components.jl")
    include("regression/test_multi_f.jl")
    include("regression/test_multi_likelihood.jl")
end

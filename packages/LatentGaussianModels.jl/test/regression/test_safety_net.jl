# Phase M PR-4 — PD-failure safety net regression suite.
#
# Verifies the four behaviours that together harden the LGM Newton hot
# path against pathological θ inputs (extreme κ on SPDE precisions,
# overflowing Gaussian τ_obs on log-precision scale, NaN warm-starts
# forwarded from a neighbouring CCD point):
#
#   1. `_is_bad_theta_failure` accepts the documented numerical
#      exceptions and rejects everything else (so genuine bugs surface
#      unmasked).
#   2. `_neg_log_posterior_θ` returns the smooth quadratic penalty
#      `1e10 + 1e3·‖θ‖²` when the inner Laplace step throws / returns
#      a non-finite log-marginal.
#   3. `laplace_mode` resets a NaN/Inf warm-start before the first
#      Newton step instead of propagating the contamination.
#   4. `_safe_inverse_hessian` falls back to the identity covariance
#      when the FD Hessian at θ̂ contains non-finite entries.
#
# Plus an end-to-end `inla()` fit started from a deliberately extreme
# `θ0` that lands inside the penalty region — LBFGS escapes via the
# smooth quadratic and the integration stage produces a finite
# `log_marginal`.

using LatentGaussianModels: GaussianLikelihood, Intercept, IID,
                            LatentGaussianModel, laplace_mode, Laplace,
                            INLA, fit
using LinearAlgebra
using SparseArrays
using Random

const LGM = LatentGaussianModels

@testset "_is_bad_theta_failure — bad-θ classification" begin
    @test LGM._is_bad_theta_failure(LinearAlgebra.PosDefException(1))
    @test LGM._is_bad_theta_failure(LinearAlgebra.SingularException(1))
    @test LGM._is_bad_theta_failure(DomainError(NaN, "test"))
    @test !LGM._is_bad_theta_failure(ArgumentError("not a numerical failure"))
    @test !LGM._is_bad_theta_failure(MethodError(sin, ()))
    @test !LGM._is_bad_theta_failure(ErrorException("oops"))
end

# Gaussian + Intercept fixture. At `log τ_obs ≈ 1500` the Hessian
# `H = Q + τ_obs A'A` blows past Float64 range and the inner-Newton
# Cholesky throws — this is the canonical bad-θ regime the safety net
# is responsible for.
function _build_safety_net_model(n::Int=30, rng_seed::Int=0)
    rng = Random.Xoshiro(rng_seed)
    y = randn(rng, n)
    c = Intercept()
    A = sparse(ones(n, 1))
    ℓ = GaussianLikelihood()
    model = LatentGaussianModel(ℓ, (c,), A)
    return model, y
end

@testset "_neg_log_posterior_θ — smooth penalty at extreme θ" begin
    model, y = _build_safety_net_model()
    f = LGM._neg_log_posterior_θ(model, y, Laplace())

    # Sane θ — finite, well below the 1e10 penalty floor.
    val_normal = f([0.0], nothing)
    @test isfinite(val_normal)
    @test val_normal < 1.0e9

    # Extreme θ — Cholesky on Q + exp(1500)·A'A fails; closure must
    # return the smooth penalty rather than raising or returning Inf.
    θ_extreme = [1500.0]
    val_extreme = f(θ_extreme, nothing)
    @test isfinite(val_extreme)
    @test val_extreme > 1.0e9
    @test val_extreme≈1.0e10 + 1.0e3 * sum(abs2, θ_extreme) atol=1.0e-6
end

@testset "laplace_mode — non-finite warm-start reset" begin
    model, y = _build_safety_net_model()
    # Pass a NaN warm-start; the inner Newton must reset to zeros and
    # converge to the same posterior mode as the cold start.
    res_nan = laplace_mode(model, y, [0.0]; x0=fill(NaN, model.n_x))
    @test res_nan.converged
    @test all(isfinite, res_nan.mode)
    @test isfinite(res_nan.log_marginal)

    res_inf = laplace_mode(model, y, [0.0]; x0=fill(Inf, model.n_x))
    @test res_inf.converged
    @test all(isfinite, res_inf.mode)

    # Reference: cold start from zeros should give the same answer.
    res_cold = laplace_mode(model, y, [0.0])
    @test res_nan.mode≈res_cold.mode rtol=1.0e-8
    @test res_inf.mode≈res_cold.mode rtol=1.0e-8
end

@testset "_safe_inverse_hessian — non-finite fallback" begin
    # NaN entry → identity fallback + warning.
    H_nan = [NaN 0.0; 0.0 1.0]
    Σ_nan = @test_logs (:warn, r"non-finite") LGM._safe_inverse_hessian(H_nan)
    @test Σ_nan == Matrix{Float64}(I, 2, 2)

    # Inf entry → identity fallback + warning.
    H_inf = [Inf 0.0; 0.0 1.0]
    Σ_inf = @test_logs (:warn, r"non-finite") LGM._safe_inverse_hessian(H_inf)
    @test Σ_inf == Matrix{Float64}(I, 2, 2)

    # Well-formed PD H goes through the eigen-floor path unchanged.
    H_pd = [4.0 0.0; 0.0 9.0]
    Σ_pd = LGM._safe_inverse_hessian(H_pd)
    @test Σ_pd≈[0.25 0.0; 0.0 1/9] atol=1.0e-12
end

@testset "INLA — recovers from extreme θ0 via smooth penalty" begin
    model, y = _build_safety_net_model(50, 20260505)
    # `θ0 = [1500.0]` lands inside the bad-θ region. The penalty is
    # quadratic in θ, so LBFGS sees a finite gradient pointing back to
    # the origin, escapes the penalty region, and converges to the
    # genuine θ-mode. Without the safety net, the first Optim step
    # would have crashed on a non-finite objective from the inner
    # Newton's Cholesky failure.
    res = fit(model, y, INLA(; θ0=[1500.0], int_strategy=:grid))
    @test isfinite(res.log_marginal)
    @test all(isfinite, res.x_mean)
    @test all(isfinite, res.x_var)
    @test all(≥(0), res.x_var)
    # The recovered θ̂ should not still be in the penalty region.
    @test abs(res.θ̂[1]) < 50.0
end

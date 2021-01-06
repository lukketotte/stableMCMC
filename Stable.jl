module AlphaStableDistribution

using Distributions, LinearAlgebra, Statistics, SpecialFunctions
using AlphaStableDistributions, QuadGK

import Distributions.pdf

function CF(d::AlphaStable, t)
    μ, σ = d.location, d.scale
    if d.α != 1
        k = d.α - 1 + sign(1-d.α)
        exp(-σ^d.α * abs(t)^d.α * (-im * d.β * sign(t) * π / 2 * k) + im * μ * t)
    else
        exp(-σ * abs(t) * (π / 2 + im*d.β*sign(t)*log(abs(t))) + im * μ * t)
    end
end

function invertCF(d::AlphaStable, t, x)
    real(CF(d, t) * exp(-im * x * t) / 2)
end

function pdf(d::AlphaStable, x)
    integral, err = quadgk(t -> invertCF(d, t, x), 2., Inf, rtol=1e-10)
    integral
end


quadgk(t -> invertCF(d, t, 0.), 0., Inf, rtol=1e-10)


CF(d, 0.)
pdf(AlphaStable(2., 0. , 1., 1.), 0.5)

d.α = 1

real(1 + 2*im)

using CSV, DataFrames, AlphaStableDistributions, Distributions
using KernelDensity, Plots, PlotThemes
theme(:juno)

y = CSV.File("data/nsa_y.csv") |> DataFrame
x = convert(Array{Float64, 1}, y.x)
nsa = x[x.>0]

y = CSV.File("data/cm_os.csv"; missingstrings = ["NA"]) |> DataFrame
x = convert(Array{Float64, 1}, y[completecases(y), :ged_best_os])
cm = x[(x .> 1) .& (x .< 10)]

maximum(cm)

d1 = fit(AlphaStable, cm)
d2 = fit(AlphaStable, nsa)

k = kde(cm)
x = range(minimum(cm)-1, maximum(cm)+1, length = 2000)
plot(x, pdf(k, x))


d = AlphaStable(0.5, 0.5 ,1., 0.)
d.α

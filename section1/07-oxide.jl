### Oxide --- nested ###
home = homedir()
datafile = "$home/work/projects/testing-random-effects/data/oxide.csv"
# From R package nlme, not included in RDatasets
resultspath = "$home/work/projects/testing-random-effects/julia/results"
savename = resultspath * "/oxideanalysis"


using 
  CSV,
  DataFrames,
  LinearAlgebra,
  MixedModels,
  CategoricalArrays,
  varcomptest,
  StatsPlots,
  Statistics,
  RDatasets

Oxide = CSV.read(datafile, DataFrame);
# Have to create the nested effect manually
Oxide.WaferLot = string.(Oxide.Lot) .* ":" .* string.(Oxide.Wafer)
Oxide.Lot = categorical(Oxide.Lot, ordered = false, levels = unique(Oxide.Lot));
Oxide.Wafer = categorical(Oxide.Wafer, ordered = false, levels = unique(Oxide.Wafer));
Oxide.WaferLot = categorical(Oxide.WaferLot, ordered = false, levels = unique(Oxide.WaferLot));
Oxide.Thickness = Float64.(Oxide.Thickness);


ff = @formula(Thickness ~ (1 | Lot) + (1 | WaferLot))
control = VarCompControl(NewtonControl(verbose = false), 100000)
A = Matrix([1., -1.]')
_ = varcompmodel(ff, Oxide, control = VarCompControl(NewtonControl(), 0), A = nothing)
vmod = varcompmodel(ff, Oxide, control = control, A = A)

out = DataFrame(
  hcat([vmod.bootresults.mle, vmod.bootresults.lrt]...),
  vcat([vmod.vr.names, ["lrt"]]...)
)

CSV.write("$savename.csv", out)

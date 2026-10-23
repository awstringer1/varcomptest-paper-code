### Analyze the Pastes data ###

home = homedir()
resultspath = "$home/work/projects/testing-random-effects/julia/results"
savename = resultspath * "/pastesanalysis"


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

Pastes = dataset("lme4", "Pastes");
Pastes.Batch = categorical(Pastes.Batch, ordered = false, levels = unique(Pastes.Batch));
Pastes.Cask = categorical(Pastes.Cask, ordered = false, levels = unique(Pastes.Cask));
Pastes.Sample = categorical(Pastes.Sample, ordered = false, levels = unique(Pastes.Sample));
Pastes.BatchCask = string.(Pastes.Batch) .* ":" .* string.(Pastes.Cask)
Pastes.BatchCask = categorical(Pastes.BatchCask, ordered = false, levels = unique(Pastes.BatchCask));


ff = @formula(Strength ~ (1 | Batch) + (1 | BatchCask))
# JIT
_ = varcompmodel(ff, Pastes, control = VarCompControl(NewtonControl(), 0), A = nothing)
control = VarCompControl(NewtonControl(), 100000)
A = Matrix([-1., 1.]')
vmod = varcompmodel(ff, Pastes, control = control, A = A)

out = DataFrame(
  hcat([vmod.bootresults.mle, vmod.bootresults.lrt]...),
  vcat([vmod.vr.names, ["lrt"]]...)
)

CSV.write("$savename.csv", out)

### Analyze the Pastes data ###

home = homedir()
resultspath = home
savename = resultspath * "/pastesanalysis"


using 
  MixedModels,
  CategoricalArrays,
  varcomptest,
  RDatasets

Pastes = dataset("lme4", "Pastes");
Pastes.Batch = categorical(Pastes.Batch, ordered = false, levels = unique(Pastes.Batch));
Pastes.Cask = categorical(Pastes.Cask, ordered = false, levels = unique(Pastes.Cask));
Pastes.Sample = categorical(Pastes.Sample, ordered = false, levels = unique(Pastes.Sample));
Pastes.BatchCask = string.(Pastes.Batch) .* ":" .* string.(Pastes.Cask)
Pastes.BatchCask = categorical(Pastes.BatchCask, ordered = false, levels = unique(Pastes.BatchCask));

ff = @formula(Strength ~ (1 | Batch) + (1 | BatchCask))
lmod = lmm(ff, Pastes)
control = VarCompControl(NewtonControl(verbose = false), 1000)
A = Matrix([-1., 1.]')
vmod = varcompmodel(ff, Pastes, control = control, A = A)

out = DataFrame(
  hcat([vmod.bootresults.mle, vmod.bootresults.lrt]...),
  vcat([vmod.vr.names, ["lrt"]]...)
)

CSV.write("$savename.csv", out)

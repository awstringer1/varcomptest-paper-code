### Analyze the Penicillin data ###

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

Penicillin = dataset("lme4", "Penicillin");
Penicillin.Plate = categorical(Penicillin.Plate, ordered = false, levels = unique(Penicillin.Plate));
Penicillin.Sample = categorical(Penicillin.Sample, ordered = false, levels = unique(Penicillin.Sample));
Penicillin.Diameter = Float64.(Penicillin.Diameter);


ff = @formula(Diameter ~ (1 | Plate) + (1 | Sample))
lmod = lmm(ff, Penicillin)
control = VarCompControl(NewtonControl(verbose = false), 1000)
A = Matrix([-1., 1.]')
vmod = varcompmodel(ff, Penicillin, control = control, A = A)
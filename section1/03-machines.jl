### Analyze the Machines data ###
home = homedir()
datafile = "$home/work/projects/testing-random-effects/data/machines.csv"
# From R package nlme, not included in RDatasets

using 
  CSV,
  DataFrames,
  MixedModels,
  CategoricalArrays,
  varcomptest

Machines = CSV.read(datafile, DataFrame);
Machines.Worker = categorical(Machines.Worker, ordered = false, levels = unique(Machines.Worker));
Machines.Machine = categorical(Machines.Machine, ordered = false, levels = unique(Machines.Machine));

ff = @formula(score ~ (1 | Machine) + (1 | Worker))
control = VarCompControl(NewtonControl(verbose = false), 100000)
A = Matrix([1., -1.]')
# JIT
_ = varcompmodel(ff, Machines, control = VarCompControl(NewtonControl(), 0), A = A)
vmod = varcompmodel(ff, Machines, control = control, A = A)
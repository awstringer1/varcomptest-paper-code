### Analyze the Machines data ###
home = homedir() # Assumed to be the root of the repo
datafile = "$home/data/machines.csv"
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
lmod = lmm(ff, Machines)
control = VarCompControl(NewtonControl(verbose = false), 1000)
A = Matrix([1., -1.]')
vmod = varcompmodel(ff, Machines, control = control, A = A)

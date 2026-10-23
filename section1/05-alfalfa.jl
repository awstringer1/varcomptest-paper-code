### Analyze the Alfalfa data ###
home = homedir()
datafile = "$home/work/projects/testing-random-effects/data/alfalfa.csv"
# From R package nlme, not included in RDatasets

using 
  CSV,
  DataFrames,
  LinearAlgebra,
  MixedModels,
  CategoricalArrays,
  varcomptest

Alfalfa = CSV.read(datafile, DataFrame);
Alfalfa.Block = categorical(Alfalfa.Block, ordered = false, levels = unique(Alfalfa.Block));
Alfalfa.Variety = categorical(Alfalfa.Variety, ordered = false, levels = unique(Alfalfa.Variety));
Alfalfa.Date = categorical(Alfalfa.Date, ordered = false, levels = unique(Alfalfa.Date));
Alfalfa.BlockDate = string.(Alfalfa.Block) .* ":" .* string.(Alfalfa.Date)
Alfalfa.BlockDate = categorical(Alfalfa.BlockDate, ordered = false, levels = unique(Alfalfa.BlockDate));


ff = @formula(Yield ~ Variety + (1 | Block) + (1 | BlockDate))
control = VarCompControl(NewtonControl(verbose = false), 100000)
A = Matrix([1., -1.]')
_ = varcompmodel(ff, Alfalfa, control = VarCompControl(NewtonControl(), 0), A = nothing)
vmod = varcompmodel(ff, Alfalfa, control = control, A = A)

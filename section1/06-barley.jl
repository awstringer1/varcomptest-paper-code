### Analyze the barley data ###

using 
  CSV,
  DataFrames,
  LinearAlgebra,
  MixedModels,
  CategoricalArrays,
  varcomptest,
  RDatasets

barley = dataset("lattice", "barley");

ff = @formula(Yield ~ (1 | Variety) + (1 | Site))
control = VarCompControl(NewtonControl(verbose = false), 100000)
A = Matrix([-1., 1.]')
_ = varcompmodel(ff, barley, control = VarCompControl(NewtonControl(), 0), A = nothing)
vmod = varcompmodel(ff, barley, control = control, A = A)

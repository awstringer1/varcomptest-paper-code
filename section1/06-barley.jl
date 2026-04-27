### Analyze the barley data ###

using 
  MixedModels,
  CategoricalArrays,
  varcomptest,
  RDatasets

barley = dataset("lattice", "barley");

ff = @formula(Yield ~ (1 | Variety) + (1 | Site))
lmod = lmm(ff, barley)
control = VarCompControl(NewtonControl(verbose = false), 1000)
A = Matrix([-1., 1.]')
vmod = varcompmodel(ff, barley, control = control, A = A)
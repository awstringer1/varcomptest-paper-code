### Analyze the Oats data ###

using 
  MixedModels,
  CategoricalArrays,
  varcomptest,
  RDatasets

oats = dataset("MASS", "oats");
oats.B = categorical(oats.B, ordered = false, levels = unique(oats.B));
oats.N = categorical(oats.N, ordered = false, levels = unique(oats.N));
oats.V = categorical(oats.V, ordered = false, levels = unique(oats.V));
oats.Y = Float64.(oats.Y);
oats.BinV = string.(oats.B) .* ":" .* string.(oats.V)
oats.BinV = categorical(oats.BinV, ordered = false, levels = unique(oats.BinV));

ff = @formula(Y ~ N + V + (1 | B) + (1 | BinV))
lmod = lmm(ff, oats)
control = VarCompControl(NewtonControl(verbose = false), 1000)
A = Matrix([1., -1.]')
vmod = varcompmodel(ff, oats, control = control, A = A)

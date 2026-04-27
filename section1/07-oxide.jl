### Oxide --- nested ###

using 
  MixedModels,
  CategoricalArrays,
  varcomptest,
  RDatasets

Oxide = CSV.read(datafile, DataFrame);
Oxide.WaferLot = string.(Oxide.Lot) .* ":" .* string.(Oxide.Wafer)
Oxide.Lot = categorical(Oxide.Lot, ordered = false, levels = unique(Oxide.Lot));
Oxide.Wafer = categorical(Oxide.Wafer, ordered = false, levels = unique(Oxide.Wafer));
Oxide.WaferLot = categorical(Oxide.WaferLot, ordered = false, levels = unique(Oxide.WaferLot));
Oxide.Thickness = Float64.(Oxide.Thickness);


ff = @formula(Thickness ~ (1 | Lot) + (1 | WaferLot))

lmod = lmm(ff, Oxide)
control = VarCompControl(NewtonControl(verbose = false), 1000)
A = Matrix([1., -1.]')
vmod = varcompmodel(ff, Oxide, control = control, A = A)

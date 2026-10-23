### Analyze the Penicillin data ###

home = homedir()
resultspath = "$home/work/projects/testing-random-effects/julia/results"
savename = resultspath * "/penicillinanalysis"


using 
  CSV,
  DataFrames,
  LinearAlgebra,
  MixedModels,
  CategoricalArrays,
  varcomptest,
  RDatasets

Penicillin = dataset("lme4", "Penicillin");
Penicillin.Plate = categorical(Penicillin.Plate, ordered = false, levels = unique(Penicillin.Plate));
Penicillin.Sample = categorical(Penicillin.Sample, ordered = false, levels = unique(Penicillin.Sample));
Penicillin.Diameter = Float64.(Penicillin.Diameter);


ff = @formula(Diameter ~ (1 | Plate) + (1 | Sample))
control = VarCompControl(NewtonControl(verbose = false), 100000)
A = Matrix([-1., 1.]')
_ = varcompmodel(ff, Penicillin, control = VarCompControl(NewtonControl(), 0), A = nothing)
vmod = varcompmodel(ff, Penicillin, control = control, A = A)

out = DataFrame(
  hcat([vmod.bootresults.mle, vmod.bootresults.lrt]...),
  vcat([vmod.vr.names, ["lrt"]]...)
)

CSV.write("$savename.csv", out)

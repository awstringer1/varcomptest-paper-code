### Simulations: Optimization, Crossed Random Effects ###

## Paths ----

resultspath = homedir() * "/"
if !isdir(resultspath)
  throw(SystemError("Could not access directory: $resultspath"))
end

if length(ARGS) == 0
  DATE = "2026"
  SIMVERSION = "1"
  nsim = 2
  bootsamps = 300
  sampling = 1
  SIMIDX = 1
  JOBID = 1
  TASKID = 1
else
  nsim = parse(Int64, ARGS[1])
  DATE = ARGS[2]
  SIMVERSION = ARGS[3]
  JOBID = ARGS[4]
  SIMIDX = ARGS[5]
  TASKID = ARGS[6]
end
simname = "sims-crossed-opt-$DATE-v$SIMVERSION-job$JOBID-sim$SIMIDX-task$TASKID"
savename = resultspath*simname
println("simname = ", simname)

## Packages ----

using 
  varcomptest,
  DataFrames,
  MixedModels,
  Distributions,
  CategoricalArrays,
  JLD2,
  StatsPlots,
  CSV,
  Base.Threads,
  LinearAlgebra

nt = Threads.nthreads();

## Parameters ----

m = [20, 50, 100]
n = [20, 50, 100]
r = [-1., 0., .5] # This is now the sampling correlation not the reps
tau1 = [0., .01, .02, .05, .1]
taudiff = [0., .01, .05, .1, 1., 2.] # tau2 = tau1 + taudiff
params = Dict(
  :beta => [0., 1., -2.],
  :sigma => 1.,
  :B => 0
)

## Functions ----

simulate_data = function(m::Int64, n::Int64, r::Float64, tau::Vector{Float64}, params::Dict)
  p = length(params[:beta])

  idx = Iterators.product(1:m, 1:n)
  N = prod(size(idx))
  u1 = rand(Normal(0, params[:sigma] * sqrt(tau[1])), m);
  u2 = rand(Normal(0, params[:sigma] * sqrt(tau[2])), n);

  if r == -1
    # Use all pairs
    idx1 = vcat([i[1] for i in idx]...)
    idx2 = vcat([i[2] for i in idx]...)
  elseif 0. <= r && r < 1.
    # Sample the inclusion indices from a correlated discrete uniform with correlation r
    cor = Float64.(ones(2, 2))
    cor[1, 2] = cor[2, 1] = r
    mvn = MvNormal([0., 0.], cor)
    mvnsamps = rand(mvn, N)
    stdnormal = Normal(0, 1)
    idx1 = Int64.([floor(m * cdf(stdnormal, mvnsamps[1, j]) + 1.) for j in 1:N])
    idx2 = Int64.([floor(n * cdf(stdnormal, mvnsamps[2, j]) + 1.) for j in 1:N])
  else
    ss = params[:sampling]
    throw("Unknown sampling strategy $ss, options are -1 or 0 <= r < 1")
  end

  
  dat = DataFrame(
    y = zeros(N),
    subject = idx1,
    item = idx2
  );
  # dat = unique(dat)

  
  X = hcat(ones(N), rand(Uniform(0, 1), N, p-1));
  mu = X * params[:beta];

  for j in 1:(p-1)
    dat[!, "x$j"] = X[:, j + 1]
  end;
  
  for i in 1:N
    dat.y[i] = rand(Normal(mu[i] + u1[dat.subject[i]] + u2[dat.item[i]], params[:sigma]))
  end;
  
  # convert to categorical
  dat.subject = categorical(dat.subject, ordered = true, levels = unique(dat.subject));
  dat.item = categorical(dat.item, ordered = true, levels = unique(dat.item));

  dat
end

# tau2 = tau1 + taudiff
allcombinations = repeat(DataFrame(
    [(; m=i[1], n=i[2], r=i[3], tau1=i[4], tau2=i[4] + i[5])
     for i in Iterators.product(m, n, r, tau1, taudiff)]
), nsim)
allcombinations.sim = 1:size(allcombinations, 1)

# df = allcombinations[1,:]

dosim = function(df::DataFrameRow{DataFrame, DataFrames.Index}; params::Dict = params) 
  dat = simulate_data(df.m, df.n, df.r, [df.tau1, df.tau2], params)
  A = zeros(1, 2)
  A[1, :] = [1., -1.]
  
  # Precompile so the timing of the first opt isn't messed up
  try
    vmod = varcompmodel(@formula(y ~ 1 + x1 + x2 + (1 | subject) + (1 | item)), dat; control = VarCompControl(NewtonControl(), 0), A = A)
    success = 1
  catch e
    vmod = nothing
  end
  # Note: we do NOT return an error if the optimization fails, we capture this
  # and report as a result.

  # Now, run it with the different combinations of starting values
  # Use ANOVA, 0, and a "large" value.
  startingvalues = [
    [0.0, 0.0],
    [0.0, 10.0],
    [10.0, 10.0],
  ]

  startingvalue_labels = [
    "ANOVA",
    "[0.0, 0.0]",
    "[0.0, 10.0]",
    "[10.0, 10.0]",
  ]
  optruns = Vector{Union{varcomptest.VarCompModel, Nothing}}()
  # Do a run with ANOVA starting values
  vmod_default = try
    varcompmodel(@formula(y ~ 1 + x1 + x2 + (1 | subject) + (1 | item)), dat; control = VarCompControl(NewtonControl(), 0), A = A)
  catch e
    nothing
  end
  push!(optruns, vmod_default)
  # Now loop over the starting values
  for start in startingvalues
    vmod = try
      control = NewtonControl(startingvalues = start)
      varcompmodel(@formula(y ~ 1 + x1 + x2 + (1 | subject) + (1 | item)), dat; control = VarCompControl(control, 0), A = A)
    catch e
      nothing
    end
    push!(optruns, vmod)
  end
  
  # Report the results
  results_normal = DataFrame(
    itr = Int64[],
    numstephalves = Int64[],
    totalstephalves = Int64[],
    numuturns = Int64[],
    numgradient = Int64[],
    value = Float64[],
    gradnorm = Float64[],
    hessmineigen = Float64[],
    startingvalue = String[],
    time = Float64[],
    m = Int64[],
    n = Int64[],
    r = Float64[],
    tau1 = Float64[],
    tau2 = Float64[],
    tau1mle = Float64[],
    tau2mle = Float64[],
    simindex = Int64[],
  )

  results_conditional = copy(results_normal)

  function push_opt_result!(results::DataFrame, opt, startingvalue_label)
      if isnothing(opt)
          push!(results, (
              itr = -1,
              numstephalves = -1,
              totalstephalves = -1,
              numuturns = -1,
              numgradient = -1,
              value = -1.0,
              gradnorm = -1.0,
              hessmineigen = -1.0,
              startingvalue = startingvalue_label,
              time = -1.,
              m = df.m,
              n = df.n,
              r = df.r,
              tau1 = df.tau1,
              tau2 = df.tau2,
              tau1mle = -1.,
              tau2mle = -1., 
              simindex = df.sim,
          ))
      else
          push!(results, (
              itr = opt.itr,
              numstephalves = count(>(0), opt.history.step),
              totalstephalves = sum(opt.history.step),
              numuturns = count(identity, opt.history.uturn),
              numgradient = count(identity, opt.history.grad),
              value = opt.val,
              gradnorm = maximum(abs, opt.derivs.gradient),
              hessmineigen = minimum(abs.(eigvals(opt.derivs.Hessian))),
              startingvalue = startingvalue_label,
              time = opt.executiontime,
              m = df.m,
              n = df.n,
              r = df.r,
              tau1 = df.tau1,
              tau2 = df.tau2,
              tau1mle = opt.par[1],
              tau2mle = opt.par[2],
              simindex = df.sim,
          ))
      end

      return results
  end

  for i in eachindex(optruns)
      run = optruns[i]
      label = startingvalue_labels[i]

      if isnothing(run)
          push_opt_result!(results_normal, nothing, label)
          push_opt_result!(results_conditional, nothing, label)
      else
          push_opt_result!(results_normal, run.opt, label)
          push_opt_result!(results_conditional, run.optcond, label)
      end
  end
  

  # Return a named tuple
  (
    results_normal = results_normal,
    results_conditional = results_conditional,
  )
end

# Allocate output
Nsimtotal = size(allcombinations, 1)
results = Vector{NamedTuple{(:results_normal, :results_conditional), Tuple{Union{Nothing, DataFrame}, Union{Nothing, DataFrame}}}}(undef, Nsimtotal)

# JIT
_ = dosim(allcombinations[1, :])

function run_sims(allcombinations)
  N = size(allcombinations, 1)
  results = Vector{Any}(undef, N)

  # Shared atomic counter
  counter = Atomic{Int}(1)

  @sync for _ in 1:nthreads()
      @spawn begin
          while true
              # Atomically claim a task index
              b = atomic_add!(counter, 1)
              if b > N
                  return
              end

              results[b] = dosim(allcombinations[b, :])
          end
      end
  end

  return results
end

results = run_sims(allcombinations)



# Write to CSV
results_normal = vcat([item.results_normal for item in results]...)
results_conditional = vcat([item.results_conditional for item in results]...)
results_normal.type = fill("normal", nrow(results_normal))
results_conditional.type = fill("conditional", nrow(results_conditional))
results_total = vcat(results_normal, results_conditional)

# Write to csv, make plots in R
CSV.write("$savename.csv", results_total)

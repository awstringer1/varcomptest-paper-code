### Simulations: Crossed Random Effects, zero component ###

## Paths ----

resultspath = homedir() * "/"
if !isdir(resultspath)
  throw(SystemError("Could not access directory: $resultspath"))
end
# Parse command line arguments
if length(ARGS) == 0
  nsim = 100
  DATE = "20260812"
  SIMVERSION = "1"
  bootsamps = 1000
  sampling = 2
  JOBID = "1"
  SIMIDX = "1"
  TASKID = "1"
  SAMPLINGTYPE = "unbalanced"
else
  nsim = parse(Int64, ARGS[1])
  DATE = ARGS[2]
  SIMVERSION = ARGS[3]
  bootsamps = parse(Int64, ARGS[4])
  sampling = parse(Int64, ARGS[5])
  JOBID = ARGS[6]
  SIMIDX = ARGS[7]
  TASKID = ARGS[8]
  SAMPLINGTYPE = sampling == 1 ? "balanced" : sampling == 2 ? "unbalanced" : "unknown"
  if SAMPLINGTYPE == "unknown"
    throw("Unknown sampling type $sampling")
  end
end
simname = "sims-crossed-zerocomponents-$DATE-v$SIMVERSION-s$SAMPLINGTYPE-job$JOBID-sim$SIMIDX-task$TASKID"
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
tau2 = [0., .01, .02, .05, .1]

params = Dict(
  :beta => [0., 1., -2.],
  :sigma => 1.,
  :B => bootsamps
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

## Simulation metadata ----

allcombinations = repeat(DataFrame(
    [(; m=i[1], n=i[2], r=i[3], tau1=i[4], tau2=i[5])
     for i in Iterators.product(m, n, r, tau1, tau2)]
), nsim)
allcombinations.sim = 1:size(allcombinations, 1)

## Perform simulations ----

dosim = function(df::DataFrameRow{DataFrame, DataFrames.Index}; params::Dict = params) 
  dat = simulate_data(df.m, df.n, df.r, [df.tau1, df.tau2], params)
  # Three hypotheses -> three model fits
  Afull = Matrix{Float64}(I, 2, 2)
  A1 = zeros(1, 2)
  A1[1, :] = [1., 0.]
  A2 = zeros(1, 2)
  A2[1, :] = [0., 1.]
  control = VarCompControl(NewtonControl(verbose=false), params[:B])
  success = -1
  vmodfull = nothing
  vmod1 = nothing
  vmod2 = nothing
  tm = @elapsed begin
    try
      vmodfull = varcompmodel(@formula(y ~ 1 + x1 + x2 + (1 | subject) + (1 | item)), dat; control = control, A = Afull)
      vmod1 = varcompmodel(@formula(y ~ 1 + x1 + x2 + (1 | subject) + (1 | item)), dat; control = control, A = A1)
      vmod2 = varcompmodel(@formula(y ~ 1 + x1 + x2 + (1 | subject) + (1 | item)), dat; control = control, A = A2)
      success = 1
    catch e
      vmodfull = nothing
      vmod1 = nothing
      vmod2 = nothing
    end
  end

  if any(isnothing, (vmodfull, vmod1, vmod2))
    return (
      success = -1,
      pvals = nothing,
      time = tm,
      simparams = NamedTuple(df)
    )
  end

  pvals = (
    m = df.m,
    n = df.n,
    r = df.r,
    tau1 = df.tau1,
    tau2 = df.tau2,
    tau1est = vmodfull.vr.tau[1],
    tau2est = vmodfull.vr.tau[2],
    tau1bootsd = std(vmodfull.bootresults.mle[:, 1]),
    tau2bootsd = std(vmodfull.bootresults.mle[:, 2]),
    tau1q025 = quantile(vmodfull.bootresults.mle[:, 1], 0.025),
    tau1q975 = quantile(vmodfull.bootresults.mle[:, 1], 0.975),
    tau2q025 = quantile(vmodfull.bootresults.mle[:, 2], 0.025),
    tau2q975 = quantile(vmodfull.bootresults.mle[:, 2], 0.975),
    pvalfull = vmodfull.bootresults.pval,
    pvaltau1 = vmod1.bootresults.pval,
    pvaltau2 = vmod2.bootresults.pval,    
    pvalf1 = vmodfull.vr.ftable[1,5],
    pvalf2 = vmodfull.vr.ftable[2,5],
    pvalfboth = vmodfull.vr.ftable[4,5]
  )
  (
    success = success,
    pvals = pvals,
    time = tm,
    simparams = NamedTuple(df),
    params = params
  )
end

# Allocate output
Nsimtotal = size(allcombinations, 1)
results = Vector{NamedTuple{(:success, :pvals, :time, :simparams, :params), Tuple{Int64, Union{Nothing, NamedTuple}, Float64, NamedTuple, Dict}}}(undef, Nsimtotal)

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
Nsimunclean = length(results)
resultsclean = filter(res -> res.success == 1, results)
Nsimtotal = length(resultsclean)
pvals = DataFrame([item.pvals for item in resultsclean])
pvals.difftau = pvals.tau1 - pvals.tau2
pvals.absdifftau = abs.(pvals.tau1 - pvals.tau2)
d = 2
r = 1

# Write to csv, make plots in R
CSV.write("$savename.csv", pvals)


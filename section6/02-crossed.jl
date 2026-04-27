### Simulations: Crossed Random Effects ###

## Paths ----

resultspath = homedir()
if !isdir(resultspath)
  throw(SystemError("Could not access directory: $resultspath"))
end

nsim = parse(Int64, ARGS[1])
DATE = ARGS[2]
SIMVERSION = ARGS[3]
bootsamps = parse(Int64, ARGS[4])
JOBID = ARGS[5]
SIMIDX = ARGS[6]
TASKID = ARGS[7]

simname = "sims-crossed-$DATE-v$SIMVERSION-job$JOBID-sim$SIMIDX-task$TASKID"
savename = resultspath*simname
println("simname = ", simname)
## Packages ----

using 
  varcomptest,
  Base.Threads,
  DataFrames,
  MixedModels,
  Distributions,
  CategoricalArrays,
  StatsPlots,
  CSV,
  Random


## Parameters ----

m = [20, 50, 100]
n = [20, 50, 100]
r = [-1., 0., .5] # This is now the sampling correlation not the reps
tau1 = [0., .01, .02, .05, .1]
tau2 = [0., .01, .02, .05, .1]

params = Dict(
  :beta => [0., 1., -2.],
  :sigma => 1.,
  :B => bootsamps,
  :sampling => -1
)
# sampling = -1: m x n fully crossed, N = mn
# sampling = r \in [0,1): sample the indices from a correlated discrete uniform without replacement with correlation r

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


## Simulation metadata ----

allcombinations = repeat(DataFrame(
    [(; m=i[1], n=i[2], r=i[3], tau1=i[4], tau2=i[5])
     for i in Iterators.product(m, n, r, tau1, tau2)]
), nsim)
allcombinations.sim = 1:size(allcombinations, 1)

dosim = function(df::DataFrameRow{DataFrame, DataFrames.Index}; params::Dict = params) 
  dat = simulate_data(df.m, df.n, df.r, [df.tau1, df.tau2], params)
  A = zeros(1, 2)
  A[1, :] = [1., -1.]
  control = VarCompControl(NewtonControl(), params[:B])
  success = -1
  vmod = nothing
  tm = @elapsed begin
    try
      vmod = varcompmodel(@formula(y ~ 1 + x1 + x2 + (1 | subject) + (1 | item)), dat; control = control, A = A)
      success = 1
    catch e
      vmod = nothing
    end
  end

  if vmod == nothing
    return (
      success = -1,
      pvals = nothing,
      time = tm,
      simparams = NamedTuple(df)
    )
  end

  # crossed model, 2 random effects, 1 hypothesis
  d = 2
  r = 1

  pvals = (
    m = df.m,
    n = df.n,
    r = df.r,
    tau1 = df.tau1,
    tau2 = df.tau2,
    tau1est = vmod.vr.tau[1],
    tau2est = vmod.vr.tau[2],
    tau1bootsd = std(vmod.bootresults.mle[:, 1]),
    tau2bootsd = std(vmod.bootresults.mle[:, 2]),
    tau1q025 = quantile(vmod.bootresults.mle[:, 1], 0.025),
    tau1q975 = quantile(vmod.bootresults.mle[:, 1], 0.975),
    tau2q025 = quantile(vmod.bootresults.mle[:, 2], 0.025),
    tau2q975 = quantile(vmod.bootresults.mle[:, 2], 0.975),
    pvalboot = vmod.bootresults.pval,
    pvaloneside = vmod.bootresults.pvaloneside,
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

# Run simulations

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
# Extract chisq and bootstrap p-values
pvals = DataFrame([item.pvals for item in resultsclean])
pvals.difftau = pvals.tau1 - pvals.tau2
pvals.absdifftau = abs.(pvals.tau1 - pvals.tau2)
# Write to csv, make plots in R
CSV.write("$savename.csv", pvals)
println("Wrote results to file: $savename.csv")
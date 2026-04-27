### Simulations: Nested Random Effects ###

## Paths ----

resultspath = homedir()
if !isdir(resultspath)
  throw(SystemError("Could not access directory: $resultspath"))
end
# Parse command line arguments
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
simname = "sims-nested-$DATE-v$SIMVERSION-s$SAMPLINGTYPE-job$JOBID-sim$SIMIDX-task$TASKID"
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
  Base.Threads

nt = Threads.nthreads();

## Parameters ----

m = [20, 50, 100]
n = [2, 5, 10]
r = [2, 5, 10]
tau1 = [0., .01, .02, .05, .1]
tau2 = [0., .01, .02, .05, .1]
params = Dict(
  :beta => [0., 1., -2.],
  :sigma => 1.,
  :B => bootsamps,
  :sampling => sampling
)

## Functions ----

simulate_data = function(m::Int64, n::Int64, r::Int64, tau::Vector{Float64}, params::Dict)
  N = n * m * r
  p = length(params[:beta])

  if params[:sampling] == 1    
    nj = repeat(n:n, inner = m)
    rj = [repeat(r:r, inner = njj) for njj in nj]
  elseif params[:sampling] == 2
    nj = rand(2:(2 * n - 2), m)
    rj = [rand(2:(2 * r - 2), njj) for njj in nj]
  elseif params[:sampling] == 3
    nj = [2 + m * (n - 2), repeat(2:2, inner = m - 1)...]
    rj = [repeat(r:r, inner = njj) for njj in nj]
  else
    ss = params[:sampling]
    throw("Unknown sampling strategy $ss, options are 1, 2, or 3")
  end

  block = vcat([vcat([fill(i, rj[i][j]) for j in 1:nj[i]]...) for i in 1:m]...)
  wholeplot = vcat([vcat([fill(j, rj[i][j]) for j in 1:nj[i]]...) for i in 1:m]...)
  blockwholeplot = string.(block) .* ":" .* string.(wholeplot)
  N = length(blockwholeplot)

  Nb = length(unique(block))
  Nw = length(unique(blockwholeplot))
  u1 = rand(Normal(0, params[:sigma] * sqrt(tau[1])), Nb);
  u2 = rand(Normal(0, params[:sigma] * sqrt(tau[2])), Nw);

  dat = DataFrame(
    y = zeros(N),
    block = block,
    wholeplot = wholeplot,
    blockwholeplot = blockwholeplot
  );
  
  X = hcat(ones(N), rand(Uniform(0, 1), N, p-1));
  mu = X * params[:beta];

  for j in 1:(p-1)
    dat[!, "x$j"] = X[:, j + 1]
  end;
  
  blockidx = 0
  plotidx = 0
  totalidx = 0
  for i in 1:Nb
    blockidx += 1
    for j in 1:nj[i]
      plotidx += 1
      for k in 1:rj[i][j]
        totalidx += 1
        dat.y[totalidx] = rand(Normal(mu[totalidx] + u1[blockidx] + u2[plotidx], params[:sigma]))
      end
    end
  end;
  
  # convert to categorical
  dat.block = categorical(dat.block, ordered = true, levels = unique(dat.block));
  dat.wholeplot = categorical(dat.wholeplot, ordered = true, levels = unique(dat.wholeplot));
  dat.blockwholeplot = categorical(dat.blockwholeplot, ordered = false, levels = unique(dat.blockwholeplot));

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
  A = zeros(1, 2)
  A[1, :] = [1., -1.]
  control = VarCompControl(NewtonControl(), params[:B])
  success = -1
  vmod = nothing
  tm = @elapsed begin
    try
      vmod = varcompmodel(@formula(y ~ 1 + x1 + x2 + (1 | block) + (1 | blockwholeplot)), dat; control = control, A = A)
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

  # Nested model, 2 random effects, 1 hypothesis
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
    pvaloneside = vmod.bootresults.pvaloneside
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


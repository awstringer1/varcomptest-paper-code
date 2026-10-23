### Simulations: Optimization, Nested Random Effects ###

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
  bootsamps = parse(Int64, ARGS[4])
  sampling = parse(Int64, ARGS[5])
  JOBID = ARGS[6]
  SIMIDX = ARGS[7]
  TASKID = ARGS[8]
end
SAMPLINGTYPE = sampling == 1 ? "balanced" : sampling == 2 ? "unbalanced" : sampling == 3 ? "veryunbalanced" : "unknown"
if SAMPLINGTYPE == "unknown"
  throw("Unknown sampling type $sampling")
end
simname = "sims-nested-opt-$DATE-v$SIMVERSION-s$SAMPLINGTYPE-job$JOBID-sim$SIMIDX-task$TASKID"
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
n = [2, 5, 10]
r = [2, 5, 10]
tau1 = [0., .01, .02, .05, .1]
taudiff = [0., .01, .05, .1, 1., 2.] # tau2 = tau1 + taudiff
params = Dict(
  :beta => [0., 1., -2.],
  :sigma => 1.,
  :B => 0, # NOT bootstrapping here, just optimizing
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
        # println("i = $i, block = $blockidx, plot = $plotidx, total = $totalidx")
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

# tau2 = tau1 + taudiff
allcombinations = repeat(DataFrame(
    [(; m=i[1], n=i[2], r=i[3], tau1=i[4], tau2=i[4] + i[5])
     for i in Iterators.product(m, n, r, tau1, taudiff)]
), nsim)
allcombinations.sim = 1:size(allcombinations, 1)

dosim = function(df::DataFrameRow{DataFrame, DataFrames.Index}; params::Dict = params) 
  dat = simulate_data(df.m, df.n, df.r, [df.tau1, df.tau2], params)
  A = zeros(1, 2)
  A[1, :] = [1., -1.]
  
  # Precompile so the timing of the first opt isn't messed up
  try
    vmod = varcompmodel(@formula(y ~ 1 + x1 + x2 + (1 | block) + (1 | blockwholeplot)), dat; control = VarCompControl(NewtonControl(), 0), A = A)
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
    varcompmodel(@formula(y ~ 1 + x1 + x2 + (1 | block) + (1 | blockwholeplot)), dat; control = VarCompControl(NewtonControl(), 0), A = A)
  catch e
    nothing
  end
  push!(optruns, vmod_default)
  # Now loop over the starting values
  for start in startingvalues
    vmod = try
      control = NewtonControl(startingvalues = start)
      varcompmodel(@formula(y ~ 1 + x1 + x2 + (1 | block) + (1 | blockwholeplot)), dat; control = VarCompControl(control, 0), A = A)
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
    r = Int64[],
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

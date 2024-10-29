
function size_structured_dde_data( yrs, aulab, fishery_data, model_variation )

  Yyrs = floor.(Int, fishery_data["Y"].yrs)
  Y = fishery_data["Y"][∈(yrs).(Yyrs), :]
  
  removalsyrs = floor.(Int, fishery_data["L"].yrs)
  removals = fishery_data["L"][∈(yrs).(removalsyrs), :]  # in numbers (not mass)
  
  MW = fishery_data["M0_W"][∈(yrs).(fishery_data["M0_W"].mw_yrs), :]
  MW.yrs = MW.mw_yrs
    
  Kmu = [5.5, 60.0, 1.25]   # 5.0, 60.0, 1.25 
  
      #=
          # alternatively, if running manually:
          # can run R-code that creates local RData file with required data
          # run in R externally or from within julia or ..
  
          # from within julia
  
          using RCall
          # type $ in Julia's command prompt starts an R session.
          # .. run below
          # type <backspace> to escape back to julia
  
          source( file.path( code_root, "bio_startup.R" )  )
          require(bio.snowcrab)   # loadfunctions("bio.snowcrab")
          fishery_model_data_inputs( year.assessment=year.assessment, type="size_structured_numerical_dynamics",  for_julia=TRUE, time_resolution=1/12)
  
          # then back in Julia, fetch data into julia's workspace (replace fndat with the  filenane printed above )
          @rget Y
          @rget Kmu
          @rget removals
          @rget ty
  
          # example line plots
          plot(  Y[:,:yrs], Y[:,:cfasouth_M0] )
       
      =#
  
  
  nT = length(yrs)
  nP = 5  # number of predictions into future (with no fishing)
  nM = nP + nT  # total number of prediction years
  
  nS = 6  # no. state variables
  
  # BLY =[8]
  BLY = [7,8,9] #  birth lag years, centered on 8 yrs offset
  
  nB = length(BLY)
  
  # "survey index"
  statevars = [
    Symbol("$aulab","_M0"),
    Symbol("$aulab","_M1"),
    Symbol("$aulab","_M2"),
    Symbol("$aulab","_M3"),
    Symbol("$aulab","_M4"),
    Symbol("$aulab","_f_mat")
  ]
  
  S = Matrix(Y[:, statevars ])
  
  # scale index where required
  Smean = [mean(skipmissing(S)) for i in 1:nS ]
  Sstd = [std( skipmissing(S)) for i in 1:nS ]
  Smin = [minimum(skipmissing(S[:,i])) for i in 1:nS ]
  Smax = [maximum(skipmissing(S[:,i])) for i in 1:nS ]
  Srange = Smax .- Smin 
  
  SminFraction = Smin ./ Srange  # used as informative prior mean in some runs
  
  # CV of Y
  statevars_sd = [
    Symbol("$aulab", "_sd", "_M0"),
    Symbol("$aulab", "_sd", "_M1"),
    Symbol("$aulab", "_sd", "_M2"),
    Symbol("$aulab", "_sd", "_M3"),
    Symbol("$aulab", "_sd", "_M4"),
    Symbol("$aulab", "_sd", "_f_mat")
  ]
  
  Ssd = Matrix(Y[:, statevars_sd ])
  
  Scv = Ssd ./ S  
  
  # scale index to min-max
  if occursin( r"unnormalized", model_variation )
    # do nothing (no scaling)
  elseif occursin( r"scale_center", model_variation ) 
    for i in 1:nS
      S[:,i] = (S[:,i] .- Smean[i] ) ./ Sstd[i]    # scale to std and center to 0 
    end
  else 
    # default is to normalize (min, max) to (0,1)
    for i in 1:nS
      S[:,i] = (S[:,i] .- Smin[i] ) ./ Srange[i]   # range from 0=min to 1=max
    end
  end
   
  i =findall(x -> !ismissing(x) && x < 1.0e-9, S)
  if length(i) > 0 
    S[i] .= 1.0e-9
  end
  
  
  # deal with missing CV's
  for i in 1:nS
    u = findall( x-> ismissing(x), Scv[:,i] )  
    if length(u) > 0
      Scv[u,i] .= S[u,i] # ie. poisson
    end
  end
  for i in 1:nS
    u = findall( x-> ismissing(x), Scv[:,i] )  
    if length(u) > 0
      Scv[u,i] .= 0.5 # ie. no information ( u is in interval 0,1 .. 0.5 covers it nicely )
    end
  end
  
  logScv = log.(Scv) # on log scale to reduce further computations
  
  
  # interpolating function for mean weight (of fb only)
  mwspline = extrapolate( interpolate( MW[:,Symbol("mw_", "$aulab") ], (BSpline(Linear()) ) ),  Interpolations.Flat() )
  mw = Interpolations.scale(mwspline, yrs )
  
  scale_factor = mw(yrs) / (1000 *1000 ) # convert numbers to kt biomass , also used in plots
   
  # convert to (biomass kt to number)
  
  # id index
  ki = aulab == "cfanorth" ? 1 :
       aulab == "cfasouth" ? 2 :
       aulab == "cfa4x"    ? 3 :
       0  # default
  
  kmu  =  Kmu[ki] / mean(scale_factor)
  logkmu = log(kmu)
  
  smallnumber = 1.0 / (kmu * 10.0) # floating point value of sufficient to assume 0 valued
       
  dt = (0.02, 0.02, 0.02)[ki]    # resolution of time (fraction of year)  (i.e. days = dt*365; weeks=dt*52); dt=0.02 is ~ weekly
  
  # spin up time of ~ 1 cycle prior to start of dymamics and project nP years into the future
  tspan = (minimum(yrs) - 11.1, maximum(yrs) + nP + 1.1 )
  
  survey_time = discretize_decimal( (Y[:,:yrs] .+ 0.999), dt)     # add 0.999 to make it the end of the year
  
  Si = findall( x-> !ismissing(x), vec(sum(S, dims=2)))  # compute data likelihoods only when data exist ... to speed up comps
  nSI = length(Si)
   
  PM = (
      nS = nS, 
      nSI = nSI,
      nB = nB,
      nG = 4,  # n transition moults .. growth
      nT = length(yrs),
      nP = 5,  # number of predictions into future (with no fishing)
      nM = nP + nT,  # total number of prediction years
      logkmu = (logkmu, 0.1),
      logScv = (logScv, 0.1),
      kmu = kmu,
      yrs = yrs,
      b = ( log(1.0), 0.1),
      d =  ( log( exp(0.25)-1.0 ), 0.1 ),
      d2 = ( log( exp(0.75)-1.0 ), 0.1 ),
      v =  ( log( exp(0.90)-1.0 ), 0.1 ),
      q1 = ( 0.5, 0.1 ),  # assume about 50% overly optimistic estimation from CARSTM (due to assumption that each areal unit is homogenous)
      q0 = ( SminFraction ./ 2.0, 0.1 ),  # lower detection limit of survey/CARSM
      BLY = BLY,
      Si = Si,
      S = S,
      data = S[Si,:],
      # datavector = vec(S[Si,:]),  # for MVN .. no advantage
      Stime = survey_time[Si],
      eps=1e-9
  ) 


  # this only adds habitat space  ... predation is also a useful one ..
  predtime =  discretize_decimal( 9.0/12.0, dt )
  prediction_time = floor.( vcat( collect(minimum(yrs) : (maximum(yrs)+nP) ) )  ) .+ predtime   # sept
  
  #  sa to fraction
  external_forcing =  reshape( [
     Y[:,Symbol("H", "$aulab","_M0")]  / maximum( Y[:,Symbol("H", "$aulab","_M0")] )
     Y[:,Symbol("H", "$aulab","_M1")]  / maximum( Y[:,Symbol("H", "$aulab","_M1")] )
     Y[:,Symbol("H", "$aulab","_M2")]  / maximum( Y[:,Symbol("H", "$aulab","_M2")] )
     Y[:,Symbol("H", "$aulab","_M3")]  / maximum( Y[:,Symbol("H", "$aulab","_M3")] )
     Y[:,Symbol("H", "$aulab","_M4")]  / maximum( Y[:,Symbol("H", "$aulab","_M4")] )
     Y[:,Symbol("H", "$aulab","_f_mat")]  / maximum( Y[:,Symbol("H", "$aulab","_f_mat")] )
    ], nT, nS )
  
  
  efc = extrapolate( interpolate( external_forcing, (BSpline(Linear()), NoInterp()) ), Interpolations.Flat() )
  hsa = Interpolations.scale(efc, yrs .+ predtime, 1:nS )
  
  ys = ( "yrs", "yrs", "yrs_4x")[ki]
  
  # fishing pattern seasonal over past 5 yrs (function) .. used for projections
  fish_time_max = discretize_decimal( maximum(removals[:,:ts]), dt)  + dt 
  fp_time = fish_time_max:dt:maximum(prediction_time)
  fish_time_project = discretize_decimal( collect(fp_time), dt )     
  
  # fish_year is fishery "year" 
  fsh = DataFrame(
    :fish_time => discretize_decimal( removals[:,:ts] , dt) ,
    :fish_year => discretize_decimal( removals[:,Symbol(ys)], dt),
    :removed => removals[:,Symbol("$aulab")]  #number
  )
  
  # keep nonzero elements
  fsh = fsh[findall( x-> x>0, fsh.removed),:]
  
  # aggregate if required
  afsh = DataFrame( combine( groupby(fsh, [:fish_time]), :removed => sum, :fish_year => unique) )
  
  # make available as a global variable
  fish_time = afsh.fish_time
  fish_year = afsh.fish_year_unique
  removed = afsh.removed_sum
   
  # model-specifics functions and data
   
  # DiffEq-model setup
  
  if model_variation=="size_structured_dde_normalized" 

    function affect_fishing!(integrator)
      i = findall(t -> t == integrator.t, fish_time)[1]
      integrator.u[1] -=  removed[ i ] / integrator.p[3][1]  # p[3] ==K divide by K[1]  .. keep unscaled to estimate magnitude of other components
    end

  elseif model_variation=="size_structured_dde_unnormalized"

    function affect_fishing!(integrator)
      i = findall(t -> t == integrator.t, fish_time)[1]
      integrator.u[1] -=  removed[ i[1] ]  # sol on same scale
    end

  else 
    error("model_variation not found")
  end
   
  # callbacks for external perturbations to the system (deterministic fishing without error)
  cb = PresetTimeCallback( fish_time, affect_fishing! )
  # alternative formulation:
  # cb = CallbackSet( PresetTimeCallback( fish_time, affect_fishing! ), PositiveDomain() );
   
  # generics to bootstrap the process; initial conditions .. "pristine", unfished condition ... 
  # however, natural mortality might have been different so be careful here
  u0 = 
    model_variation=="size_structured_dde_unnormalized" ? ones(nS) .* kmu .* 1.0 :
    model_variation=="size_structured_dde_normalized"   ? ones(nS) .* 1.0 :
    1.0
  
  # history function (prior to start)  defaults to values of kmu / 1 before t0;  
  # h(p, t; idxs=nothing) = typeof(idxs) <: Number ? u0 : u0 
  h(p, t) = u0 
  
  tau = BLY  # delay resolution
  p = dde_parameters( PM ) # dummy values needed to bootstrap DifferentialEquations/Turing initialization
  prob = DDEProblem{true}( size_structured_dde!, u0, h, tspan, p, constant_lags=tau  )  #  , neutral=true create container for problem definition 
  
  # choose DiffEq solver:
  # stiff solvers: Rodas4()  ; Rosenbrock23()
  # diffeq_solver = MethodOfSteps(Rosenbrock23()) # slow but good for stiff problems
  # diffeq_solver = MethodOfSteps(AutoTsit5(Rodas5()))
  # diffeq_solver = MethodOfSteps(AutoTsit5(KenCarp47()))
  # diffeq_solver = MethodOfSteps(KenCarp47())  
  # diffeq_solver = MethodOfSteps(AutoTsit5())
  diffeq_solver = MethodOfSteps(Tsit5())   # faster
  # diffeq_solver = MethodOfSteps(Vern6() )  # better resolution of error
  # diffeq_solver = MethodOfSteps(Rodas5())  # safer
  # diffeq_solver =  LSODA.lsoda()  #ode
  # diffeq_solver = MethodOfSteps(QNDF())  # stiff
  # diffeq_solver = MethodOfSteps(TRBDF2())  # stiff
   
  
  solver_params = (
    prob=prob,
    abstol = 1.0e-12,    
    reltol = 1.0e-9,
    maxiters=1.0e6, 
    dt = dt,
    saveat = collect(tspan[1]:dt:tspan[2]),
    cb = cb,
    tspan = tspan,
    h=h,
    hsa=hsa,
    solver = diffeq_solver
  )
  
  #=  lognormal rates:
      # Can go higher than 100% .. as the size-based categories are imperfect and there is also input from other groups
      log( exp(0.1)-1.0 ) = log(0.1052 ) = -2.252  .. 10% mortality (mode)
      log( exp(0.2)-1.0 ) = log(0.2214 ) = -1.508  .. 20% mortality (mode)
      log( exp(0.3)-1.0 ) = log(0.3499 ) = -1.050  .. 30% 
      log( exp(0.4)-1.0 ) = log(0.4918 ) = -0.7096 .. 40%
      log( exp(0.5)-1.0 ) = log(0.6487 ) = -0.4328 .. 50%
    
      # Can go higher than 100% .. as the size-based categories are imperfect and there is also input from other groups
      log( exp(0.90)-1.0) = log(1.46)  = 0.3782  .. ~90% (moult) transition rate per year (mode)
      log( exp(0.95)-1.0) = log(1.586) = 0.461   .. ~95% (moult) transition rate per year (mode) 
      log( exp(0.99)-1.0) = log(1.6912 = 0.5254   .. ~99% (moult) transition rate per year (mode) 
  =#
  
  # pl = plot(x->pdf(LogNormal(log(10), 1.0), x), xlim=(0,10)) #
  
  # choose model and over-rides if any
   
  # the following tweak Lognormal priors by area  
  
    
  if model_variation=="size_structured_dde_normalized"
    
    fmod = size_structured_dde_turing( PM, solver_params )

  elseif  model_variation=="size_structured_dde_unnormalized"

    print( "warning: model needs some updating .. do not use until it is checked" )

    fmod = size_structured_dde_turing( S, kmu, tspan, prob, nS )
  end

  return PM, fmod, solver_params 
  
end


 


function size_structured_dde!( du, u, h, p, t )
  # here u, du are actual numbers .. not normalized by K due to use of callbacks
  b5, b6, K, d, d2, v  = p  # hsa is global 
    begin
      ux =  max.(solver_params.abstol, u)
      trg = [h(p, t-j)[6] for j in PM.BLY ]
      tr = v  .* h(p, t-1)[2:5] 
      dh = d2 .* ux ./ solver_params.hsa(t, 1:6) 
      dr = d  .* ux .+  dh .* ux
      du[1] = tr[1] * K[2] / K[1]              - dr[1]       # note:
      du[2] = tr[2] * K[3] / K[2]     - tr[1]  - dr[2]
      du[3] = tr[3] * K[4] / K[3]     - tr[2]  - dr[3]
      du[4] = tr[4] * K[5] / K[4]     - tr[3]  - dr[4]
      du[5] = b5' * trg * K[6] / K[5] - tr[4]  - dr[5]
      du[6] = b6' * trg                        - dr[6]      # fem mat simple logistic with lag tau and density dep on present numbers
    end
end


function dde_parameters(PM)
    # these are dummy initial values .. just to get things started
    b5= repeat( [1.0], PM.nB )
    b6= repeat( [1.0], PM.nB )
    K = repeat( [PM.kmu], PM.nS )  
    d=[0.15, 0.11, 0.14, 0.17, 0.16, 0.19];
    d2 = [0.4, 0.4, 0.4, 0.4, 0.4, 0.4] 
    v=[0.65, 0.68, 0.61, 0.79];
    params = ( b5, b6, K, d, d2, v)
    return params
end

Turing.@model function size_structured_dde_turing( PM, solver_params )
 
  K ~ filldist( LogNormal( PM.logkmu[1], PM.logkmu[2]), PM.nS )  # kmu already on log scale # is max of a multiyear group , serves as upper bound for all
  q1 ~ filldist( Normal( PM.q1[1], PM.q1[2] ), PM.nS )

  q0 ~ arraydist([Normal( PM.q0[1][i], PM.q0[2] ) for i in 1:PM.nS])  # informative prior on relative height 
  model_sd ~  arraydist(LogNormal.( PM.logScv[1], PM.logScv[2]) ) # working: β(0.1, 10.0);  plot(x->pdf(β(0.3,12), x), xlim=(0,1)) # uniform 
  
  # minor speed improvement (if any) with LazyArray
  # q0 ~ arraydist(LazyArray(@~ Normal.( PM.q0[1], PM.q0[2] ) ) )  # informative prior on relative height 
  # model_sd ~  arraydist(LazyArray(@~ LogNormal.( PM.logScv[1], PM.logScv[2]) ) ) # working: β(0.1, 10.0);  plot(x->pdf(β(0.3,12), x), xlim=(0,1)) # uniform 
  
  # BroadcastArray minor (if any) speed improvement
  # q0 ~ arraydist(BroadcastArray(Normal, PM.q0[1], PM.q0[2]  ) )  # informative prior on relative height 
  # model_sd ~  arraydist(BroadcastArray(LogNormal, PM.logScv[1], PM.logScv[2] ) ) # working: β(0.1, 10.0);  plot(x->pdf(β(0.3,12), x), xlim=(0,1)) # uniform 

  b5 ~  filldist( LogNormal( PM.b[1],  PM.b[2] ), PM.nB )   # centered on 1; plot(x->pdf(LogNormal(log(10), 1.0), x), xlim=(0,10)) # mode of 5
  b6 ~  filldist( LogNormal( PM.b[1],  PM.b[2] ), PM.nB )   # centered on 1; plot(x->pdf(LogNormal(log(10), 1.0), x), xlim=(0,10)) # mode of 5
  d ~   filldist( LogNormal( PM.d[1],  PM.d[2] ), PM.nS ) # plot(x->pdf(LogNormal(0.2, 1.0), x), xlim=(0, 2)) 
  d2 ~  filldist( LogNormal( PM.d2[1], PM.d2[2]), PM.nS ) # plot(x->pdf(LogNormal(-0.7096, 0.25 ), x), xlim=(0, 2)) 
  v ~   filldist( LogNormal( PM.v[1],  PM.v[2] ), PM.nG ) # transition rates # plot(x->pdf(LogNormal( 0.3782, 0.5 ), x), xlim=(0,1))  
  u0 ~  filldist( Beta(5, 2), PM.nS )  # plot(x->pdf(Beta(5, 2), x), xlim=(0,1)) # uniform 
     
    #likelihood
    # Like this, samples in 50 s
    # y .~ BernoulliLogit.(α .+ X * β)

    # Like this, samples in 45 s
    # y ~ arraydist(BernoulliLogit.(α .+ X * β))
 
    # Like this, samples in 15 s
    # y ~ arraydist(LazyArray(@~ BernoulliLogit.(α .+ X * β)))

  # process model
  prob_new = remake( solver_params.prob; u0=u0, h=solver_params.h, tspan=solver_params.tspan, p=( b5, b6, K, d, d2, v ) )

  msol = solve(
      prob_new,
      solver_params.solver, 
      callback=solver_params.cb,
      abstol=solver_params.abstol, 
      # reltol=solver_params.reltol, 
      # maxiters=solver_params.maxiters,
      # isoutofdomain=(y,p,t)->any(x -> x < solver_params.abstol, y),  # permit exceeding K, only check lower
      saveat=solver_params.saveat
      # saveat=solver_params.dt
      # alg_hints=[:stiff]
      # ,
      # saveat=solver_params.saveat
  )

  # @show msol.retcode
  if !SciMLBase.successful_retcode(msol)  
    Turing.@addlogprob! -Inf
    return nothing
  end

  ii = indexin(PM.Stime,  msol.t)
  ii = ii[ .!(isnothing.(ii)) ]


  if length(ii) != PM.nSI
    Turing.@addlogprob! -Inf
    return nothing
  end
  
  # map S <=> m  where S = observation index on unit scale; m = latent, scaled abundance on unit scale
  # observation model: S = (m - q0)/ q1   <=>   m = S * q1 + q0  
  A = ( Array(msol)[:, ii] .- q0 ) ./ q1

  if any( x -> !isfinite(x), A)  
    Turing.@addlogprob! -Inf
    return nothing
  end
 
    
  # ia = findall( x -> x < solver_params.abstol, A )
  # if length(ia) > 0
  #   Turing.@addlogprob! -Inf
  #   A[ia] .= solver_params.abstol
  #   return nothing
  # end

  sigma = view( model_sd, PM.Si, : )
  
  # rand.( Normal.( A', sigma )   )
  @. PM.data ~ Normal( A', sigma )  
  # PM.data ~ arraydist( LazyArray(@~ Normal.(A', sigma) ) )
  # PM.data ~ arraydist( BroadcastArray(Normal, A', sigma ) )
   
 
  # equivalent representations:
  # PM.datavector ~ MvNormal( vec(A'), Diagonal(vec(sigma).^2.0)  )  # no real speed gain using MVN .. JC: Feb 2023

end
 
# ------------------------------
  

 

# -----------
  


function fishery_model_predictions_timeseries( num; prediction_time, plot_k )
  gk = num[:,plot_k,:,1]
  pl = plot()
  pl = plot!(pl, prediction_time, gk;  alpha=0.02, color=:lightslateblue)
  pl = plot!(pl, prediction_time, mean(gk, dims=2);  alpha=0.8, color=:darkslateblue, lw=4)
  pl = plot!(pl; legend=false )
  pl = plot!(pl; ylim=(0, maximum(gk)*1.01 ) )

  S_K_sims = abundance_from_index( S, res; k=plot_k )
  S_K = mean(S_K_sims, dims=2)  # average by year

  pl = plot!(pl, survey_time, S_K, color=:gray, lw=2 )
  pl = scatter!(pl, survey_time, S_K, markersize=4, color=:grey)
  pl = plot!(pl; legend=false )
  return (gk, pl)
end

# -----------





function fishing_mortality_instantaneous( removed, abundance )
  -log(  1.0 - (removed  / abundance)  )  ;
end



function removals_aggregate( removed, fish_year )
  landings_aggregated = DataFrame( yr=floor.(fish_year), remNum = removed );
  landings_aggregated.rem = landings_aggregated.remNum .* mw(fish_time) 
  out = combine(groupby(landings_aggregated,:yr),[:rem ] .=> sum )
  stimes = DataFrame( yr=floor.(survey_time) )
  out = leftjoin(stimes, out, on=:yr)
  sort!(out, :yr)
  oo = findall(x->ismissing(x), out[:,:rem_sum])
  if length(oo) > 0
    out[ oo[1], :rem_sum ] = 0.0
  end
  return(out)
end



function plots_diagnostic( res; vn="K", i=-1, toplot="" ) 
  gr()
  pl = plot()
  if toplot == "" 
    if i != -1 
      toplot = Symbol(vn,"[$i]")
    else 
      toplot = Symbol(vn)
    end
  end
  pl = density!(pl, res[ toplot ])
  return pl
 end


# ----------

  

# -------------------


function fishery_model_predictions( res; prediction_time=prediction_time, solver_params=solver_params, PM=PM, 
  n_sample=-1, lower_bound=-0.001, override_negative_solution=true, ntries_mult=5 )

  nchains = size(res)[3]
  nsims = size(res)[1]
  
  if n_sample == -1
    # do all
    n_sample = nchains * nsims
  end

  md = zeros(nM, nS, n_sample, 2)  # number normalized
  mn = zeros(nM, nS, n_sample, 2)  # numbers
  mb = mn[:,1,:,:]  # biomass of first class
 
  trace_time = collect( solver_params.tspan[1]:solver_params.dt:solver_params.tspan[2] )
  ntt = length(trace_time)
   
  trace_bio = zeros(3, ntt, n_sample ) 
  trace = zeros( nS, ntt, n_sample )

  ntries = 0
  z = 0

  while z <= n_sample 
    ntries += 1
    (ntries > ntries_mult * n_sample) && break 
    (z >= n_sample) && break

    j = rand(1:nsims)  # nsims
    l = rand(1:nchains) # nchains

    b5 = [ res[j, Symbol("b5[$k]"), l] for k in 1:PM.nB]
    b6 = [ res[j, Symbol("b6[$k]"), l] for k in 1:PM.nB]
    K = [ res[j, Symbol("K[$k]"), l] for k in 1:PM.nS]
    v = [ res[j, Symbol("v[$k]"), l] for k in 1:4]
    d = [ res[j, Symbol("d[$k]"), l] for k in 1:PM.nS]
    d2= [ res[j, Symbol("d2[$k]"), l] for k in 1:PM.nS]
    u0 = [ res[j, Symbol("u0[$k]"), l] for k in 1:PM.nS]

    prb = remake( solver_params.prob; u0=u0, h=solver_params.h, tspan=solver_params.tspan, p=( b5, b6, K, d, d2, v ) )
    msol1 = solve( prb, solver_params.solver, saveat=solver_params.dt, callback=solver_params.cb ) #, # with callback == with fishing
    #  isoutofdomain=(y,p,t)->any(x -> x<lower_bound, y) )  # permit exceeding K, only check lower

    msol0 = solve( prb, solver_params.solver, saveat=solver_params.dt ) # ,  # no callback == no fishing
    #  isoutofdomain=(y,p,t)->any(x -> x<lower_bound, y) )  # permit exceeding K, only check lower

    if SciMLBase.successful_retcode(msol1) && SciMLBase.successful_retcode(msol0)
      
      ii0 = firstindexin( prediction_time, msol0.t )
      ii1 = firstindexin( prediction_time, msol1.t )
      
      jj0 = firstindexin( trace_time, msol0.t)
      jj1 = firstindexin( trace_time, msol1.t)
     
      if length(ii0) != PM.nM | length(ii1) != PM.nM
        continue
      end
      
      if length(jj0) != ntt | length(jj1) != ntt
        continue
      end

      # likelihood of the data
      MS0 = Array(msol0) 
      MS1 = Array(msol1) 
      
      iv0 = findall(x -> x<0.0, skipmissing(MS0) )
      iv1 = findall(x -> x<0.0, skipmissing(MS1) )

      if override_negative_solution
        # arguably, negative values suggest extirpation and is something to track ..
        MS0[ iv0 ] .= 0.0
        MS1[ iv1 ] .= 0.0
      else 
        if length(iv0) > 0 | length(iv1) > 0 
          continue
        end
      end

      z += 1

      sf  = nameof(typeof(mw)) == :ScaledInterpolation ? mw(msol1.t[ii1])  ./ 1000.0 ./ 1000.0  :  scale_factor   # n to kt
      
      md[:,:,z,1] = MS1[:,ii1]'  # with fishing
      md[:,:,z,2] = MS0[:,ii0]'  # witout fishing
      mn[:,:,z,1] = MS1[:,ii1]'  .* K' # with fishing scaled to K
      mn[:,:,z,2] = MS0[:,ii0]'  .* K' # without fishing scaled to K
      mb[:,z,1] = mn[:,1,z,1]  .* sf  # biomass of state var 1 
      mb[:,z,2] = mn[:,1,z,2]  .* sf

      # traces for plotting, etc 
      sft  = nameof(typeof(mw)) == :ScaledInterpolation ? mw(trace_time)  ./ 1000.0 ./ 1000.0 :  scale_factor
      trace_bio[1,:,z] = MS1[1,jj1] .* K[1] .* sft  # with 
      trace_bio[2,:,z] = MS0[1,jj0] .* K[1] .* sft  # without
      trace_bio[3,:,z] = ( trace_bio[2,:,z] .- trace_bio[1,:,z] ) ./ trace_bio[2,:,z] 
 
      trace[:,:,z] = MS1[:,jj1] .* K
    end # if
  end  # while
 
  if z < n_sample 
    @warn  "Insufficient number of solutions" 
  end
 
  sols= findall( x -> x!=0, vec( sum(mb, dims=(1,3,4) )) )
  if length(sols) > 1 
    md = md[:,:,sols,:]
    mn = mn[:,:,sols,:]
    mb = mb[:,sols,:,:]
  end
 
  return (md, mn, mb, trace, trace_bio, trace_time )

end



# -----------


function fishery_model_mortality( ; removed=removed, bio=bio, survey_time=survey_time, fish_year=fish_year )   
  fb = bio[1:length(survey_time),:,1]  # the last 1 is for size struct; no effect in discrete 
  removed_annual_kt = removals_aggregate( removed, fish_year )
  Fkt = removed_annual_kt[:,:rem_sum] ./1000.0 ./ 1000.0  # removal in kg -> kt
  FR =  Fkt ./ ( Fkt .+  fb )  # relative F
  FM = -1 .* log.(  1.0 .- min.( FR, 0.99) )  # instantaneous F
  FM[ FM .< eps(0.0)] .= zero(eltype(FM))
  return ( Fkt, FR, FM  )
end


# -----------




function abundance_from_index( S, res; k=1, model_variation="size_structured_dde_normalized"  )
  # map S <=> m  where S = observation index on unit scale; m = latent, scaled abundance on unit scale
  # observation model: S = (m - q0)/ q1   <=>   m = S * q1 + q0  

  if model_variation=="size_structured_dde_normalized"  
    q0 =  vec(res[:,Symbol("q0[$k]"),:])'
    q1 =  vec(res[:,Symbol("q1[$k]"),:])'
    K =  vec(res[:,Symbol("K[$k]"),:])'
    S_m = (( S[:,k]  .* q1) .+ q0 ) .*  K   # abundance_from_index  now on latent scale
  elseif model_variation=="size_structured_dde_unnormalized"  
    q0 =  vec(res[:,Symbol("q0[$k]"),:])'
    q1 =  vec(res[:,Symbol("q1[$k]"),:])' 
    S_m = (( S[:,k]  .* q1) .+ q0 )    # abundance_from_index  now on latent scale
  end
  return S_m
end


# -----------


function fishery_model_plot(; toplot=("fishing", "nofishing", "survey"), n_sample=min(250, size(bio)[2]),
  res=res, bio=bio, num=num, trace=trace, trace_bio=trace_bio, FM=FM, 
  S=S, si=1, scale_factor=scale_factor, 
  prediction_time=prediction_time, survey_time=survey_time, yrs=yrs, 
  alphav=0.05, labelsize=16,
  pl= Plots.plot(), 
  time_range=(floor(minimum(survey_time))-0.25, ceil(maximum(survey_time))+0.25 ),
  time_range_predictions=(floor(minimum(survey_time))-1.0, ceil(maximum(prediction_time)) )
)

  trace_time = collect( solver_params.tspan[1]:solver_params.dt:solver_params.tspan[2] )
   
  nsims = size(bio)[2]
  ss = rand(1:nsims, n_sample)  # sample index

  if any(isequal.("nofishing", toplot))  
    g = bio[:,:,2]   # [ yr,  sim, (with fishing=1; nofishing=2) ]
    pl = plot!(pl, prediction_time, g[:,ss] ;  alpha=alphav, color=:lime)
    pl = plot!(pl, prediction_time, mean(g, dims=2);  alpha=0.8, color=:limegreen, lw=4)
    pl = plot!(pl; legend=false )
    pl = plot!(pl; ylim=(0, maximum(g)*1.01 ) )
    pl = plot!(pl; xlim=time_range )
  end

  if any(isequal.("fishing", toplot))  
    g = bio[:,:,1]   # [ yr,  sim, (with fishing=1; nofishing=2) ]
    pl = plot!(pl, prediction_time, g[:,ss] ;  alpha=alphav, color=:orange)
    pl = plot!(pl, prediction_time, mean(g, dims=2);  alpha=0.8, color=:darkorange, lw=4)
    pl = plot!(pl; legend=false )
    pl = plot!(pl; ylim=(0, maximum(g)*1.01 ) )
    pl = plot!(pl; xlim=time_range )
  end

  if any(isequal.("footprint", toplot))  
    g1 = bio[:,:,1]   # [ yr,  sim, (with fishing=1; nofishing=2) ]
    g2 = bio[:,:,2]   # [ yr,  sim, (with fishing=1; nofishing=2) ]
    g = ( g2 - g1 ) ./ g2
    pl = plot!(pl, prediction_time, g[:,ss] ;  alpha=alphav, color=:lightslateblue)
    pl = plot!(pl, prediction_time, mean(g, dims=2);  alpha=0.8, color=:darkslateblue, lw=4)
    pl = plot!(pl; legend=false )
    pl = plot!(pl; ylim=(0, maximum(g)*1.01 ) )
    pl = plot!(pl; xlim=time_range )
    
  end

  if any(isequal.("survey", toplot))  
 
    # map S <=> m  where S = observation index on unit scale; m = latent, scaled abundance on unit scale
    # observation model: S = (m - q0)/ q1   <=>   m = S * q1 + q0  
    S_K_sims = abundance_from_index( S, res; k=1 )
 
    if nameof(typeof(mw)) == :ScaledInterpolation
      S_K_sims = S_K_sims .* mw(yrs) ./ 1000.0  ./ 1000.0
    else
      S_K_sims = S_K_sims .* scale_factor
    end

    S_K = mean(S_K_sims, dims=2)  # average by year
    pl = plot!(pl, survey_time, S_K, color=:gray, lw=2 )
    pl = scatter!(pl, survey_time, S_K, markersize=4, color=:darkgray)
    pl = plot!(pl; legend=false )
    pl = plot!(pl; xlim=time_range )
  end
  

  if any(isequal.("number", toplot))  

    gk = num[:,si,:,1]
    pl = plot!(pl, prediction_time, gk[:,ss];  alpha=alphav, color=:lightslateblue)
    pl = plot!(pl, prediction_time, mean(gk, dims=2);  alpha=0.8, color=:darkslateblue, lw=4)
    pl = plot!(pl; legend=false )
    pl = plot!(pl; ylim=(0, maximum(gk)*1.01 ) )

    S_K_sims = abundance_from_index( S, res; k=si )
    S_K = mean(S_K_sims, dims=2)  # average by year
 
    pl = plot!(pl, survey_time, S_K, color=:gray, lw=2 )
    pl = scatter!(pl, survey_time, S_K, markersize=4, color=:grey)
    pl = plot!(pl; legend=false )
    pl = plot!(pl; xlim=time_range )
  end

  if any(isequal.("number_trace", toplot))  

    gk = trace_bio[si,:,:]
    pl = plot!(pl, trace_time, gk[:,ss];  alpha=alphav, color=:lightslateblue)
    pl = plot!(pl, trace_time, mean(gk, dims=2);  alpha=0.8, color=:darkslateblue, lw=4)
    pl = plot!(pl; legend=false )
    pl = plot!(pl; ylim=(0, maximum(gk)*1.01 ) )

    S_K_sims = abundance_from_index( S, res; k=si )
    S_K = mean(S_K_sims, dims=2)  # average by year

    pl = plot!(pl, survey_time, S_K, color=:gray, lw=2 )
    pl = scatter!(pl, survey_time, S_K, markersize=4, color=:grey)
    pl = plot!(pl; legend=false )
    pl = plot!(pl; xlim=time_range )
  end


  if any(isequal.("trace", toplot))  
    pl = plot!( pl, trace_time, trace_bio[2,:,ss], alpha=alphav, lw=1, color=:lime )
    pl = plot!( pl, trace_time, trace_bio[1,:,ss], alpha=alphav, lw=1, color=:orange )
    pl = plot!( pl; legend=false )
    pl = plot!( pl; xlim=time_range )
  end

  if any(isequal.("trace_projections", toplot))  
    pl = plot!( pl, trace_time, trace_bio[2,:,ss], alpha=alphav, lw=1, color=:lime )
    pl = plot!( pl, trace_time, trace_bio[1,:,ss], alpha=alphav, lw=1, color=:orange )
    pl = plot!( pl, trace_time, mean(trace_bio[2,:,ss], dims=2);  alpha=0.8, color=:limegreen, lw=4)
    pl = plot!( pl, trace_time, mean(trace_bio[1,:,ss], dims=2);  alpha=0.8, color=:darkorange, lw=4)
    pl = plot!( pl; legend=false )
    pl = plot!( pl; xlim=time_range_predictions )
    pl = vline!( pl, year_assessment .+ [1,2],  alpha=0.5, color=:lightslategray, line=:dash, lw=1.5 )
  end

  if any(isequal.("trace_nofishing", toplot))  
    pl = plot!( pl, trace_time, trace_bio[2,:,ss], alpha=alphav, lw=1, color=:lime )
    pl = plot!(pl; legend=false )
    pl = plot!(pl; xlim=time_range )
  end

  if any(isequal.("trace_fishing", toplot))  
    pl = plot!( pl, trace_time, trace_bio[1,:,ss], alpha=alphav, lw=1, color=:orange )
    pl = plot!(pl; legend=false )
    pl = plot!(pl; xlim=time_range )
  end
 
  if any(isequal.("trace_footprint", toplot))  
    pl = plot!( pl, trace_time, trace_bio[3,:,ss], alpha=alphav, lw=1, color=:lightslateblue )
    pl = plot!(pl; legend=false )
    pl = plot!(pl; xlim=time_range )
  end

  if any(isequal.("trace_footprint_projections", toplot))  
    pl = plot!( pl, trace_time, trace_bio[3,:,ss], alpha=alphav, lw=1, color=:lightslateblue )
    pl = plot!(pl; legend=false )
    pl = plot!(pl; xlim=(floor(minimum(survey_time))-1.0, ceil(maximum(prediction_time)) ) )
    pl = vline!( pl, year_assessment .+ [1,2],  alpha=0.5, color=:lightslategray, line=:dash, lw=1.5 )
    pl = plot!(pl; xlim=time_range_predictions )
  end

  if any(isequal.("fishing_mortality", toplot))  
    FMmean = mean( FM, dims=2)
    FMmean[isnan.(FMmean)] .= zero(eltype(FM))
    ub = maximum(FMmean) * 1.1
    pl = plot!(pl, survey_time, FM[:,ss] ;  alpha=0.02, color=:lightslateblue)
    pl = plot!(pl, survey_time, FMmean ;  alpha=0.8, color=:slateblue, lw=4)
    pl = plot!(pl, ylim=(0, ub ) )
    pl = plot!(pl ; legend=false )
    pl = plot!(pl; xlim=time_range )
  end


  if any(isequal.("fishing_mortality_vs_footprint", toplot))  
    FMmean = mean( FM, dims=2)
    FMmean[isnan.(FMmean)] .= zero(eltype(FM))
    ub = maximum(FMmean) * 1.1
    g1 = bio[:,:,1]   # [ yr,  sim, (with fishing=1; nofishing=2) ]
    g2 = bio[:,:,2]   # [ yr,  sim, (with fishing=1; nofishing=2) ]
    g = ( g2 - g1 ) ./ g2
    g = g[1:length(survey_time),:]
    pl = scatter!(pl, FM[:,ss], g[:,ss];  alpha=alphav, color=:lightslateblue)
    pl = scatter!(pl, FMmean, mean(g, dims=2);  
      alpha=0.8, color=:darkslateblue, lw=4, markersize=4, markerstrokewidth=0,
      series_annotations = text.(trunc.(Int, survey_time), :top, :left, pointsize=4))
    pl = plot!(pl ; legend=false )
  end


  if any(isequal.("harvest_control_rule_footprint", toplot))  
    fb = bio[1:length(survey_time),:,1] 
 
    # mean weight by year
    sf = nameof(typeof(mw)) == :ScaledInterpolation ?  mw(yrs) ./ 1000.0  ./ 1000.0 : scale_factor
  
    # sample and plot posterior K
    K = vec( Array(res[:, Symbol("K[1]"), :]) ) .* mean(sf)  # convert to biomass 
  
    pl = vline!(pl, K[ss];  alpha=0.05, color=:limegreen )
    pl = vline!(pl, K[ss]./2;  alpha=0.05, color=:darkkhaki )
    pl = vline!(pl, K[ss]./4;  alpha=0.05, color=:darkred )
  
    pl = vline!(pl, [mean(K)];  alpha=0.6, color=:chartreuse4, lw=5 )
    pl = vline!(pl, [quantile(K, 0.975)];  alpha=0.5, color=:chartreuse4, lw=2, line=:dash )
    pl = vline!(pl, [quantile(K, 0.025)];  alpha=0.5, color=:chartreuse4, lw=2, line=:dash )
  
    pl = vline!(pl, [mean(K)/2.0];  alpha=0.6, color=:darkkhaki, lw=5 )
    pl = vline!(pl, [quantile(K, 0.975)]/2.0;  alpha=0.5, color=:darkkhaki, lw=2, line=:dash )
    pl = vline!(pl, [quantile(K, 0.025)]/2.0;  alpha=0.5, color=:darkkhaki, lw=2, line=:dash )
  
    pl = vline!(pl, [mean(K)/4.0];  alpha=0.6, color=:darkred, lw=5 )
    pl = vline!(pl, [quantile(K, 0.975)]/4.0;  alpha=0.5, color=:darkred, lw=2, line=:dash )
    pl = vline!(pl, [quantile(K, 0.025)]/4.0;  alpha=0.5, color=:darkred, lw=2, line=:dash )
  
    nt = length(survey_time)
    colours = get(ColorSchemes.tab20c, 1:nt, :extrema )[rand(1:nt, nt)]
  
    # scatter!( fb, FM ;  alpha=0.3, color=colours, markersize=4, markerstrokewidth=0)
  
    fb_mean = mean(fb, dims=2)

    g1 = bio[:,:,1]   # [ yr,  sim, (with fishing=1; nofishing=2) ]
    g2 = bio[:,:,2]   # [ yr,  sim, (with fishing=1; nofishing=2) ]
    g = ( g2 - g1 ) ./ g2
  
    g = g[1:length(survey_time),:]
  
    for ir in 1:size(g)[1]
      gm = mean(Iterators.filter(!isnan, g[ir,:]))
      g[ir, isnan.(g[ir,:])] .= gm
    end 

    g_mean = mean( g, dims=2)
  
    fbbb = [quantile(fb[nt,:], 0.025), quantile(fb[nt,:], 0.975) ]

    gbb = [quantile(g[nt,:], 0.975), quantile(g[nt,:], 0.025) ]
    
    pl = scatter!(pl, [fb[nt,:]], [g[nt,:]] ;  alpha=0.01, color=:lightslayeblue, markersize=2.5, markerstrokewidth=0)
    pl = scatter!(pl, fbbb, gbb;  alpha=0.5, color=:lightslategrey, markershape=:star, markersize=6, markerstrokewidth=1)
     
    pl = scatter!(pl,  [fb_mean[nt]], [g_mean[nt]] ;  alpha=0.9, color=:gold, markersize=8, markerstrokewidth=1)
  
    pl = plot!(pl, fb_mean, g_mean ;  alpha=0.8, color=:slateblue, lw=3)
  
    pl = scatter!(pl,  fb_mean, g_mean ;  alpha=0.8, color=colours,  markersize=4, markerstrokewidth=0,
      series_annotations = text.(trunc.(Int, survey_time), :top, :left, pointsize=6) )
    pl = scatter!(pl,  [fb_mean[nt]], [g_mean[nt]] ;  alpha=0.8, color=:yellow, markersize=10, markerstrokewidth=2)
    
    ub = max( quantile(K, 0.75), maximum( fb_mean ) ) * 1.05
    pl = plot!(pl; legend=false, xlim=(0, ub ), ylim=(0, maximum(g_mean ) * 1.05  ) )
 
  end
   

  if any(isequal.("harvest_control_rule", toplot))  
    fb = bio[1:length(survey_time),:,1] 
 
    # mean weight by year
    sf = nameof(typeof(mw)) == :ScaledInterpolation ?  mw(yrs) ./ 1000.0  ./ 1000.0 : scale_factor
  
    # sample and plot posterior K
    K = vec( Array(res[:, Symbol("K[1]"), :]) ) .* mean(sf)  # convert to biomass 
  
    pl = vline!(pl, K[ss];  alpha=0.05, color=:limegreen )
    pl = vline!(pl, K[ss]./2;  alpha=0.05, color=:darkkhaki )
    pl = vline!(pl, K[ss]./4;  alpha=0.05, color=:darkred )
  
    pl = vline!(pl, [mean(K)];  alpha=0.6, color=:chartreuse4, lw=5 )
    pl = vline!(pl, [quantile(K, 0.975)];  alpha=0.5, color=:chartreuse4, lw=2, line=:dash )
    pl = vline!(pl, [quantile(K, 0.025)];  alpha=0.5, color=:chartreuse4, lw=2, line=:dash )
  
    pl = vline!(pl, [mean(K)/2.0];  alpha=0.6, color=:darkkhaki, lw=5 )
    pl = vline!(pl, [quantile(K, 0.975)]/2.0;  alpha=0.5, color=:darkkhaki, lw=2, line=:dash )
    pl = vline!(pl, [quantile(K, 0.025)]/2.0;  alpha=0.5, color=:darkkhaki, lw=2, line=:dash )
  
    pl = vline!(pl, [mean(K)/4.0];  alpha=0.6, color=:darkred, lw=5 )
    pl = vline!(pl, [quantile(K, 0.975)]/4.0;  alpha=0.5, color=:darkred, lw=2, line=:dash )
    pl = vline!(pl, [quantile(K, 0.025)]/4.0;  alpha=0.5, color=:darkred, lw=2, line=:dash )
  
    nt = length(survey_time)
    colours = get(ColorSchemes.tab20c, 1:nt, :extrema )[rand(1:nt, nt)]
  
    fb_mean = mean(fb, dims=2)
    fm_mean = mean(FM, dims=2)
  
    fbbb = [quantile(fb[nt,:], 0.025), quantile(fb[nt,:], 0.975) ]

    FMbb = [quantile(FM[nt,:], 0.975), quantile(FM[nt,:], 0.025) ]
    
    pl = scatter!(pl, [fb[nt,:]], [FM[nt,:]] ;  alpha=0.01, color=:magenta, markersize=2.5, markerstrokewidth=0)
    pl = scatter!(pl, fbbb, FMbb;  alpha=0.5, color=:magenta, markershape=:star, markersize=6, markerstrokewidth=1)

    pl = scatter!(pl,  [fb_mean[nt]], [fm_mean[nt]] ;  alpha=0.9, color=:gold, markersize=8, markerstrokewidth=1)
    
    pl = plot!(pl, fb_mean, fm_mean ;  alpha=0.8, color=:slateblue, lw=3)
    pl = scatter!(pl,  fb_mean, fm_mean;  alpha=0.8, color=colours,  markersize=4, markerstrokewidth=0  )
    pl = scatter!(pl,  fb_mean .+0.051, fm_mean .-0.0025;  alpha=0.8, color=colours,  markersize=0, markerstrokewidth=0,
      series_annotations = text.(trunc.(Int, survey_time), :top, :left, pointsize=8) )

    ub = max( quantile(K, 0.75), maximum( fb_mean )  ) * 1.05
    pl = plot!(pl; legend=false, xlim=(0, ub ), ylim=(0, maximum(fm_mean ) * 1.05  ) )
 
  end
  pl = plot!(pl, xguidefontsize=labelsize, yguidefontsize=labelsize, xtickfontsize=labelsize,ytickfontsize=labelsize )
   
  return(pl)

end
 

function fishing_pattern_from_data(  fish_time, removed, ny=5 ) 
  # choose last n years of fishing and model attack rate 
  # scale to 1 then rescale to % of catch by time

  f = DataFrame( fish_time=fish_time, removals=removed )  # number
  f.yr = floor.( f.fish_time)
  f = f[(f[!,:yr] .> maximum(f[!,:yr]) - ny ),:]
  
  annual = DataFrame(
    fsum = [sum(x[!,:removals]) for x in groupby(f, :yr)],
    yr =   [ x[1,:yr] for x in groupby(f, :yr)]
  )

  f = innerjoin( f, annual, on=:yr)
  f.sy = f.fish_time - f.yr
  dtt = 12
  f.sy = round.( floor.( f.sy * dtt ) / dtt, digits=3) 
  f.rem = f.removals ./ f.fsum
 
  allowmissing!(f)
  f[findall(x->x>0.99, f.rem),:rem] .= missing

  ff = DataFrame(
    sy = [ x[1,:sy] for x in groupby(f, :sy)],
    fa = [mean(x[!,:rem ]) for x in groupby(f, :sy)]
  )
  sort!(ff, [order(:sy)])
  replace!(ff.fa, missing=> 0)
    
  ff.fa = ff.fa / sum(ff.fa)

  gg = DataFrame( mon=0:1:12 )
  gg.sy = round.( floor.( gg.mon / 12 * dtt ) / dtt, digits=3) 
  
  fs = outerjoin( gg, ff, on=:sy)
  fs[findall(x->ismissing(x), fs.fa),:fa] .= 0
  sort!(fs, [order(:sy)])
  
  ifs = extrapolate( interpolate( fs.fa, (BSpline(Linear())) ), Interpolations.Flat() )
  fishing_pattern_seasonal_interpolation_function = Interpolations.scale(ifs, 0:1/12:1 )

  return fishing_pattern_seasonal_interpolation_function
end

function impute_mean!(v)
  m = mean(Iterators.filter(!isnan, v))
  v[isnan.(v)] .= m
  v
end

function project_with_constant_catch( res; solver_params=solver_params, PM=PM, Catch=0, ny_fishing_pattern=5 )
  
  # forward project assuming constant fishing pattern
  # callbacks for external perturbations to the system (deterministic fishing without error)
   
  sf = nameof(typeof(mw)) == :ScaledInterpolation ?  mw(yrs) ./ 1000.0  ./ 1000.0 : scale_factor
  
  # sample and plot posterior K
  K = vec( Array(res[:, Symbol("K[1]"), :]) ) .* mean(sf)  # convert to biomass 
  ER =  Catch / mean(K) 

  exploitationrate = exp(ER)-1.0  # relative to K

  fishing_pattern_seasonal = fishing_pattern_from_data(fish_time, removed, ny_fishing_pattern ) # fraction of annual total .. fishing_pattern(0.1) gives fraction captured on average by 0.1 * 365 days 
  
  function condition_fp(u, t, integrator )
    t in fish_time_project
  end

  function affect_fishing_project!(integrator)
    k = integrator.t - floor(integrator.t)
    integrator.u[1] -=  exploitationrate / solver_params.hsa(integrator.t,1) * fishing_pattern_seasonal(k)  # p[2] ==K divide by K[1]  .. keep unscaled to estimate magnitude of other components
  end
 
  sp = deepcopy(solver_params)
  sp = @set sp.cb =  CallbackSet(
    PresetTimeCallback( fish_time, affect_fishing! ),
    DiscreteCallback( condition_fp, affect_fishing_project!, save_positions=(true, true) )
  );
  
  # n scaled, n unscaled, biomass of fb with and without fishing, model_traces, model_times 
  # override_negative_solution=true makes it more permissive .. but tuncates at 0, so be careful
  m, num, bio, trace, trace_bio, trace_time = fishery_model_predictions(res, solver_params=sp, PM=PM, ntries_mult=10 )

  return m, num, bio, trace, trace_bio, trace_time
end

 


function plot_prior_posterior( vn, prior, posterior; bw=0.02 )
  pri =  vec(collect( prior[:,Symbol(vn),:] ))
  pos =  vec(collect( posterior[:,Symbol(vn),:] ))
  pl = plot(ylabel="Density", xlabel=vn ) 
  pl = density!(pl, pri,  fill=true, color = :slateblue, fillalpha=0.25, bandwidth = bw, lw=0, label="Prior")
  pl = density!(pl, pos,  fill=true, color = :purple, fillalpha=0.5, bandwidth = bw, lw=0, label="Posterior")
  return(pl)
end

 
 

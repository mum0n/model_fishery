



function logistic_discrete_data( yrs, aulab, fishery_data, model_variation )
 
  Y = fishery_data["Y"][∈(yrs).(fishery_data["Y"].yrs), :]
  removals = fishery_data["L"][∈(yrs).(fishery_data["L"].yrs), :]
  Kmu = [5.0, 65.0, 1.5]
  nT = length(yrs)
  nP = 5  # number of predictions into future (with no fishing)
  nM = nP + nT  # total number of prediction years

  dt = 1  # time resolution of solutions .. discrete annual so 1 year
  nS = 1 # no state variables, not used 

  no_digits = 3  # time floating point rounding 
  smallnumber = 1.0e-9 # floating point value sufficient to assume 0 valued
  
  # "survey index"
  S = Y[:,Symbol("$aulab"  )]

  # scale index where required
  Smean = mean(skipmissing(S))
  Sstd = std( skipmissing(S))
  Smin = minimum(skipmissing(S))
  Smax = maximum( skipmissing(S))
  Srange = Smax - Smin 

  SminFraction = Smin ./ Srange  # used as informative prior mean in some runs

  if model_variation=="logistic_discrete_historical"
  # do nothing (no scaling)
  elseif occursin( r"scale_center", model_variation ) 
  S = (S .- Smean ) ./ Sstd  # scale to std and center to 0 
  else 
  S = (S .- Smin ) ./ Srange    # default is to scale (min,max) -> (0,1)
  end


  # id index
  ki = aulab=="cfanorth" ? 1 :
      aulab=="cfasouth" ? 2 :
      aulab=="cfa4x"    ? 3 :
      0  # default

  kmu = Kmu[ki]  

  smallnumber = 1.0 / kmu / 10.0  # floating point value of sufficient to assume 0 valued
  no_digits = 3  # time floating point rounding

  survey_time = Y[:,:yrs]    # time of observations for survey .. for plotting
  # adjust survey time to be approximately in spring or autumn
  spring = findall( x -> x < 2004, survey_time )
  fall   = findall( x -> x >= 2004, survey_time )
  survey_time[spring] = survey_time[spring] .+ 5.0/12.0  # arbitrary ..  "spring"
  survey_time[fall]   = survey_time[fall] .+ 10.0/12.0   # time of survey in "fall" 
  survey_time =  round.(  survey_time ; digits=no_digits)    # time of observations for survey

  removed = removals[:,Symbol("$aulab")]

  predtime =  4.0/12.0  # predictions ("m") are "prefishery" .. arbitrarily setting to 4/12
  prediction_time = floor.( vcat( collect(minimum(yrs) : (maximum(yrs)+nP) ) )  ) .+  round( predtime/dt ; digits=no_digits)   # april (m== prefishery)
  # prediction_time = floor.(vcat( survey_time, collect(1:nP) .+ maximum(survey_time) ) )     
  yrs_pred_report = findall( x -> (x >= 1999.0) & (x <= (Real(maximum(yrs))+1.0)), prediction_time )
  prediction_time_ss = prediction_time[yrs_pred_report]

  iok = findall( !ismissing, S )

  # basic params for "logistic_discrete"
  PM = (
      yrs = yrs,
      nS = nS, 
      nT = length(yrs),
      nP = nP,  # number of predictions into future (with no fishing)
      nM = nM,  # total number of prediction years
      kmu = kmu,
      K = (kmu, 0.25*kmu, kmu/5.0, kmu*5.0 ),
      r = (1.0, 0.1, 0.5, 1.5),
      bpsd = ( 0.1, 0.05, 0.01, 0.5 ),
      bosd = ( 0.1, 0.05, 0.01, 0.5 ),
      q1 = (  1.0, 0.2,  0.01, 10.0 ),
      q0 = ( SminFraction, 0.1, -1.0, 1.0 ),
      m0 = ( 0.9, 0.2, 0.1, 1.0), 
      mlim =(0.0, 1.0),
      removed=removed,
      S = S,
      iok = iok,
      yeartransition = 0
  ) 

  # translate model-specific functions, etc to generics
  if model_variation=="logistic_discrete_historical"

    if (aulab=="cfanorth")  
        PM = @set PM.yeartransition = 6
        PM = @set PM.K = PM.kmu .* (1.0, 0.1, 1.0/5.0, 5.0 )
        PM = @set PM.r = (1.0, 0.1, 0.25, 3.0)
        PM = @set PM.bpsd = ( 0.1, 0.1, 0.01, 0.5 )
        PM = @set PM.bosd = PM.kmu .* ( 0.1, 0.1, 0.01, 0.75 )
        PM = @set PM.q1 = ( 1.0, 0.1, 0.01, 2.0 )
        PM = @set PM.mlim =( 0.01, 1.25 )
        PM = @set PM.m0 = (0.5, 0.25, 0.0, 1.25 )
    end

    if (aulab=="cfasouth")
        PM = @set PM.yeartransition = 6
        PM = @set PM.K = PM.kmu .* (1.0, 0.1, 1.0/5.0, 5.0 )
        PM = @set PM.r = (1.0, 0.1, 0.25, 3.0)
        PM = @set PM.bpsd = ( 0.1, 0.1, 0.01, 0.5 )
        PM = @set PM.bosd = PM.kmu .* ( 0.1, 0.1, 0.01, 0.75 )
        PM = @set PM.q1 = ( 1.0, 0.1, 0.01, 2.0 )
        PM = @set PM.mlim =( 0.01, 1.25 )
        PM = @set PM.m0 = (0.5, 0.25, 0.0, 1.25 )
     end

    if (aulab=="cfa4x")  
        PM = @set PM.yeartransition = 0
        PM = @set PM.K = PM.kmu .* (1.0, 0.1, 1.0/5.0, 5.0 )
        PM = @set PM.r = (1.0, 0.1, 0.25, 3.0)
        PM = @set PM.bpsd = ( 0.1, 0.1, 0.01, 0.5 )
        PM = @set PM.bosd = PM.kmu .* ( 0.1, 0.1, 0.01, 0.75 )
        PM = @set PM.q1 = ( 1.0, 0.1, 0.01, 2.0 )
        PM = @set PM.mlim =( 0.01, 1.25 )
        PM = @set PM.m0 = (0.5, 0.25, 0.0, 1.25 )
     end
 
  elseif model_variation=="logistic_discrete_basic"
  # same as historical but normalize to reduce influce of large magnitudes and faster convergence

      if (aulab=="cfanorth")  
          PM = @set PM.yeartransition = 6
          PM = @set PM.K = (PM.kmu, 0.25*PM.kmu, PM.kmu/5.0, PM.kmu*5.0 )
          PM = @set PM.r = (1.0, 0.1, 0.25, 3.0)
          PM = @set PM.bpsd = ( 0.1, 0.1, 0.01, 0.5 )
          PM = @set PM.bosd = ( 0.1, 0.1, 0.01, 0.75 )
          PM = @set PM.q1 = (  1.0, 0.1,  0.25, 1.25 )
          PM = @set PM.mlim =( 0.1, 1.25 )
          PM = @set PM.m0 = ( 0.9, 0.2, 0.2, 1.25)
      end

      if (aulab=="cfasouth")
          PM = @set PM.yeartransition = 6
          PM = @set PM.K = (PM.kmu, 0.25*PM.kmu, PM.kmu/5.0, PM.kmu*5.0 )
          PM = @set PM.r = (1.0, 0.1, 0.25, 3.0)
          PM = @set PM.bpsd = ( 0.1, 0.1, 0.001, 0.5 )
          PM = @set PM.bosd = ( 0.1, 0.1, 0.001, 0.75 )
          PM = @set PM.q1 = (  1.0, 0.1,  0.25, 1.25 )
          PM = @set PM.mlim =( 0.2, 1.25 )
          PM = @set PM.m0 = ( 0.5, 0.2, 0.0, 1.25)
      end

      if (aulab=="cfa4x")  

          PM = @set PM.yeartransition = 0
          PM = @set PM.K = (PM.kmu, 0.2*PM.kmu, PM.kmu/5.0, PM.kmu*5.0 )
          PM = @set PM.r = (1.0, 0.2, 0.25, 3.0)
          PM = @set PM.bpsd = ( 0.1, 0.1, 0.001, 0.5 )
          PM = @set PM.bosd = ( 0.1, 0.1, 0.001, 0.75 )
          PM = @set PM.q1 = ( 1.0, 0.2, 0.25, 1.25 )
          PM = @set PM.mlim =( 0.1, 1.25 )
          PM = @set PM.m0 = ( 0.5, 0.2, 0.25, 1.25)
      end
 
  elseif model_variation=="logistic_discrete"

      if (aulab=="cfanorth") | (aulab=="cfasouth")
          PM = @set PM.yeartransition = 6
      end
 
  elseif model_variation=="logistic_discrete_map"

      # not used .. just for testing
      if (aulab=="cfanorth") | (aulab=="cfasouth")
          PM = @set PM.yeartransition = 6
      end
      PM = @set PM.K = (PM.kmu, 0.2*PM.kmu, PM.kmu/5.0, PM.kmu*5.0 )
      PM = @set PM.r = (1.0, 0.1, 0.5, 3.0)
      PM = @set PM.bpsd = (0, 0.05, 0.01, 0.5) 
      PM = @set PM.bosd = (0, 0.05, 0.01, 0.75)
      PM = @set PM.q1 = (1.0, 0.1,  0.01, 10.0)
      PM = @set PM.q0 = (0.0, 0.1, -1.0, 1.0)
   end
 
  # translate model-specific functions, etc to generics
  if model_variation=="logistic_discrete_historical"

    fmod = logistic_discrete_turing_historical(PM)  # q1 only

  elseif model_variation=="logistic_discrete_basic"
  # same as historical but normalize to reduce influce of large magnitudes and faster convergence

    fmod = logistic_discrete_turing_basic(PM)

  elseif model_variation=="logistic_discrete"

    fmod =  logistic_discrete_turing(PM)

  elseif model_variation=="logistic_discrete_map"

      # not used .. just for testing
      fmod = logistic_discrete_map_turing(PM)
  end
 
  return PM, fmod
end

 


Turing.@model function logistic_discrete_turing( PM )
  # biomass process model: dn/dt = r n (1-n/K) - removed ; b, removed are not normalized by K  

  K ~ truncated(Normal( PM.K[1], PM.K[2]), PM.K[3], PM.K[4])  
  r ~  truncated(Normal( PM.r[1], PM.r[2]), PM.r[3], PM.r[4])   # (mu, sd)
  bpsd ~  truncated(Normal( PM.bpsd[1], PM.bpsd[2]), PM.bpsd[3], PM.bpsd[4] )  ;  # slightly informative .. center of mass between (0,1)
  bosd ~  truncated(Normal( PM.bosd[1], PM.bosd[2]), PM.bosd[3], PM.bosd[4] )  ;  # slightly informative .. center of mass between (0,1)
  q1 ~ truncated(Normal( PM.q1[1], PM.q1[2]), PM.q1[3], PM.q1[4] )    
  q0 ~ truncated(Normal(PM.q0[1], PM.q0[2]), PM.q0[3], PM.q0[4]  ) 

  m = tzeros( PM.nM )
  m[1] ~  truncated(Normal( PM.m0[1], PM.m0[2]), PM.m0[3], PM.m0[4] )  ; # starting b prior to first catch event

  for i in 2: PM.nT
    m[i] ~ truncated(Normal( m[i-1] + r * m[i-1] * ( 1.0 - m[i-1] ) -  PM.removed[i-1]/K, bpsd), PM.mlim[1], PM.mlim[2])  ;
  end

  for i in ( PM.nT+1): PM.nM
    m[i] ~ truncated(Normal( m[i-1] + r * m[i-1] * ( 1.0 - m[i-1] ), bpsd), PM.mlim[1], PM.mlim[2])  ;  # predict with no removals
  end

  if any( x -> x < 0.0 || x >1.0, m)
    Turing.@addlogprob! -Inf
    return nothing
  end
 
  # likelihood
    # map S <=> m  where S = observation index on unit scale; m = latent, scaled abundance on unit scale
    # observation model: S = (m - q0)/ q1   <=>   m = S * q1 + q0  
    # see function: abundance_from_index      

  # m = abundance (prefishery)
  # survey timing: spring before 2004 and fall afterwards 
  # fishery operates in winter for 4x and spring and summer for N and Sens
  #   4x:  (m[i-1] - removals[i-1]) = abundance (post fishery) upon which dynamics is applied, to give m[i], the abundance (prefishery), so S[i] ~ m[i] - rem[i]
  #   n and sens: post 2004 .. same as above, , so S[i] ~ m[i] - rem[i] (surveys post fishery)
  #               pre 2004 ... same as above but.. S[i] ~ m[i]  (no removals due to survey being prefishery)
  if PM.yeartransition == 0
    # 4X
    for i in PM.iok
      PM.S[i] ~ Normal( (m[i] -  PM.removed[i]/K - q0 ) / q1, bosd )  ;
    end 
  else 
    # NENS, SENS
    for i in PM.iok
      if i < PM.yeartransition
        PM.S[i] ~ Normal( ( m[i] - q0 ) / q1, bosd )  ;  # spring survey
      elseif i == PM.yeartransition
        PM.S[i] ~ Normal( ( m[i] - (PM.removed[i-1] + PM.removed[i]) / (2.0*K ) - q0 ) / q1  , bosd )  ;  # transition year  .. averaging should be done before .. less computation
      else
        PM.S[i] ~ Normal( ( m[i] - PM.removed[i]/K - q0 ) / q1, bosd )  ; # fall survey
      end
    end
  
  end

end
 

Turing.@model function logistic_discrete_turing_basic( PM )
  # biomass process model: dn/dt = r n (1-n/K) - removed ; b, removed are not normalized by K  
  # priors 

  K ~ truncated(Normal( PM.K[1], PM.K[2]), PM.K[3], PM.K[4])  
  r ~  truncated(Normal( PM.r[1], PM.r[2]), PM.r[3], PM.r[4])   # (mu, sd)
  bpsd ~  truncated(Normal( PM.bpsd[1], PM.bpsd[2]), PM.bpsd[3], PM.bpsd[4] )  ;  # slightly informative .. center of mass between (0,1)
  bosd ~  truncated(Normal( PM.bosd[1], PM.bosd[2]), PM.bosd[3], PM.bosd[4] )  ;  # slightly informative .. center of mass between (0,1)
  q1 ~ truncated(Normal( PM.q1[1], PM.q1[2]), PM.q1[3], PM.q1[4] )    

  m = tzeros( PM.nM )
  m[1] ~  truncated(Normal( PM.m0[1], PM.m0[2]), PM.m0[3], PM.m0[4]  )  ; # starting b prior to first catch event

  for i in 2:PM.nT
    m[i] ~ truncated(Normal( m[i-1] + r * m[i-1] * ( 1.0 - m[i-1] ) - PM.removed[i-1]/K, bpsd), PM.mlim[1], PM.mlim[2])  ;
  end
  
  for i in (PM.nT+1):PM.nM
    m[i] ~ truncated(Normal( m[i-1] + r * m[i-1] * ( 1.0 - m[i-1] ), bpsd), PM.mlim[1], PM.mlim[2])  ; # predict with no removals
  end
    
  if any( x -> x < 0.0 || x >1.0, m)
    Turing.@addlogprob! -Inf
    return nothing
  end
   

  # likelihood
    # map S <=> m  where S = observation index on unit scale; m = latent, scaled abundance on unit scale
    # observation model: S = (m - q0)/ q1   <=>   m = S * q1 + q0  
    # see function: abundance_from_index      
  if PM.yeartransition == 0
    # 4X
    for i in PM.iok
      PM.S[i] ~ Normal( (m[i] - PM.removed[i]/K) / q1, bosd )  ;
    end 
  else
    for i in PM.iok
      if i < PM.yeartransition
        PM.S[i] ~ Normal( ( m[i] ) / q1 , bosd )  ;  # spring survey
      elseif i == PM.yeartransition
        PM.S[i] ~ Normal( ( m[i] - (PM.removed[i-1] + PM.removed[i]) / (2.0*K ) ) / q1 , bosd )  ;  # transition year  .. averaging should be done before .. less computation
      else
        PM.S[i] ~ Normal( ( m[i] - PM.removed[i]/K ) / q1 , bosd )  ; # fall survey
      end
    end
  
  end

end
  
    


Turing.@model function logistic_discrete_turing_historical( PM )
  # biomass process model: dn/dt = r n (1-n/K) - removed ; b, removed are not normalized by K  
  # priors 
  K ~ truncated(Normal( PM.K[1], PM.K[2]), PM.K[3], PM.K[4])  
  r ~  truncated(Normal( PM.r[1], PM.r[2]), PM.r[3], PM.r[4])   # (mu, sd)
  bpsd ~  truncated(Normal( PM.bpsd[1], PM.bpsd[2]), PM.bpsd[3], PM.bpsd[4] )  ;  # slightly informative .. center of mass between (0,1)
  bosd ~  truncated(Normal( PM.bosd[1], PM.bosd[2]), PM.bosd[3], PM.bosd[4] )    ;  # slightly informative .. center of mass between (0,1)
  q1 ~ truncated(Normal( PM.q1[1], PM.q1[2]), PM.q1[3], PM.q1[4] )    

  # m's are "total avaialble for fishery" (latent truth)
  m = tzeros( PM.nM )
  m[1] ~ truncated(Normal(PM.m0[1], PM.m0[2]), PM.m0[3], PM.m0[4])   ; # starting b prior to first catch event

  for i in 2:PM.nT
    m[i] ~ truncated(Normal( m[i-1] + r * m[i-1] * ( 1.0 - m[i-1] ) - PM.removed[i-1]/K, bpsd), PM.mlim[1], PM.mlim[2])  ;
  end

  for i in (PM.nT+1):PM.nM
    m[i] ~ truncated(Normal( m[i-1] + r * m[i-1] * ( 1.0 - m[i-1] ), bpsd), PM.mlim[1], PM.mlim[2])  ; # predict with no removals (prefishery)
  end
 
  if any( x -> x < 0.0, m)  # permit overshoot
    Turing.@addlogprob! -Inf
    return nothing
  end

  # likelihood
    # map S <=> m  where S = observation index on unit scale; m = latent, scaled abundance on unit scale
    # observation model: S = (m - q0)/ q1   <=>   m = S * q1 + q0  
    # see function: abundance_from_index      

  if PM.yeartransition == 0
    # 4X
    for i in PM.iok
      PM.S[i] ~ Normal(  (K * m[i] - PM.removed[i]) / q1, bosd )  ; # fall survey
    end

  else

    # spring to fall survey: transition year = 2004
    # spring = 1:5
    # fall = 6:last
    
    for i in PM.iok
      if  i < PM.yeartransition
        PM.S[i] ~ Normal(  K * m[i] / q1, bosd )  ;  # spring survey
      elseif i == PM.yeartransition
        PM.S[i] ~ Normal(  ( K * m[i] - (PM.removed[i-1] + PM.removed[i]) / 2.0) / q1, bosd )  ;  # transition year  .. averaging should be done before .. less computation 
      else
        PM.S[i] ~ Normal(  ( K * m[i] - PM.removed[i] ) / q1 , bosd )  ; # fall survey
      end
    end
  end

end
  


Turing.@model function logistic_discrete_map_turing( PM )
  # biomass process model: n(t+1) = r n (1-n/K) - removed ; b, removed are not normalized by K  
  # priors 

  K ~ truncated(Normal( PM.K[1], PM.K[2]), PM.K[3], PM.K[4])  
  r ~  truncated(Normal( PM.r[1], PM.r[2]), PM.r[3], PM.r[4])   # (mu, sd)
  bpsd ~  truncated(Normal( PM.bpsd[1], PM.bpsd[2]), PM.bpsd[3], PM.bpsd[4] )  ;  # slightly informative .. center of mass between (0,1)
  bosd ~  truncated(Normal( PM.bosd[1], PM.bosd[2]), PM.bosd[3], PM.bosd[4] )  ;  # slightly informative .. center of mass between (0,1)
  q1 ~ truncated(Normal( PM.q1[1], PM.q1[2]), PM.q1[3], PM.q1[4] )    
  q0 ~ truncated(Normal(PM.q0[1], PM.q0[2]), PM.q0[3], PM.q0[4]  ) 

  m = tzeros( PM.nM )
  m[1] ~  truncated(Normal( PM.m0[1], PM.m0[2]), PM.m0[3], PM.m0[4] )  ; # starting b prior to first catch event

  for i in 2:PM.nT
    m[i] ~ truncated(Normal(  r * m[i-1] * ( 1.0 - m[i-1] ) - PM.removed[i-1]/K, bpsd), PM.mlim[1], PM.mlim[2])  ;
  end

  for i in (PM.nT+1):PM.nM
    m[i] ~ truncated(Normal(  r * m[i-1] * ( 1.0 - m[i-1] ), bpsd), PM.mlim[1], PM.mlim[2])  ;  # predict with no removals
  end

  if any( x -> x < 0.0, m)  # permit overshoot
    Turing.@addlogprob! -Inf
    return nothing
  end

  # likelihood
    # map S <=> m  where S = observation index on unit scale; m = latent, scaled abundance on unit scale
    # observation model: S = (m - q0)/ q1   <=>   m = S * q1 + q0  
    # see function: abundance_from_index      
  if PM.yeartransition == 0

    for i in PM.iok
      if i == 1
        PM.S[i] ~ Normal( (m[i] - q0 ) / q1, bosd )  ;
      else
        PM.S[i] ~ Normal( (m[i] - PM.removed[i-1]/K - q0 ) / q1, bosd )  ;
      end
    end

  else 
    # NENS, SENS
    for i in PM.iok
      if i < PM.yeartransition
        PM.S[i] ~ Normal( ( m[i] - q0 ) / q1, bosd )  ;  # spring survey
      elseif i == PM.yeartransition
        PM.S[i] ~ Normal( ( m[i] - (PM.removed[i-1] + PM.removed[i]) / (2.0*K ) - q0 ) / q1, bosd )  ;  # transition year 
      else
        PM.S[i] ~ Normal( ( m[i] - PM.removed[i]/K - q0 ) / q1, bosd )  ; # fall survey
      end
    end
  end
  
end

 

  
function fishery_model_test( test=("basic" ) )

  ## test model by sampling from random priors 
  gr()
  pl = plot()

  if any( occursin.( r"basic", test )  )

    res = sample( fmod, Prior(), 100, nwarmup = 100, nchains =1 )
 
    for l in 1:size(res)[3]
      for i in 1:length(res)  
          w = zeros(nM)
          for j in 1:nM
              w[j] = res[i, Symbol("K"),l] * res[i, Symbol("m[$j]"),l] 
          end
          pl = plot!(pl, prediction_time, w;  alpha=0.1, color=:orange)
      end
    end
    pl = plot!(pl; legend=false, title="basic prior check" )
 
  end
 
  return (res, pl) 

end


 
# -------------------

# function expand_grid(; iters...)
#     var_names = collect(keys(iters))
#     var_itr = [1:length(x) for x in iters.data]
#     var_ix = vcat([collect(x)' for x in Iterators.product(var_itr...)]...)
#     out = DataFrame()
#     for i = 1:length(var_names)
#         out[:,var_names[i]] = collect(iters[i])[var_ix[:,i]]
#     end
#     return out
# end
#  expand_grid(a=1:2, b=1.0:5.0, c=["one", "two", "three", "four"])

function expand_grid(; kws...)
  names, vals = keys(kws), values(kws)
  return DataFrame(NamedTuple{names}(t) for t in Iterators.product(vals...))
end
 

function discretize_decimal( x, delta=0.01 ) 
  num_digits = Int(ceil( log10(1.0 / delta)) )   # time floating point rounding
  out = round.( round.( x ./ delta; digits=0 ) .* delta; digits=num_digits)
  return out
end



function fishery_model_predictions( res; prediction_time=prediction_time, n_sample=-1 )

  nchains = size(res)[3]
  nsims = size(res)[1]
 
  if n_sample == -1
    # do all
    n_sample = nchains * nsims
    oo = expand_grid( sims=1:nsims, chains=1:nchains)
  else
    oo = DataFrame( sims=rand(1:nsims, n_sample), chains=rand(1:nchains, n_sample) )
  end

  md = zeros(nM, n_sample) 
  mb = zeros(nM, n_sample)
   
  z = 0

  while z <= n_sample 
    z += 1
    z > n_sample && break
    j = oo[z, :sims]  # nsims
    l = oo[z, :chains] # nchains
    for i in 1:nM
      md[i,z] = res[j, Symbol("m[$i]"), l]
      mb[i,z] = md[i,z] * res[j, Symbol("K"), l]
    end
  end

  # additional nothings to keep same expectations as continuous models
  return (md, nothing, mb, nothing, nothing, nothing )  

end



# -----------


function fishery_model_predictions_timeseries(nothing; prediction_time=prediction_time, plot_k=1)
  print("Nothing to do in model with a single state variable")
  return (nothing, nothing)
end


# ----------

function fishery_model_predictions_trace( res; n_sample=10, plot_k=1, alpha=0.01, plot_only_fishing=true )
  print("Nothing to do in model with a single state variable")
  return (nothing, nothing, nothing)
end


# ----------

function logistic_discrete_reference_points(r, K)
  expK = exp.(K) 
  msy   = r .* expK ./ 4.0 ; # maximum height of of the latent productivity (yield)
  bmsy  = expK ./ 2.0 ; # biomass at MSY
  fmsy  = 2.0 .* msy ./ expK ; # fishing mortality at MSY
  return (msy, bmsy, fmsy)
end

      # plot!(prediction_time, u, lw=2, color=:orangered )
      # scatter!(prediction_time, u, markersize=4, color=:goldenrod1 )

      # plot!(survey_time, yhat, color=:purple2, lw=2 )
      # scatter!(survey_time, yhat, markersize=4, color=:purple4)
  

function fishing_mortality_instantaneous( removed, abundance )
  -log(  1.0 - (removed  / abundance)  )  ;
end

 

function showall( x )
    # print everything to console
    show(stdout, "text/plain", x) # display all estimates
end


function plots_diagnostic( res, vn="K" ) 
  gr()
  pl = plot()
  pl = density!(pl, res[ Symbol(vn) ])
  return pl
 
end



# ----------
    

function fishery_model_mortality(; removed=removed, bio=bio, survey_time=survey_time )    
  fb = bio[1:length(survey_time),:,1]  # the last 1 is for size struct; no effect in discrete 
  Fkt = removed
  FR =  Fkt ./ ( Fkt .+  fb )  # relative F
  FM = -1 .* log.(  1.0 .- min.( FR, 0.99) )  # instantaneous F
  # FM[ FM .< eps(0.0)] .= zero(eltype(FM))
  return ( Fkt, FR, FM  )
end



# -----------


function abundance_from_index( Sai, res, model_variation="logistic_discrete_historical" )
  # map S <=> m  where S = observation index on unit scale; m = latent, scaled abundance on unit scale
  # observation model: S = (m - q0)/ q1   <=>   m = S * q1 + q0  

  K =  vec( res[:,Symbol("K"),:] )'
  
  # if nameof(typeof(mw)) == :ScaledInterpolation
  #   Sbacktransf = Sbacktransf .* mw(yrs) ./ 1000.0  ./ 1000.0

  if model_variation=="logistic_discrete_basic"  
    q1 =  vec( res[:,Symbol("q1"),:] )'
    S_m = Sai .* q1 .* K
  elseif model_variation=="logistic_discrete_historical"  
    q1 =  vec( res[:,Symbol("q1"),:] )'
    S_m = Sai .* q1   # already on K scale
  elseif model_variation=="logistic_discrete"  
    q0 =  vec( res[:,Symbol("q0"),:] )'
    q1 =  vec( res[:,Symbol("q1"),:] )'
    S_m = ( Sai .* q1 .+ q0 )  .* K
  elseif model_variation=="logistic_discrete_map"  
    q0 =  vec( res[:,Symbol("q0"),:] )'
    q1 =  vec( res[:,Symbol("q1"),:] )'
    S_m = (Sai .* q1 .+ q0 ) .* K
  end

  return S_m
end

# -----------



function fishery_model_plot(; toplot=("fishing", "survey"), n_sample=min(250, size(bio)[2]),
  res=res, bio=bio, FM=FM, 
  S=S,
  prediction_time=prediction_time, prediction_time_ss=prediction_time_ss, survey_time=survey_time, yrs=yrs, 
  alphav=0.075, pl= plot(), time_range=(floor(minimum(prediction_time_ss))-1.0, ceil(maximum(prediction_time_ss))+0.5 )
)
 
  nsims = size(bio)[2]
  ss = rand(1:nsims, n_sample)  # sample index

  if any(isequal.("trace", toplot))  
    @warn "trace is not valid for a discrete model"
    
  end 

  if any(isequal.("nofishing", toplot))  
    @warn "nofishing not implemented"
    
  end 

  # extract sims (with fishing)
  # plot biomass
  if any(isequal.("fishing", toplot))   # this means there is fishing occuring ( and plot biomass )
    g = bio   # [ yr,  sim ]
    pl = plot!(pl, prediction_time, g[:,ss] ;  alpha=alphav, color=:orange)
    pl = plot!(pl, prediction_time, mean(g, dims=2);  alpha=0.8, color=:darkorange, lw=4)
    pl = plot!(pl; legend=false )
    pl = plot!(pl; ylim=(0, maximum(g)*1.01 ) )
    pl = plot!(pl; xlim=time_range )
  end
 

  if any(isequal.("footprint", toplot))  
    @warn "footprint not implemented"
    
  end
   

  if any(isequal.("survey", toplot))  
    # map S -> m and then multiply by K
    # where S=observation on unit scale; m=latent, scaled abundance on unit scale
    S_m = abundance_from_index( S, res  )
    S_K = mean(S_m, dims=2)  # average by year
    pl = plot!(pl, survey_time, S_K, color=:gray, lw=2 )
    pl = scatter!(pl, survey_time, S_K, markersize=4, color=:darkgray)
    pl = plot!(pl; legend=false )
    pl = plot!(pl; xlim=time_range )

  end
   

  if any(isequal.("fishing_mortality", toplot))  
    FMmean = mean( FM, dims=2)
    FMmean[isnan.(FMmean)] .= zero(eltype(FM))
    ub = maximum(FMmean) * 1.1
    pl = plot!(pl, prediction_time_ss, FM[:,ss] ;  alpha=0.02, color=:lightslateblue)
    pl = plot!(pl, prediction_time_ss, FMmean ;  alpha=0.8, color=:slateblue, lw=4)
    pl = plot!(pl, ylim=(0, ub ) )
    pl = plot!(pl ; legend=false )
    pl = plot!(pl; xlim=time_range )
  end


  if any(isequal.("fishing_mortality_vs_footprint", toplot))  
    @warn "footprint not implemented"
    

  end


  if any(isequal.("harvest_control_rule_footprint", toplot))  
    @warn "footprint not implemented"
    

  end
   

  if any(isequal.("harvest_control_rule", toplot))  

    r = vec( Array(res[:, Symbol("r"), :]) )
    K = vec( Array(res[:, Symbol("K"), :]) ) 
    (msy, bmsy, fmsy) = logistic_discrete_reference_points(r, K)

    pl = hline!(pl, fmsy[ss]; alpha=0.01, color=:lightgray )
    pl = hline!(pl, [mean(fmsy)];  alpha=0.6, color=:darkgray, lw=5 )
    pl = hline!(pl, [quantile(fmsy, 0.975)];  alpha=0.5, color=:gray, lw=2, line=:dash )
    pl = hline!(pl, [quantile(fmsy, 0.025)];  alpha=0.5, color=:gray, lw=2, line=:dash )
  
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
  
    nt = length(prediction_time_ss)
    colours = get(ColorSchemes.tab20c, 1:nt, :extrema )[rand(1:nt, nt)]
  
    # scatter!( fb, FM ;  alpha=0.3, color=colours, markersize=4, markerstrokewidth=0)
    fb = bio[1:length(prediction_time_ss),:]
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
      series_annotations = text.(trunc.(Int, prediction_time_ss), :top, :left, pointsize=8) )

    ub = max( quantile(K, 0.75), maximum( fb_mean ), maximum(fmsy) ) * 1.05
    pl = plot!(pl; legend=false, xlim=(0, ub ), ylim=(0, maximum(fm_mean ) * 1.05  ) )
    # TODO # add predictions ???
  
  end
   
 

  if any(isequal.("prediction_space", toplot))  
    # b(t+1) vs b(t) * K / q1
    
    K = vec( Array(res[:, Symbol("K"), :]) )[ss] 
    q1 = vec( Array(res[:, Symbol("q1"), :]) )[ss] 

    t1 = 2:length(survey_time)
    t0 = 1:(length(survey_time)-1)

    b1 = bio[t1,ss] ./ K' .* q1'  # b(t+1)
    b0 = S[t0] ./ K' .* q1'   # b(t)

    nt = length(survey_time)
    colours = get(ColorSchemes.tab20c, 1:nt, :extrema )[rand(1:(nt-1), (nt-1) )]
    pl = scatter!(pl, b0, b1;  alpha=0.2, color=colours, markersize=2.5, markerstrokewidth=0)

  end


  if any(isequal.("state_space", toplot))  
    # b(t+1) vs b(t) * K / q1
    K = vec( Array(res[:, Symbol("K"), :]) )[ss]
    q1 = vec( Array(res[:, Symbol("q1"), :]) )[ss]

    t1 = 2:length(survey_time)
    t0 = 1:(length(survey_time)-1)

    b = bio[:,ss]  ./ K' .* q1' 
    b1 = b[ t1,:]  # b(t+1)
    b0 = b[ t0,:]  # b(t)

    nsims = size(bio)[2]
    ss = rand(1:nsims, 1000)  # sample index
  
    nt = length(survey_time)
    colours = get(ColorSchemes.tab20c, 1:nt, :extrema )[rand(1:(nt-1), (nt-1) )]
    pl = scatter!(pl, b0, b1;  alpha=0.2, color=colours, markersize=2.5, markerstrokewidth=0)
  end

  return(pl)

end


function plot_prior_posterior( vn, prior, posterior; bw=0.02 )
  pri =  vec(collect( prior[:,Symbol(vn),:] ))
  pos =  vec(collect( posterior[:,Symbol(vn),:] ))
  pl = plot(ylabel="Density", xlabel=vn ) 
  pl = density!(pl, pri,  fill=true, color = :slateblue, fillalpha=0.25, bandwidth = bw, lw=0, label="Prior")
  pl = density!(pl, pos,  fill=true, color = :purple, fillalpha=0.5, bandwidth = bw, lw=0, label="Posterior")
  return(pl)
end

  
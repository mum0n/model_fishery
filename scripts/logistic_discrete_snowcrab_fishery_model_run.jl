

# script to be run at top (main level) for a given combination of aulab and  model_variation
# only by "logistic_discrete_snowcrab_fishery_model.md"

   

# set up data and parameters
PM, M = logistic_discrete_data( yrs, aulab, fishery_data, model_variation )

Random.seed!(year_assessment);

# Turing.setprogress!(false);

rand(M) # above instantiates model as M: make sure M works
     
# MCMC samples

# Logging.disable_logging(Logging.Warn) # or e.g. Logging.Info

# NOTE:  make sure to check autocorrelation where required and thin


# adjust parameters for specific model_variation and get model setup
# model-specific sampler tweaks/over-rides
if model_variation=="logistic_discrete_historical"   # pre-2022, mimic STAN defaults
    max_depth =  7   
end

if model_variation=="logistic_discrete_basic"
    # same as historical but normalize to reduce influce of large magnitudes and faster convergence
    target_acceptance_rate, max_depth, init_ϵ = 0.65, 7, 0.25
end

 
prior = sample( M, Prior(), n_samples ) # 1- 10 sec

# Enzyme.API.runtimeActivity!(true);
# Enzyme.API.typeWarning!(false);

# by default use NUTS sampler ... SMC is another good option if NUTS is too slow
# turing_sampler = Turing.NUTS(n_samples, target_acceptance_rate; max_depth=max_depth, init_ϵ=init_ϵ )

turing_sampler = Turing.NUTS( n_samples, target_acceptance_rate; max_depth=max_depth, init_ϵ=init_ϵ, adtype = ADTypes.AutoForwardDiff() )

res = sample( M, turing_sampler, MCMCThreads(), n_samples, n_chains ; init_ϵ=init_ϵ, drop_warmup=true, progress=true ) # sample in parallel
 
## Save results: to (outputs_dir) as a hdf5  .. can also read back in R as:  h5read( res_fn, "res" )
# directory location is created in environment
res_fn = joinpath( outputs_dir, string("results_turing", "_", aulab, ".hdf5" ) )   # already defined in *.environment.jl .. here as a reminder
@save res_fn res

#=  to reload a save file:

  @load res_fn res
  
=#
  
## Using model solutions, compute other variables of interest via posterior simulation and load into main memory:

# generate predicted trajectories and params of interest from posterior samples 
# n scaled, n unscaled, biomass of fb with and without fishing, model_traces, model_times 
m, num, bio, trace, trace_bio, trace_time = fishery_model_predictions(res ) 
Fkt, FR, FM = fishery_model_mortality()   # fishing (kt), relative Fishing mortality, instantaneous fishing mortality:

# quick check of fits and solutions:
# fishery_model_plot( toplot=("survey", "fishing"), alphav=0.025 )

# save a few data files as semicolon-delimited CSV's for use outside Julia
summary_fn = joinpath( outputs_dir, string("results_turing", "_", aulab, "_summary", ".csv" ) )  
CSV.write( summary_fn,  summarize( res ), delim=";" )  # use semicolon as , also used in parm names
  
bio_fn1 = joinpath( outputs_dir, string("results_turing", "_", aulab, "_bio_fishing", ".csv" ) )  
CSV.write( bio_fn1,  DataFrame(bio[:,:], :auto), delim=";" )  # use semicolon as , also used in parm names
  
fm_fn = joinpath( outputs_dir, string("results_turing", "_", aulab, "_fm", ".csv" ) )  
CSV.write( fm_fn,  DataFrame( FM, :auto), delim=";" )  # use semicolon as , also used in parm names
 
## Figures

# pl = plot(pl, ylim=(0, 2.5))  # to rescale axes

# annual snapshots of biomass (kt) 
pl = fishery_model_plot( toplot=("survey", "fishing" ) )
# pl = plot(pl, ylim=(0, 4.0))
savefig(pl, joinpath( outputs_dir, string("plot_predictions_", aulab, ".", outformat ) )  ) 

# plot fishing mortality
pl = fishery_model_plot( toplot="fishing_mortality" )
pl = plot(pl, ylim=(0, 0.6))
savefig(pl, joinpath( outputs_dir, string("plot_fishing_mortality_", aulab, ".", outformat ) )  ) 

# HCR plot
pl = fishery_model_plot( toplot="harvest_control_rule", n_sample=1000 ) #, alphav=0.01 )  # hcr
pl = plot(pl, ylim=(0, 0.5))
savefig(pl, joinpath( outputs_dir, string("plot_hcr_", aulab, ".", outformat ) )  ) 
 
# grey is prior, purple is posterior 
# prior = sample( M, Prior(), n_samples)
L = truncated( Normal( PM.r[1], PM.r[2]), PM.r[3], PM.r[4])
pl = plot_prior_posterior( "r", prior, res )
pl = plot(pl, x->pdf(L, x))
save(joinpath( outputs_dir, string("plot_prior_r_", aulab, ".", outformat ) ), pl  )  

L = truncated( Normal( PM.K[1], PM.K[2]), PM.K[3], PM.K[4]) 
pl = plot_prior_posterior( "K", prior, res, bw=0.04)
pl = plot(pl, x->pdf(L, x))
save(joinpath( outputs_dir, string("plot_prior_K_", aulab, ".", outformat ) ), pl  )

L = truncated( Normal( PM.q1[1], PM.q1[2]), PM.q1[3], PM.q1[4] )    
pl = plot_prior_posterior( "q1", prior, res)
pl = plot(pl, x->pdf(L, x))
save(joinpath( outputs_dir, string("plot_prior_q1_", aulab, ".", outformat ) ), pl  )


L = truncated( Normal( PM.bpsd[1], PM.bpsd[2]), PM.bpsd[3], PM.bpsd[4] )
pl = plot_prior_posterior( "bpsd", prior, res)
pl = plot(pl, x->pdf(L, x))
save(joinpath( outputs_dir, string("plot_prior_bpsd_", aulab, ".", outformat ) ), pl  )


L = truncated( Normal( PM.bosd[1], PM.bosd[2]), PM.bosd[3], PM.bosd[4] ) 
pl = plot_prior_posterior( "bosd", prior, res)
pl = plot(pl, x->pdf(L, x))
save(joinpath( outputs_dir, string("plot_prior_bosd_", aulab, ".", outformat ) ), pl  )


pl = fishery_model_plot( toplot="state_space", n_sample=1000 ) #, alphav=0.01 )  # hcr
savefig(pl, joinpath( outputs_dir, string("state_space_", aulab, ".", outformat ) )  )

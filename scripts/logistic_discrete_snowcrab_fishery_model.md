
# Logistic discrete biomass dynamics models

This model form is used operationally for the snow crab fisheries of the Maritimes Region of Canada. It is implemented in Julia and Turing library to permit Bayesian parameter estimation. 

## Background

Population dynamics can be studied in many ways. Arguably, the simplest is to see it as a single system, such as the total number of individuals. Such single component models have a long tradition (Verhulst 1845, McKendrick & Pai 1912, Pearl & Reed 1920, Lotka 1925, Bacaër 2011). This is often called a *phenomenological* perspective in that it abstracts a multitude of processes and interactions into very few externally facing parameters and one state variable (population size). Often used for such cases, the logistic model specifies dynamics with two parameters: a maximum specific rate of change, and some upper limit of system state. Their utility and properties are well known and used in many domains.

Verhulst's original formulation of the logistic equation was used to model Belgium's population growth (Verhulst 1845, Bacaër 2011). His arguments were also phenomenological, alluding to Malthusian geometric (exponential) growth with some upper limit. The intrinsic rate of increase (r) represented the exponential rate of increase, while carrying capacity (K) represented this upper limit. Subsequently, McKendrick & Pai (1912) used it for bacterial growth and then (Pearl & Reed 1920, Lotka 1925) used it for the US population. It continues to be used today in numerous fields due to its underlying sigmoid shape (Beverton & Holt 1957, Pianka 1970, Bohner & Chieochan 2013).

Looking for a stronger justification, Lotka (1925) argued that the logistic equation can be seen as a 2nd order Taylor (Maclaurin) series approximation of any general function $\frac{dN}{dt}=G(N)$:

$$\frac{dN}{dt} = rN -(rK^{-1}) N^{2} + O(N^{h}),$$  

where, $O(N^{h})$ are higher order terms, $h>2.$

Pianka (1970, Fig. 9.5) suggested a slightly different approach. He begins with birth rate $\beta = \beta_0 - eN$ and death rate $\delta = \delta_0 + fN$ that are each linear functions of abundance. At equilibrium $N=K$, the lines intersect with $\beta = \delta$: 

$$\beta_0 - eK = \delta_0 + fK$$

$$r_0 = \beta_0 - \delta_0 = (e-f)K.$$

Away from equilibrium, after some algebra, we obtain the logistic equation:

$$\frac{dN}{dt} = rN = [(\beta_0 - eN) - (\delta_0 + fN)] N = r_0 N - r_0 N^2.$$

Interestingly, the sigmoid Beverton-Holt difference equation has also been shown to be related to the logistic equation (Beverton & Holt 1957, Bohner & Chieochan 2013):

$${n_{t+1}=\frac{\nu Kn_t}{{K_t + {({\nu-1})}}n_t}}.$$

It is discrete in time and also sigmoid in shape with asymptotic limit (carrying capacity that varies with time):

$K_{t} > 0$ and growth rate $\nu>1$, for number $n_{t}>0$. 

If we define: $r={{({\nu-1})}/\nu}$ then $\nu={1/{({1-r})}}$, 

then substitution and solving for $n_{t}$ and then ${\Delta n}_{t}$ gives:

$${{\Delta n}_{t}=r}n_{t+1}{{({K_{t}-n_{t+1}})}/{({{K_{t}-r}n_{t+1}})}},$$

which simplifies further to the discrete logistic difference equation:

$${{\Delta n}_{t}=r}n_{t+1}{({1-{n_{t}/K_{t}}})}.$$ 

In the limit as $\Delta t\rightarrow0$, this becomes the differential form of the
logistic equation:

$${{{dn}/{dt}}={rn}}{({1-{n/K_{t}}})}.$$


Finally, another similar model, the "logistic discrete map", could be used for parameter estimation. The "discrete logistic map" is often
represented as a Euler approximation of the logistic differential equation:

$${\frac{{dZ}}{{dt}}=r}Z{({{1-Z}K^{-1}})}.$$

Normalization of $Z$ to $z={Z/K}$ allows simplification to:

$$\frac{dz}{dt}=r z ( 1-z).$$ 

The Euler approximation for a small increment of time, ${\tau_{t+1}={\tau_{t}+\Delta}}\tau$, $z$ in the next time interval is estimated as:

$$z_{t+1} = z_{t} + \Delta \tau\cdot\frac{dz}{dt} z_{t},$$ 

which upon substitution of the differential equation and simplification gives:

$$z_{t+1} = ( 1 + \Delta \tau r) z_{t} - \Delta \tau r z_{t}^{2}.$$

The identities: ${\rho={1+\Delta}}\tau r$ and ${n_{t}=\frac{\Delta\tau rz_{t}}{\rho}}$ help to simplify further:

$$n_{t+1} = \frac{ \Delta \tau r z_{t+1}}{\rho}.$$ 

Substitution of the identity of $z_{t+1}$ to gives: 

$$n_{t+1} = \frac{\Delta \tau r (( 1 + \Delta \tau r) z_{t}-\Delta \tau r{z_t}^2)}{\rho} = \rho n_{t} (1-n_{t}).$$

Note that $\rho$ is ${1+\textnormal{rescaled}}{(r)}$, 

where we recall that it means a fractional change relative to the previous *time step*. Similarly, $n$ are rescaled values of $z$, more meaningful with a time step of $\Delta t$. Importantly, this derivation is for infinitesimal changes in time. If the increment in time is relatively large, such as $\Delta t = 1 \; \textnormal{year}$, as it is commonly used in fisheries applications, the "infinitesimal" change assumption is not biologically reasonable.

### References

Bacaër N (2011) Verhulst and the logistic equation (1838). In: *A Short History of Mathematical Population Dynamics*. Bacaër N (ed), Springer, London, p 35--39.

Bertalanffy L von (1950) The theory of open systems in physics and biology. Science 111:23--29.

Beverton R, Holt S (1957) On the dynamics of exploited fish populations, Volume 19 of Fishery investigations (Great Britain, Ministry of Agriculture, Fisheries, and Food). HM Station Off Lond.

Bohner M, Chieochan R (2013) The Beverton--Holt q-difference equation. J Biol Dyn 7:86--95. 

Bolker BM (2008) Ecological Models and Data in R. Princeton University Press.

Lotka AJ (1920) Analytical note on certain rhythmic relations in organic systems. Proc Natl Acad Sci 6:410--415. 

Lotka AJ (1925) The elements of physical biology. Williams and Wilkins.

McKendrick AG, Pai MK (1912) XLV.---The Rate of Multiplication of
Micro-organisms: A Mathematical Study. Proc R Soc Edinb
31:649--655.
 
Pearl R, Reed LJ (1920) On the Rate of Growth of the Population of the United States since 1790 and Its Mathematical Representation. Proc Natl Acad Sci U S A 6:275--288.

Pianka ER (1970) On r- and K-selection. Am Nat 104:592--597.

Quinn TJ (2003) Ruminations on the Development and Future of Population Dynamics Models in Fisheries. Nat Resour Model 16:341--392.

Schaefer MB (1954) Some Aspects of the Dynamics of Populations Important to the Management of the Commercial Marine Fisheries.

Verhulst PF (1845) Resherches mathematiques sur la loi d'accroissement de la population. Nouv Memoires Acad R Sci 18:1--41.

Wangersky PJ (1978) Lotka-volterra population models. Annu Rev Ecol Syst 9:189--218.


## Fishery application

Also known as surplus production models, these are simple models that treat the overall population as one homogenous entity with a simple perspective: increase and decrease are governed as a first order input-output processes that when combined provide a second order model. It is called "surplus production" historically as the perspective was that a population produces more than is necessary to maintain/replace itself (due to natural mortality). This "surplus production" can be diverted to human consumption/exploitation/harvest. From this perspective, maximum surplus production was the target, commonly known as "Maximum Sustainable Yield" (MSY).   

Further abstraction is made by focussing upon biomass $B$ rather than numbers $N$:

$$B_{t+1} = B_t + g(B_t) - C_t,$$

where the term $g(B_t)$ is the "surplus production" as a function of biomass at time $t$. 

If:

- $g(B_t) = rB (1-B/K)$ ; then it is a "Schaefer" 1954 form

- $g(B_t) = rB(1-log(B)/log(K))$ ; then it is a "Fox" 1970 form

- $g(B_t) = \frac{r}{p} \times B(1-(B/K)^p)$ ; then it is known as a "Pella-Tomlinson" 1969 form. 

Note that the Pella-Tomlinson form == Schaefer form when $p=1$
and the Pella-Tomlinson form $\rightarrow$ Fox form as $p \rightarrow 0$.


Here, $q$ is called a catachability coefficient, $r$ is the intrinsic rate of increase
and $K$ is the carrying capacity (or ~ sometimes the average biomass prior to exploitation). Catachability scales some index of abundance I to the "true" biomass, B: $I_t = q B_t$. 


# Methods of parameter estimation

Reference: Hilborn & Walters (1992) - Chapter 8

The population process discretized on an annual basis t, is: 

biomass(t+1) = biomass(t) + recruitment(t) + growth(t) - catch(t) - natural mortality(t).

Usually, growth and reproduction are grouped together as production, and:

surplus production = production - natural mortality

that is, in the absence of fishing. This is the amount that can be caught while still maintaining biomass at equilibrium:

B(t+1) = B(t) + S(t) - C(t) 

or,

S(t) = B(t+1) - B(t) + C(t).


In Schaefer's formulation, $dB/dt = rB(1- B/K) - C$. So, $C(t) = q * E(t) * B(t)$ ; where q = "catchability", E is fishing effort and:

$O \approx U = C/E = qB$,

where O = observed index of abundance, and U is catch per unit effort or any observed index of abundance. Alternatively: $B = C/E/q$ 

Pella-Thomlinson form adds an exponent (m) to the biomass term (B -> B^m) 
such that if m=2 it becomes the Shaeffer form. 

Estimation possible via:

  - equilibrium assumptions/methods (unstable)
  
  - Regression (Schnute 1977; erratic) 
  
  - Observation error or "time-series method" (recommended by Hilborn & Walters), and currently known as state space models.

Using cpue or biomass index: 

$B(t) = B(t-1) + r * B(t-1) * ( 1 - B  (t-1) / K ) - C(t-1)$

$B(t) = B(t-1) * ( 1 + r * ( 1- B(t-1)/K )) - C(t-1)$ 

$O(t) = q * B(t)$

$e(t) = (O(t) - E[ O(t) ] )^2$

where E[.] is the expected/estimated operator, of CPUE or other abundance metric (eg survey), "observation error".

minimize( sum(e(t)) ) in least squares assumption.

Or, using effort only

$C(t) = B(t) * q * E(t)$

$e(t) = (C(t) - E[C(t)] ) ^2$

B(0) must be estimated or one can assume:

$B(1) = C(t) / ( E(1) * q )$ ; or 

$B(1) = Bmean$, a running average of the start of the C(1..5, etc); or 

$B(0) == K$, if starting with an unfished population .. 

etc ...

---
 

## The main modelling tools used:
  
  * [Julia](https://julialang.org/) 
  * [Turing](https://turing.ml/stable/) 
  * Other libraries used can be found in the set up environment (logistic_discrete_environment.jl)  

---

## Data requirements (Maritimes snow crab): 

If data files have not already been created, then do so. Data files are generated by:

  - [bio.snowcrab::03.snowcrab_index_carstm.r](https://github.com/jae0/bio.snowcrab/blob/f2e0722de0bfa08611caa775cf6d912807bd4d91/inst/scripts/03.biomass_index_sizestructured_carstm.r)

and then copy the files manually to appropriate locations or run:

```r

  # R-code snippet show data formatting step
  source( file.path( code_root, "bio_startup.R" )  )
  loadfunctions("bio.snowcrab")
  year.assessment = 2023

  modeldir = file.path( homedir, "bio.data", "bio.snowcrab", "modelled" )
  carstm_model_label = "default_fb"
 
  # flag to switch between production and testing data 
  # run_is_operational = TRUE  
  run_is_operational = FALSE

  fishery_model_data_inputs( 
    year.assessment=year.assessment,   
    type="biomass_dynamics", 
    snowcrab_filter_class="fb",
    modeldir= ifelse( run_is_operational, pN$modeldir, file.path( homedir, "projects", "dynamical_model", "snowcrab", "data" ) ),
    carstm_model_label=carstm_model_label,
    for_julia=TRUE,
    fishery_model_label="turing1"  
  ) 

```

The remaining code are in Julia.


---

## Define key directories and main run level options: 

These directories need to be defined. 

```julia

using DrWatson
quickactivate( joinpath("projects", "model_fishery") )
 
# load libs and local functions
include( scriptsdir( "logistic_discrete_environment.jl" ))

# install_required_packages(pkgs)  # in case you need to install packages

# run-time parameters 
year_assessment = 2023   # <<<<<<<<-- change
yrs = 1999:year_assessment  


# model-type
model_variation = "logistic_discrete_historical"   
 
  #= currently implemented models:
    model_variation = "logistic_discrete_historical"      # pre-2022 ("historical") method  :: ~ 1 hr # pre-2022, no normalization, q-based observation model
    model_variation = "logistic_discrete_basic"           # q catchability only for observation model
    model_variation = "logistic_discrete"                 # q and intercept for observation model
    model_variation = "logistic_discrete_map"             # logistic map ... more extreme fluctuations
  =#
  

aulab ="choose one from below"    

# area to model .. choose a region of interest:
  #= currently implemented regions:
    aulab ="cfanorth"    
    aulab ="cfasouth"    
    aulab ="cfa4x"       
  =#


# input/output

# outputs_dir = projectdir( homedir(), "bio.data", "bio.snowcrab",  "fishery_model", "outputs", string(year_assessment), model_variation )
outputs_dir = projectdir( "outputs", string(year_assessment), model_variation )
mkpath(outputs_dir)
print( "outputs_dir: ", outputs_dir, "\n\n" )

# mcmc save file name and location
res_fn = joinpath( outputs_dir, string("results_turing", "_", aulab, ".hdf5" ) )  
print( "results file:",  res_fn )
 

# fishery_data_dir = joinpath( project_directory, "data", "default_fb" )  # remember to save/copy data to this location (manually or in R)
fishery_data_dir = joinpath( homedir(), "bio.data", "bio.snowcrab", "modelled", "default_fb"   )

fndat  = joinpath( fishery_data_dir, "biodyn_biomass.RData" )

fishery_data = load( fndat, convert=true)
  
  #= 
    # alternatively, if running manually, run R-code that creates local RData file with required data
    # run in R externally or from within julia or from within julia:
    using RCall
    # type $ in Julia's command prompt starts an R session.
    # .. run below
    # type <backspace> to escape back to julia

    #```r
    source( file.path( code_root, "bio_startup.R" )  )
    require(bio.snowcrab)   # loadfunctions("bio.snowcrab")
    fishery_model_data_inputs( year.assessment=2021, type="biomass_dynamics",  for_julia=TRUE, time_resolution=1/12)
    #```

    # then back in Julia, fetch data into julia's workspace (replace fndat with the  filenane printed above )
    @rget Y
    @rget Kmu
    @rget removals
    @rget ty
 
    plot(  Y[:,:yrs], Y[:,:cfasouth] )    # example line plots
  =#

# include( scriptsdir( "logistic_discrete_environment.jl" ))

# set up data and parameters
PM, M = logistic_discrete_data( yrs, aulab, fishery_data, model_variation )

Random.seed!(year_assessment);

# Turing.setprogress!(false);

rand(M) # above instantiates model as M: make sure M works
     
# MCMC samples

# Logging.disable_logging(Logging.Warn) # or e.g. Logging.Info

# NOTE:  make sure to check autocorrelation where required and thin

n_adapts, n_samples, n_chains, thin = 10000, 40000, 4, 0
target_acceptance_rate, max_depth, init_ϵ = 0.65, 8, 0.01

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
savefig(pl, joinpath( outputs_dir, string("plot_predictions_", aulab, ".png") )  )
savefig(pl, joinpath( outputs_dir, string("plot_predictions_", aulab, ".pdf") )  )
savefig(pl, joinpath( outputs_dir, string("plot_predictions_", aulab, ".svg") )  )

# plot fishing mortality
pl = fishery_model_plot( toplot="fishing_mortality" )
pl = plot(pl, ylim=(0, 0.6))
savefig(pl, joinpath( outputs_dir, string("plot_fishing_mortality_", aulab, ".png") )  )
savefig(pl, joinpath( outputs_dir, string("plot_fishing_mortality_", aulab, ".pdf") )  )
savefig(pl, joinpath( outputs_dir, string("plot_fishing_mortality_", aulab, ".svg") )  )

# HCR plot
pl = fishery_model_plot( toplot="harvest_control_rule", n_sample=1000 ) #, alphav=0.01 )  # hcr
pl = plot(pl, ylim=(0, 0.5))
savefig(pl, joinpath( outputs_dir, string("plot_hcr_", aulab, ".png") )  )
savefig(pl, joinpath( outputs_dir, string("plot_hcr_", aulab, ".pdf") )  )
savefig(pl, joinpath( outputs_dir, string("plot_hcr_", aulab, ".svg") )  )
 
# grey is prior, purple is posterior 
# prior = sample( M, Prior(), n_samples)
L = truncated( Normal( PM.r[1], PM.r[2]), PM.r[3], PM.r[4])
pl = plot_prior_posterior( "r", prior, res )
pl = plot(pl, x->pdf(L, x))
save(joinpath( outputs_dir, string("plot_prior_r_", aulab, ".png") ), pl  )
save(joinpath( outputs_dir, string("plot_prior_r_", aulab, ".pdf") ), pl  )
save(joinpath( outputs_dir, string("plot_prior_r_", aulab, ".svg") ), pl  )

L = truncated( Normal( PM.K[1], PM.K[2]), PM.K[3], PM.K[4]) 
pl = plot_prior_posterior( "K", prior, res, bw=0.04)
pl = plot(pl, x->pdf(L, x))
save(joinpath( outputs_dir, string("plot_prior_K_", aulab, ".png") ), pl  )
save(joinpath( outputs_dir, string("plot_prior_K_", aulab, ".pdf") ), pl  )
save(joinpath( outputs_dir, string("plot_prior_K_", aulab, ".svg") ), pl  )

L = truncated( Normal( PM.q1[1], PM.q1[2]), PM.q1[3], PM.q1[4] )    
pl = plot_prior_posterior( "q1", prior, res)
pl = plot(pl, x->pdf(L, x))
save(joinpath( outputs_dir, string("plot_prior_q1_", aulab, ".png") ), pl  )
save(joinpath( outputs_dir, string("plot_prior_q1_", aulab, ".pdf") ), pl  )
save(joinpath( outputs_dir, string("plot_prior_q1_", aulab, ".svg") ), pl  )

L = truncated( Normal( PM.bpsd[1], PM.bpsd[2]), PM.bpsd[3], PM.bpsd[4] )
pl = plot_prior_posterior( "bpsd", prior, res)
pl = plot(pl, x->pdf(L, x))
save(joinpath( outputs_dir, string("plot_prior_bpsd_", aulab, ".png") ), pl  )
save(joinpath( outputs_dir, string("plot_prior_bpsd_", aulab, ".pdf") ), pl  )
save(joinpath( outputs_dir, string("plot_prior_bpsd_", aulab, ".svg") ), pl  )

L = truncated( Normal( PM.bosd[1], PM.bosd[2]), PM.bosd[3], PM.bosd[4] ) 
pl = plot_prior_posterior( "bosd", prior, res)
pl = plot(pl, x->pdf(L, x))
save(joinpath( outputs_dir, string("plot_prior_bosd_", aulab, ".png") ), pl  )
save(joinpath( outputs_dir, string("plot_prior_bosd_", aulab, ".pdf") ), pl  )
save(joinpath( outputs_dir, string("plot_prior_bosd_", aulab, ".svg") ), pl  )


pl = fishery_model_plot( toplot="state_space", n_sample=1000 ) #, alphav=0.01 )  # hcr
savefig(pl, joinpath( outputs_dir, string("state_space_", aulab, ".png") )  )
savefig(pl, joinpath( outputs_dir, string("state_space_", aulab, ".pdf") )  )
savefig(pl, joinpath( outputs_dir, string("state_space_", aulab, ".svg") )  )

```

# End 

  
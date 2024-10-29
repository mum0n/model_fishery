
# Age-structured cohort-based population  

Age-structured that track cohorts (age-groups) over time is/was used commonly in fisheries applications to estimate population structure/demographics from fisheries catch data as that is the simplest data to gather. Generally, catch is assumed to be known without error. Cohort populations were estimated via back-calculation/propagation/projection ("Cohort Analysis") with simple assumptions (catch + cohort) to help re-create a "virtual" population (Fry 1949, based on Derzhavin 1922):

$N_{a+1,y+1} = N_{a,y} - C_{a,y} - Z_{a,y}$,

for a given age $a$ and year $y$, in discrete form. 

**For each cohort**, refining of the above to instantaneous total mortality ($Z$) as a combination of natural ($M$) and fishing ($F$) mortality and catachability ($q$) is attributed to Ricker, Beverton, Holt, and Paloheimo. Gulland (1965) is attributed this "backwards" mode which provided more convergent solutions (stable); Murphy (1965) estimated fishing mortality at early ages in forward mode which was found to give more divergent solutions (unstable). As such backwards mode was generally more favoured.

Beginning with the Baranov catch equation:

$C_t = [F_t / (F_t + M) ] \; N_t (1 - e^{-Z_t} )$

and a first order, nonlinear, decay equation:

$N_{t+1} = N_t \; e^{-Z_t} \; \; <=> \; \; N_t = N_{t+1} \; e^{Z_t}$,

gives in terms of catch, "terminal" fishing mortality and natural mortality:

$N_t = C_t / [ ( F_t / Z_t ) (1 - e^{-Z_t}) ]$ 

at $y=t$ the terminal year and Terminal fishing mortality $F_t$ and total mortality $Z_t=F_t + M_{(t)}$, where M is usually a constant but can be variable with time or age.

A further simplification, to avoid nonlinear regression, is the use of a step function for fishing mortality that is assumed to occur in one instant at exactly 1/2 through a time cycle (usually 1 year) and background natural mortality at a constant rate; this approximation is attributed to Pope (1972):

$N_{t+0.5} = N_t \; e^{-M/2} - C_t  \;$

$N_{t+1} = N_{t+0.5} \; e^{-M/2}$

$\quad \quad \; = [N_t \; e^{-M/2} - C_t ] \; e^{-M/2}$

which essentially amounts to [natural mortality before fishing event - catch] * natural mortality of remainder. Generally results are more stable with Z < 1 and M < F. It simplifies to:

$N_{t+1}  = N_t \; e^{-M} - C_t \; e^{-M/2}$

Solving for $N_t$ gives:

$N_{t} = N_{t+1} \; e^M + C_t \; e^{M/2}$, 

and fishing mortality can be estimated as:

$F_t = \log( | N_t / N/{t+1}| ) - M$

Various methods of estimating M exist. A useful one is to express Z vs effort. The Y-intercept provides the magnitude of background natural mortality and the slope represents the catchability coefficient q.

Additional structure in the form of distributional constraints that were "tuned" to an independent index of abundance, were known as ADAPT-VPA's and with multiple indices, as Integrated Catch Analysis. With the use of priors (e.g., Bayesian modelling) continued evolution/refinement of numerical methods and assumptions on age and maturity could be made, see eg, Millar and Meyer (2000) CJFAS 57:43-50 (Sequential Population Analysis).

All these methods assume age structure is reasonably well represented and errors in ageing of animals and usually a fixed natural mortality assumption and that catch is observed without error, no movement and no discards. They are important sources of uncertainty.


### Good resources: 

https://people.uncw.edu/scharff/courses/458/Lecture%2013%20-%20cohort%20models.pdf



---

## The main modelling tools used here:
  
  * [Julia](https://julialang.org/) 
  * [Turing](https://turing.ml/stable/) 
  * Other important libraries are found in the set up "environment" files

---


## Define key directories and main run level options and create input data: 

To begin, define starting directories and the example data (source: FAO tutorial's VPA data for Baltic Cod: Div 25-32): 


```julia

using DrWatson
quickactivate( joinpath("projects", "model_fishery") )
 
# load libs and local functions
include( scriptsdir( "age_structured_cohort_environment.jl" ))

# install_required_packages(pkgs)  # in case you need to install packages

ages, yrs, C, A_yrs, A, body_weight_kg, maturity_ogive = data_baltic_cod_fao()

  #=
    C = catch, numbers at age x year
    A = abundance index, numbers age x year
    the remaining variables should be self-explanitory 
  =#


```
 
## Cohort Analysis

Cohort Analysis is the core concept of most age-structured methods (attributed to Pope): 
- make assumptions on mortality (fishing and natural), usually of older animals (better captured) and plus groups 
  - Natural mortality assumed at 0.2.
  - Fishing_mortality at age is assumed constant for 7 and 8+.
 - Work backwards using numbers of these groups 
 - Exceptions: Sequential Pop. Anal., Stock synthesis (forward methods).
 

```julia
 
# dimensions and initial conditions
na, ny = size(C)

a = 1:(na - 2 )    # age indices for main age groups
ap = [ na-1, na ]  # age indices of the plus age classes (fixed)

y = 1:(ny-1)      # year indices from 1 to next to last year

# create F0 based on following conditions
Fyp = repeat([0.55], ny)  # fishing mortality across yr for the plus groups 
Fat = [ 0.1, 0.2, 0.4, 0.5, 0.5, 0.5, 0.5 ] # Fishing mortality across age for the terminal year
F0 = fishing_mortality_setup_vpa(na, ny, ap, Fyp, Fat)

# guess of natural mortality
M = repeat([0.2], na)  

# re-estimate N from F0 and M 
N0 = N_setup_vpa(na, ny, ap, C, F0, M)

# deterministically back-solve with above starting conditions
N = cohort_analysis( N0, C, M; Mtype="backward" ) 

B, Bmat, Btot, Bssb, pl = compute_biomass( N, body_weight_kg[:,1], maturity_ogive )
pl

F, Fi, pl = compute_fishing_mortality(F0, N, M, a; i=1:5)  # update F with finalized N's
pl


```

## Virtual Population Analysis (VPA)

VPA's tune Cohort Analysis results by assuming separability of fishing mortality by age and year and estimating these parameters.

Tuning is/was accomplished by least squares or Maximum-likelihood methods (respecting distributional forms of data and parameters). In the basic VPA, fishery independent estimates (A) are tuned to abundance estimates from Cohort Analysis by fitting a catachability coefficient (q) and Fishing mortality estimates of the most recent year (Fyp) and oldest groups (Fat).
 
Tuning to surveys and other sources of information are called ADAPT and Integrated Catch Analysis. 

Tuning in the reverse direction is called Sequential Population Analysis.


```julia

commonyears = collect(intersect( yrs, A_yrs))
yi = findall( x-> x in commonyears, yrs)

# objective function:  add error from CPUE indices
Ndiff = log.( N[:, yi] ) .- log.( A )
q0 = q = mean( Ndiff, dims=2 ) # across age  
objfun = Ndiff .- q
 
q0 = vec(q)  # use as prior for Bayesian methods (below)
 
```

Wrap the above into a Turing model to flexibly estimate parameters:

```julia

Fyp0 = repeat([0.55], ny)
Fat0 = [ 0.1, 0.2, 0.4, 0.5, 0.5, 0.5, 0.5 ] 
 
Random.seed!(1);

Fishing_mortality_assumption = "nonseparable"
Fishing_mortality_assumption = "separable"

fmod = Virtual_Population_Analysis( Fishing_mortality_assumption, M, C  )  # not exactly a VPA as it is not a recursive form
fmod = Adaptive_Virtual_Population_Analysis( Fishing_mortality_assumption, M )
fmod = Integrated_Catch_Analysis( Fishing_mortality_assumption, M, C  )
fmod = Sequential_Population_Analysis( Fishing_mortality_assumption, M, C  ) # untested

# above instantiates model as fmod: make sure fmod works
rand(fmod) 
 
# find estimates (various methods)
res = optimize(fmod, MLE(), NelderMead(); autodiff = :forward)  # MLE estimate .. does not converge

Optim.Options(iterations=100)  # ParticleSwarm does not check for convrgence ... control with Optim.Options(iterations=...)  
res = optimize(fmod, MAP(), ParticleSwarm() ; autodiff = :forward)  # MAP estimate... not converging ... 
  
    # extraction:
    Fyp_est = [ res.values[Symbol("Fyp[$i]")] for i in 1:ny ]
    Fat_est = [ res.values[Symbol("Fat[$i]")] for i in 1:na ]
    q_est = [ res.values[Symbol("q[$i]")] for i in 1:na ] # only in Adaptive_Virtual_Population_Analysis and ica
    S_est =  [ res.values[Symbol("S[$i, $j]")] for i in 1:na for j in 1:ny ]

    N = reshape( S_est, na, ny )

    B, Bmat, Btot, Bssb, pl = compute_biomass( N, body_weight_kg[:,1], maturity_ogive )
    pl
    
    F, Fi, pl = compute_fishing_mortality(F0, N, M, a; i=1:2)  # update F with finalized N's
    pl


# Generate a Variational Inference estimate...  Not working ? 2024-10-26 
using Turing: Variational
advi = ADVI(10, 1000, AutoForwardDiff() ) # samples_per_step (gradient estimation), max_iters (no gradient steps)
res = vi(fmod, advi )

n_samples = 100
samples = rand(res, n_samples);

    # extraction:
    Fyp_vi_ = [ res.values[Symbol("Fyp[$i]")] for i in 1:ny ]
    Fat_vi_ = [ res.values[Symbol("Fat[$i]")] for i in 1:na ]
    q_vi_ = [ res.values[Symbol("q[$i]")] for i in 1:na ]
    S_vi_ =  [ res.values[Symbol("S[$i, $j]")] for i in 1:na for j in 1:ny ]

    N = reshape( S_vi_, na, ny )

    F, Fi, pl = compute_fishing_mortality(F0, N, M, a; i=1:2)  # update F with finalized N's
    pl

    B, Bmat, Btot, Bssb, pl = compute_biomass( N, body_weight_kg[:,1], maturity_ogive )
    pl


# MCMC 

## Prior
n_samples = 100
prior = sample( fmod, Prior(), n_samples )  
showall(prior)

## SMC 
n_samples = 100
res = sample( fmod, Turing.SMC(),  n_samples   ) # sample in parallel


## NUTS
n_adapts, n_samples, n_chains, thin = 1000, 4000, 4, 0
target_acceptance_rate, max_depth, init_ϵ = 0.65, 8, 0.01

turing_sampler = Turing.NUTS( n_samples, target_acceptance_rate; max_depth=max_depth, init_ϵ=init_ϵ )

turing_sampler = Turing.NUTS(; adtype=AutoForwardDiff())

res = sample( fmod, turing_sampler,  n_samples   ) # sample in parallel

  # extraction:

    k = rand(1:nsims)  # nsims
    l = rand(1:nchains) # nchains
    k=1:n_samples
    l=1

    Fyp_est = [ res[k, Symbol("Fyp[$i]"), l] for i in 1:ny ]
    Fat_est = [ res[k, Symbol("Fat[$i]"), l] for i in 1:na ]
    q_est = [ res[k, Symbol("q[$i]"), l] for i in 1:na ]
    S_est =  [ res[k, Symbol("S[$i, $j]"), l] for i in 1:na for j in 1:ny ]

    S = vec( mean( hcat( S_est...), dims=1))

    N = reshape( S, na, ny )

    B, Bmat, Btot, Bssb, pl = compute_biomass( N, body_weight_kg[:,1], maturity_ogive )
    pl

    F, Fi, pl = compute_fishing_mortality(F0, N, M, a; i=1:2)  # update F with finalized N's
    pl

   

```

# End 

  
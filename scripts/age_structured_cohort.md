
# Age-structured cohort-based population  

Age-structured that track cohorts (age-groups) over time is/was used commonly in fisheries applications to estimate population structure/demographics from fisheries catch data as that is the simplest data to gather.  See [writeup here](../docs/cohort_models.md).

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

quickactivate( joinpath(homedir(), "projects", "model_fishery") )
 
# load libs and local functions
include( scriptsdir( "age_structured_cohort_environment.jl" ))

# install_required_packages(pkgs)  # in case you need to install packages

ages, yrs, C, A_yrs, A, body_weight_kg, maturity_ogive = data_baltic_cod_fao()

Cobs = Float64.(C) # used where Integer is not appropriate

  #=
    C = catch, numbers at age x year
    A = abundance index, numbers age x year
    the remaining variables should be self-explanatory 
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
iplus = [na-1, na]  # indices of plus size/age

# guess of natural mortality
M = 0.2   #  can be a vector or matrix too .. e.g., repeat([0.2], na)  

# create F based on following conditions
Fp = repeat([0.5], ny)  # fishing mortality across yr for the plus groups 
Ft = [ 0.1, 0.2, 0.4, 0.5, 0.5, 0.5, 0.5 ] # Fishing mortality across age for the terminal year
F = setup_F( Fp, Ft, type="nonseparable", iplus=iplus )

# estimate N from F and M 
N = setup_N( C, F, M; iplus=iplus )

# deterministically (without error processes) back-solve with above starting conditions
N = cohort_analysis!( N, C, M; iplus=iplus  ) 

F = compute_fishing_mortality!(F, N, M; iplus=iplus)  # update F with finalized N's
Fall, fsubset, pl = summary_fishing_mortality(F, subset=1:5)
pl

B = N .* body_weight_kg[:,1]
Bmat, Btot, Bssb, pl = compute_biomass( B, maturity_ogive )
pl

Cpred = catch_cohort_estimate( N, F, M ) 

  
```

## Virtual Population Analysis (VPA)

VPA's tune Cohort Analysis results by assuming separability of fishing mortality by age and year and estimating these parameters.

Tuning is/was accomplished by least squares or Maximum-likelihood methods (respecting distributional forms of data and parameters). In the basic VPA, fishery independent estimates (A) are tuned to abundance estimates from Cohort Analysis by fitting a catachability coefficient (q) and Fishing mortality estimates of the most recent year (Fp) and oldest groups (Ft).
 
Tuning to surveys and other sources of information are called ADAPT and Integrated Catch Analysis. 

Tuning in the reverse direction is called Sequential Population Analysis.

So, approach is to wrap Cohort Analysis with random initial F estimates:

```julia

Fp0 = repeat([0.5], ny)
Ft0 = [ 0.1, 0.2, 0.4, 0.5, 0.5, 0.5, 0.5 ] 
 
Random.seed!(1);

  #=
    test = rand(arraydist( LogNormal.( log.(Fp0), 0.25) ))
    test = rand(arraydist( LogNormal.( log.(Ft0), 0.25) ))

    histogram(test)
  =#


Fishing_mortality_assumption = "nonseparable"
Fishing_mortality_assumption = "separable"


  #=
    # implemented models .. add more using the same template
    fmod = Virtual_Population_Analysis_basic( Fishing_mortality_assumption, M, C , iplus )  # not exactly a classical VPA as Ft and Fp are random effects
    fmod = Virtual_Population_Analysis( Fishing_mortality_assumption, M, C , iplus ) 
    fmod = Adaptive_Virtual_Population_Analysis( Fishing_mortality_assumption, M, C , iplus)
    fmod = Integrated_Catch_Analysis( Fishing_mortality_assumption, M, C , iplus )
    fmod = Sequential_Population_Analysis( Fishing_mortality_assumption, M, C, iplus  ) # untested
  =#


 
# choose approach:

if approach =="Maximum_likelhood_basic"
  Optim.Options(iterations=1000)  # change options as required  
  fmod = Virtual_Population_Analysis_basic( Fishing_mortality_assumption, M, Cobs, iplus )  # not exactly a classical VPA as Ft and Fp are random effects
  res = optimize(fmod, MLE(), NelderMead())  # alt engines: BFGS(), Newton(), etc .. 
    Fp_hat = [ res.values[Symbol("Fp[$i]")] for i in 1:ny ]
    Ft_hat = [ res.values[Symbol("Ft[$i]")] for i in 1:na ]
    F = setup_F( Fp_hat, Ft_hat, type=Fishing_mortality_assumption, iplus=iplus )
    N = setup_N( C, F, M; iplus=iplus )
    N = cohort_analysis!( N, C, M; iplus=iplus  ) 
    F = compute_fishing_mortality!(F, N, M; iplus=iplus)  # update F with finalized N's
    Fall, fsubset, pl = summary_fishing_mortality(F, subset=1:5)
    pl
    B = N .* body_weight_kg[:,1]
    Bmat, Btot, Bssb, pl = compute_biomass( B, maturity_ogive )
    pl
    Cpred = catch_cohort_estimate( N, F, M ) # Note some negative predicted catches 
end


if approach =="Maximum_likelhood"
  Optim.Options(iterations=10000)  # change options as required  
  fmod = Virtual_Population_Analysis( Fishing_mortality_assumption, M, Cobs, iplus )  # not exactly a classical VPA as Ft and Fp are random effects
  res = optimize(fmod, MLE(), NelderMead() )  # alt engines: NelderMead(), BFGS(), Newton(), etc .. 
    # note addition of the data likelihood results in model divergence   
    Fp_hat = [ res.values[Symbol("Fp[$i]")] for i in 1:ny ]
    Ft_hat = [ res.values[Symbol("Ft[$i]")] for i in 1:na ]
    F = setup_F( Fp_hat, Ft_hat, type=Fishing_mortality_assumption, iplus=iplus )
    N = setup_N( C, F, M; iplus=iplus )
    N = cohort_analysis!( N, C, M; iplus=iplus  ) 
    F = compute_fishing_mortality!(F, N, M; iplus=iplus)  # update F with finalized N's
    Fall, fsubset, pl = summary_fishing_mortality(F, subset=1:5)
    pl
    B = N .* body_weight_kg[:,1]
    Bmat, Btot, Bssb, pl = compute_biomass( B, maturity_ogive )
    pl
    Cpred = catch_cohort_estimate( N, F, M ) # Note some negative predicted catches 
end



if approach =="Maximum_aposteriori"
  Optim.Options(iterations=10000)  # change options as required  
  fmod = Virtual_Population_Analysis( Fishing_mortality_assumption, M, Cobs, iplus )  # not exactly a classical VPA as Ft and Fp are random effects
  res = optimize(fmod, MAP(), NelderMead() )  # alt engines: NelderMead(), BFGS(), Newton(), etc .. 
    # note addition of the data likelihood results in model divergence   
    Fp_hat = [ res.values[Symbol("Fp[$i]")] for i in 1:ny ]
    Ft_hat = [ res.values[Symbol("Ft[$i]")] for i in 1:na ]
    F = setup_F( Fp_hat, Ft_hat, type=Fishing_mortality_assumption, iplus=iplus )
    N = setup_N( C, F, M; iplus=iplus )
    N = cohort_analysis!( N, C, M; iplus=iplus  ) 
    F = compute_fishing_mortality!(F, N, M; iplus=iplus)  # update F with finalized N's
        Fall, fsubset, pl = summary_fishing_mortality(F, subset=1:5)
    pl
    B = N .* body_weight_kg[:,1]
    Bmat, Btot, Bssb, pl = compute_biomass( B, maturity_ogive )
    pl
    Cpred = catch_cohort_estimate( N, F, M ) # Note some negative predicted catches 
end



q_hat = [ res.values[Symbol("q_hat[$i]")] for i in 1:na ]
S_hat = [ res.values[Symbol("S[$i, $j]")] for i in 1:na for j in 1:ny ]
N = reshape( S_hat, na, ny )

F = compute_fishing_mortality!(F, N, M; iplus=iplus)  # update F with finalized N's
Fall, fsubset, pl = summary_fishing_mortality(F, subset=1:5)
pl

B = N .* body_weight_kg[:,1]
Bmat, Btot, Bssb, pl = compute_biomass( B, maturity_ogive )
pl

Cpred = catch_cohort_estimate( N, F, M ) # Note some negative predicted catches 



# Generate a Variational Inference estimate...  Not working ? 2024-10-26 
using Turing: Variational
advi = ADVI(10, 1000, AutoForwardDiff() ) # samples_per_step (gradient estimation), max_iters (no gradient steps)
res = vi(fmod, advi )

n_samples = 100
samples = rand(res, n_samples);

  # extraction:
  Fp_hat = [ res.values[Symbol("Fp[$i]")] for i in 1:ny ]
  Ft_hat = [ res.values[Symbol("Ft[$i]")] for i in 1:na ]
  q_hat = [ res.values[Symbol("q[$i]")] for i in 1:na ]
  S_hat =  [ res.values[Symbol("S[$i, $j]")] for i in 1:na for j in 1:ny ]

  N = reshape( S_hat, na, ny )

  F = compute_fishing_mortality!(F, N, M; iplus=iplus)  # update F with finalized N's
  Fall, fsubset, pl = summary_fishing_mortality(F, subset=1:5)
  pl

  B = N .* body_weight_kg[:,1]
  Bmat, Btot, Bssb, pl = compute_biomass( B, maturity_ogive )
  pl

  Cpred = catch_cohort_estimate( N, F, M ) # Note some negative predicted catches 



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

    Fp_hat = [ res[k, Symbol("Fp[$i]"), l] for i in 1:ny ]
    Ft_hat = [ res[k, Symbol("Ft[$i]"), l] for i in 1:na ]
    q_hat = [ res[k, Symbol("q[$i]"), l] for i in 1:na ]
    S_hat =  [ res[k, Symbol("S[$i, $j]"), l] for i in 1:na for j in 1:ny ]

    S = vec( mean( hcat( S_hat...), dims=1))

    N = reshape( S, na, ny )

    F = compute_fishing_mortality!(F, N, M; iplus=iplus)  # update F with finalized N's
    Fall, fsubset, pl = summary_fishing_mortality(F, subset=1:5)
    pl

    B = N .* body_weight_kg[:,1]
    Bmat, Btot, Bssb, pl = compute_biomass( B, maturity_ogive )
    pl

    Cpred = catch_cohort_estimate( N, F, M ) # Note some negative predicted catches 

   

```

# End 

  
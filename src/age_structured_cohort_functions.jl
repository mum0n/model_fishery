

function data_baltic_cod_fao()
  
  # Data (example) come from FAO tutorial's VPA data for Baltic Cod: Div 25-32.
    
  ages = 2:8 # the last age is 8+
  na = length(ages)

  yrs = 1966:1996
  ny = length(yrs)

  A_yrs = 1981:1996  # S = survey
  A_ny = length(A_yrs)

# catch in fishery in numbers at age x year
  C = transpose( reshape( [ 
    29009, 36345, 30047, 37488, 24764, 20026, 32990, 28190, 23242, 21923, 8541, 6102, 40773, 50464, 20285, 21825, 37718, 36901, 22161, 32236, 11798, 41863, 21428, 4526, 16974, 5902, 7619, 8909, 9363, 7554, 10646, 91231, 132299, 107518, 102497, 80159, 47606, 69783, 64705, 93448, 135537, 62790, 42363, 72169, 143205, 111509, 52669, 119568, 108259, 103672, 62632, 40119, 55696, 84337, 59026, 30322, 39296, 14556, 17344, 36380, 22030, 30447, 34141, 50122, 70398, 49728, 51571, 38631, 42759, 42779, 49776, 70356, 86375, 46109, 46560, 76948, 160282, 101889, 81157, 122532, 144308, 74177, 58468, 44253, 52156, 69293, 39775, 31400, 13343, 8306, 23876, 30966, 26148, 11946, 14395, 23022, 23944, 22735, 21669, 18437, 18295, 16163, 18751, 38326, 35583, 22749, 26478, 62205, 83499, 53499, 46146, 83822, 57154, 43020, 23504, 17299, 23041, 20868, 13406, 6103, 3403, 10454, 16944, 21467, 2658, 5014, 4787, 5570, 5248, 4606, 6003, 5637, 4663, 5327, 13043, 13214, 7982, 8697, 18267, 25397, 25122, 19729, 21340, 23214, 34113, 11387, 6741, 8703, 6850, 3854, 1906, 818, 2256, 3817, 5912, 922, 1330, 653, 1190, 1439, 948, 1841, 1138, 1510, 1256, 1936, 3250, 2467, 3857, 4990, 7393, 7759, 8458, 7373, 5643, 9802, 5836, 3166, 2414, 2095, 1732, 982, 217, 483, 1138, 1416, 103, 580, 85, 213, 512, 405, 636, 1387, 1325, 1096, 1339, 1409, 867, 2006, 2206, 3364, 4465, 2399, 5553, 3696, 4180, 3102, 2158, 1677, 1642, 1586, 525, 138, 317, 722, 383 
    ], (ny, na)
  ) )
 
  # Abundance index from survey, numbers at age x year 
  A = transpose( reshape( [ 
    39252.75, 48310.15, 14746.33, 19850.89, 8138.48, 12688.36, 11826.61, 6203.89, 8257.22, 3893.25, 1569.87, 5046.63, 11037.2, 7779.2, 13105.01, 4428.91, 12428.39, 30706.3, 19493.08, 10909.37, 4311.93, 6102, 5337.9, 5631.81, 2923.19, 3005.84, 475, 989.61, 4952.48, 3699.85, 6487.94, 1922.5, 8265.87, 8295.23, 11169.09, 6952.84, 2911.54, 1623.12, 1977.26, 2119.12, 1866.32, 2115.07, 251.33, 549.04, 2164.04, 1781.99, 3598.3, 1473.02, 5981.27, 3932.19, 4236.12, 3159.43, 2379.53, 367.88, 619.04, 674.88, 526.44, 521.16, 143.15, 226.6, 463.47, 670.45, 1035.57, 493.51, 1999.34, 1729.14, 1276.54, 920.81, 713.98, 294.48, 271.89, 295.51, 239.32, 174.53, 72.62, 39.24, 175.02, 383.66, 331.54, 164.73, 783.08, 441.28, 965.76, 447.51, 240.69, 136.69, 141, 126, 76.02, 101, 26.61, 18.6, 45.01, 91.23, 131.94, 73.44, 254.82, 216.26, 390.27, 218.91, 56.18, 101, 154.6, 56.26, 17.24, 32.61, 4.63, 2.32, 34.06, 81.66, 52.97, 23.21 
    ], (A_ny, na) 
  ) )
  A = hcat(A)
  
# average body weight in kg for catch and stock by age
  body_weight_kg = transpose( reshape( [ 
    0.734, 0.287, 0.930, 0.659, 1.189, 1.150, 1.540, 1.881, 2.273, 2.669, 3.470, 3.418, 4.901, 4.630 
    ], (2, na)
  ))  # col1 = catch, col2=survey
  
# maturity ogive by age (fraction mature)
  maturity_ogive = [ 0.15, 0.44, 0.86, 0.96, 0.95, 0.98, 0.98 ]

  return ages, yrs, C, A_yrs, A, body_weight_kg, maturity_ogive

end

 
function data_atlantic_cod()

# Data (example) Myers and Cadigan (1995a) , analyzed by Millar and Meyer (2000)
  
ages = 1:108 # the last age is 8+
na = length(ages)

yrs = 1:16
ny = length(yrs)

A_yrs = 1981:1996  # S = survey
A_ny = length(A_yrs)

# catch in fishery in numbers at age x year
C = transpose( reshape( [ 
  1323, 1152, 2554, 2185, 1702, 2585, 
 782, 650, 831, 2329, 2779, 1696, 7693, 3111, 430, 940, 17556, 12361, 
 12025, 7172, 31286, 13616, 14871, 14824, 15219, 9217, 14651, 17639, 
 40557, 31654, 3860, 4993, 39206, 37439, 28814, 13191, 19003, 42602, 
 31760, 36614, 44168, 32340, 20184, 21150, 36410, 53805, 14535, 3343, 
 20319, 29202, 30016, 24800, 14297, 19028, 38624, 33922, 45869, 49061, 
 47917, 25212, 22695, 29553, 12211, 1940, 7711, 10982, 18017, 22014, 
 25435, 12044, 12503, 28006, 26025, 28469, 45725, 38708, 16390, 9064, 
 4526, 700, 3078, 3460, 4830, 11848, 16930, 14701, 7246, 7050, 14722, 
 19505, 18608, 28499, 17940, 6164, 1372, 147, 1530, 1300, 1217, 3175, 
 11936, 8934, 8910, 3836, 3104, 5818, 9026, 8696, 9156, 4745, 376, 21, 
 1083, 757, 520, 779, 1923, 6341, 4227, 5162, 2000, 1346, 4337, 3640, 
 2865, 1696, 199, 0.1, 437, 560, 232, 309, 338, 1018, 2536, 2905, 1977, 
 676, 774, 1695, 1084, 641, 104, 0.1, 324, 299, 285, 320, 246, 338, 597, 
 1935, 1675, 1264, 788, 816, 581, 338, 27, 0.1
  ], (ny, na)
) )

# Abundance index from survey, numbers at age x year 
A = transpose( reshape( [ 
 5.39, 1.94, 2.48, 5.12, 5.87, 12.22, 10.79, 7.27, 
 4.77, 2.04, 3.93, 8.98, 10.93, 3.35, 1.78, 0.6, 11.51, 11.78, 3.83, 
 2.74, 5.92, 10.62, 15.23, 12.35, 20.7, 4.03, 3.2, 8.3, 12.95, 13.97, 
 2.3, 0.83, 13.95, 16.79, 13.23, 3.26, 3.83, 10.83, 11.34, 10.01, 31.29, 
 13.23, 5.29, 6.2, 8.61, 9, 2.72, 0.34, 5.51, 10.53, 13.31, 9.67, 2.79, 
 3.87, 9.59, 7.28, 21.28, 11.61, 10.57, 6.52, 5.64, 3.31, 1.42, 0.22, 
 1.62, 2.27, 4.99, 8.79, 5.82, 2.43, 2.3, 4.24, 10.14, 4.38, 10.13, 8.23,
 3.9, 1.1, 0.35, 0.04, 0.63, 0.92, 1.19, 3.66, 5.31, 5.33, 1.37, 0.92, 
 5.26, 2.67, 2.58, 4.84, 3.98, 0.5, 0.04, 0.01, 0.47, 0.31, 0.37, 0.74, 
 2.59, 2.93, 2.09, 0.78, 1.37, 1.38, 1.55, 1.62, 1.68, 0.35, 0.02, 0.001, 
 0.33, 0.26, 0.23, 0.23, 0.57, 1.42, 1.3, 0.67, 0.58, 0.34, 0.79, 0.98, 
 0.55, 0.16, 0.01, 0.001, 0.12, 0.19, 0.11, 0.1, 0.16, 0.36, 0.54, 0.41, 
 0.68, 0.17, 0.15, 0.43, 0.23, 0.04, 0.001, 0.001, 0.09, 0.06, 0.16, 0.11,
 0.09, 0.14, 0.28, 0.15, 0.42, 0.19, 0.11, 0.16, 0.12, 0.02, 0.01, 0.001
  ], (A_ny, na) 
) )
A = hcat(A)


LogRecruits=c(13,12,12,13,13,13,13,13,12,12,12,12,12,11,11,11,11)
LogInitials=c(13,12,12,11,10,9,9,8,8,7)
PrecLogI= 10 

return ages, yrs, C, A, LogRecruits, LogInitial, PrecLogI

end



function setup_F( Fp, Ft; type="nonseparable", iplus=na )
  # starting assumption for terminal year and oldest age groups
  na, ny = size(C)
  F = tzeros( na, ny )   
  if type == "nonseparable"
    for i in iplus  # age indices of the plus age classes  
      F[ i, :] .= Fp  # by year in plus groups
    end
    F[ :, ny ]  = Ft   # by age in last/terminal year
  elseif type == "separable"  # ie.  multplicative
    F = Ft * Fp'
  end
  return F
end



function setup_N( C, F, M; iplus=na )
  # starting assumption for terminal year and oldest age groups
  # fill data with terminal years and ages
  na, ny = size(C)
  N = tzeros( na, ny )   # same aszeros but for Turing 
   
  Mdim = length(size(M))
  if Mdim==0

    for i in iplus  # age indices of the plus age classes  
      N[i,:] .=  C[i,:] ./ (F[i,:] ./ (F[i,:] .+ M) .* (1.0 .- exp.( .-(F[i,:] .+ M) )))
    end
    N[:, ny] = C[:,ny] ./ (F[:,ny] ./ (F[:,ny] .+ M) .* (1.0 .- exp.( -(F[:,ny] .+ M) )) )
 
  elseif Mdim==1 

    for i in iplus  # age indices of the plus age classes  
      N[i,:] .=  C[i,:] ./ (F[i,:] ./ (F[i,:] .+ M[i]) .* (1.0 .- exp.( .-(F[i,:] .+ M[i]) )))
    end
    N[:, ny] = C[:,ny] ./ (F[:,ny] ./ (F[:,ny] .+ M) .* (1.0 .- exp.( -(F[:,ny] .+ M) )) )
  
  elseif Mdim==2
 
    for i in iplus  # age indices of the plus age classes  
      N[i,:] .=  C[i,:] ./ (F[i,:] ./ (F[i,:] .+ M[i]) .* (1.0 .- exp.( .-(F[i,:] .+ M[i,:]) )))
    end
    N[:, ny] = C[:,ny] ./ (F[:,ny] ./ (F[:,ny] .+ M) .* (1.0 .- exp.( -(F[:,ny] .+ M[:,ny]) )) )
  
  end

  return N
end



function compute_fishing_mortality!(F, N, M; iplus=na)
  na, ny = size(N)
  
  Mdim = length(size(M))
  i_ages = reverse(setdiff( 1:na, iplus ) )

  ia = findall( x -> x < 0 || !isfinite(x), N )
  if length(ia) > 0
    N[ia] .= 1.0  # on log scale this becomes 0 .. continue but discourage
  end

  if Mdim==0
    for i in i_ages
      for j in 1:(ny-1)
        F[i,j] = log( N[i,j] /  N[i+1, j+1]) - M # update F estimates after convergence
      end
    end
  elseif Mdim==1 
    for i in i_ages
      for j in 1:(ny-1)
        F[i,j] = log( N[i,j] /  N[i+1, j+1]) - M[i] # update F estimates after convergence
      end
    end
  elseif Mdim==2
    for i in i_ages
      for j in 1:(ny-1)
        F[i,j] = log( N[i,j] /  N[i+1, j+1]) - M[i,j] # update F estimates after convergence
      end
    end
  end
  
  ib = findall( x -> x < 0 || !isfinite(x), F )
  if length(ib) > 0
    F[ib] .= 0.0  # on log scale this becomes 0 .. continue but discourage
  end

  return F

end


function summary_fishing_mortality(F; subset=1)

  Fall =  1.0 .- exp.(- mean( F, dims=1 ) ) 
  Fsubset = 1.0 .- exp.(- mean( F[subset,:], dims=1 ) )  # subset
  pl = plot( collect(yrs), vec(Fall), label="Mean" )
  pl = plot( pl, collect(yrs), vec(Fsubset), label="Subset" )
  return Fall, Fsubset, pl

end


function compute_biomass( B, maturity_ogive ) 
  Bmat  = B .* maturity_ogive
  Btot = sum( B, dims=1 )
  Bssb = sum( Bmat, dims=1 )  
  pl = plot( collect(yrs), vec(Btot), label="Total" )
  pl = plot!( pl, collect(yrs), vec(Bssb), label="SSB" ) 
  return Bmat, Btot, Bssb, pl
end


function cohort_analysis!( N, C, M; iplus=na ) 

  # back-propagate numbers of oldest groups based upon assumed natural and computed fishing mortalties 
  na, ny = size(C)
  Mdim = length(size(M))

  # iteratively fill in the matrix from older and later years
  # $N_{t} = N_{t+1} \; e^M + C_t \; e^{M/2}$     "Pope's" approximation
  i_ages = reverse(setdiff( 1:na, iplus ))

  if Mdim==0
    for i in i_ages
      for j in 1:(ny-1)
        N[i,j] = N[i+1, j+1] * exp(M) + C[i,j] * exp(M / 2)  
      end
    end
  elseif Mdim==1 
    for i in i_ages
      for j in 1:(ny-1)
        N[i,j] = N[i+1, j+1] * exp(M[i]) + C[i,j] * exp(M[i] / 2)  
      end
    end

  elseif Mdim==2
    for i in i_ages
      for j in 1:(ny-1)
        N[i,j] = N[i+1, j+1] * exp(M[i,j]) + C[i,j] * exp(M[i,j] / 2)  
      end
    end
  end
  
  return N
end



function catch_cohort_estimate( N, F, M ) 

  # back-propagate numbers of oldest groups based upon assumed natural and computed fishing mortalties 
  na, ny = size(C)
  Mdim = length(size(M))

  Cpred = tzeros(na,ny) 
  
  if Mdim==0
    for i in 1:na
      for j in 1:ny
        Cpred[i,j] = F[i,j] ./ ( F[i,j] .+ M ) .* N[i,j] .* (1.0 .- exp.(- (F[i,j] .- M) ) )
      end
    end
  elseif Mdim==1 
    for i in 1:na
      for j in 1:ny
        Cpred[i,j] = F[i,j] ./ ( F[i,j] .+ M[i] ) .* N[i,j] .* (1.0 .- exp.(- (F[i,j] .- M[i]) ) )
      end
    end

  elseif Mdim==2
    for i in 1:na
      for j in 1:ny
        Cpred[i,j] = F[i,j] ./ ( F[i,j] .+ M[i,j] ) .* N[i,j] .* (1.0 .- exp.(- (F[i,j] .- M[i,j]) ) )  
      end
    end
  end

  ia = findall( x -> x < 0 || !isfinite(x), Cpred )
  
  if length(ia) > 0
    Cpred[ia] .= 0.0  # on log scale this becomes 0 .. continue but discourage
  end

  return Cpred

end



function catchability_cohort_estimate( N, C ) 
 
  na, ny = size(N)
  Q = dQ = tzeros(na,ny)

  Q = N ./ C  # catchability
  ia = findall( x -> x >= 0 && isfinite(x), Q )
  ib = findall( x -> x < 0 || !isfinite(x), Q )
  
  if length(ib) > 0
    Q[ib] .= 1.0  # on log scale this becomes 0 .. continue but discourage
  end
  Q = log.(Q)
  
  q = tzeros(na)
  for k in 1:na
    j = findall( x -> x >= 0 && isfinite(x), Q[k,:] )
    q[k] = mean( Q[k,j] )  # across age
  end
  dQ = Q .- q
  sumdQ = sum( dQ[ia] )  
  return sumdQ, q, ia

end



Turing.@model function Virtual_Population_Analysis_basic( 
  Ftype="separable", M = repeat([0.2], na), Cobs=Cobs, iplus=na, ::Type{T} = Float64 ) where T

  Fp ~ arraydist( LogNormal.( log.(Fp0), 0.25) ) 
  Ft ~ arraydist( LogNormal.( log.(Ft0), 0.25) ) 
  Fs = setup_F( Fp, Ft; type=Ftype, iplus=iplus ) 
  Ns = setup_N( Cobs, Fs, M; iplus=iplus ) 
  Ns = cohort_analysis!( Ns, Cobs, M; iplus=iplus  ) 
  Fs = compute_fishing_mortality!(Fs, Ns, M; iplus=iplus )  # update F with finalized N's
  
  j = findall( x -> x <= 0, Fs )
  
  if length(j) > 0
    Fs[j] .= 0.001
    Turing.@addlogprob! -Inf
  end

  Cpred = catch_cohort_estimate( Ns, Fs, M ) 
 
  sumdQ, q, _ = catchability_cohort_estimate( Ns, Cpred ) 
  sumdQ ~ Normal(0.0, 0.001 * T(na*ny) );  # equivalent to mean(dQ) ~ normal(0, 0.001)

  # at present saving a result is not possible unless it is sampled
  S ~ arraydist( Normal.( Ns, 0.0001 ) )  # force a sample for all data to be kept
end



Turing.@model function Virtual_Population_Analysis( 
  Ftype="separable", M = repeat([0.2], na), Cobs=Cobs, iplus=na, ::Type{T} = Float64 ) where T
  Fp ~ arraydist( LogNormal.( log.(Fp0), 0.25) ) 
  Ft ~ arraydist( LogNormal.( log.(Ft0), 0.25) ) 
  Fs = setup_F( Fp, Ft; type=Ftype, iplus=iplus ) 
  Ns = T.(setup_N( Cobs, Fs, M; iplus=iplus ))  
  Ns = cohort_analysis!( Ns, Cobs, M; iplus=iplus  ) 
  Fs = compute_fishing_mortality!(Fs, Ns, M; iplus=iplus )  # update F with finalized N's
  
  j = findall( x -> x <= 0, Fs )
  
  if length(j) > 0
    Fs[j] .= 0.001
    Turing.@addlogprob! -Inf
  end

  Cpred = catch_cohort_estimate( Ns, Fs, M ) 
 
  sumdQ, q, ia = catchability_cohort_estimate( Ns, Cpred ) 
  sumdQ ~ Normal(0.0, 0.001 * T(na*ny) );  # equivalent to mean(dQ) ~ normal(0, 0.001)

  # at present saving a result is not possible unless it is sampled
  S ~ arraydist( Normal.( Ns, 0.0001 ) )  # force a sample for all data to be kept
  
  # the followng is an additional likelihood relative to the smple VPA form: 
  # due to negative values, simply using a "Poisson"  is more justified
  Cobs[ia] ~ arraydist( LogNormal.( log.(Cpred[ia]), 0.25) )  # # likelihood .. Broadcast not working?
end



Turing.@model function Adaptive_Virtual_Population_Analysis( Ftype="separable", M = repeat([0.2], na), iplus=na ) 
  Fp ~ arraydist( LogNormal.( log.(Fp0), 0.25) ) 
  Ft ~ arraydist( LogNormal.( log.(Ft0), 0.25) ) 
  F = setup_F( Fp, Ft; type=Ftype, iplus=iplus ) 
  N = setup_N( C, F, M; iplus=iplus ) 
  N = cohort_analysis!( N, C, M; iplus=iplus  ) 
  q ~ arraydist( LogNormal.( log.(q0), 0.25) ) 
  mu = log.(q .* N) 
  muyi = view( mu, :, yi ) # overlapping subset
  A ~ arraydist( LogNormal.( muyi, 0.25 ) )  # # likelihood .. Broadcast not working?
  S ~ arraydist( LogNormal.( mu , 0.25 ) )  # force a sample for all data to be kept
end



Turing.@model function Integrated_Catch_Analysis( Ftype="separable", M = repeat([0.2], na), C=C, iplus=na ) 
  Fp ~ arraydist( LogNormal.( log.(Fp0), 0.25) ) 
  Ft ~ arraydist( LogNormal.( log.(Ft0), 0.25) ) 
  F = setup_F( Fp, Ft; type=Ftype, iplus=iplus ) 
  N = setup_N( C, F, M; iplus=iplus ) 
  N = cohort_analysis!( N, C, M; iplus=iplus  ) 
  F[a, y] = log.( N[a, y] ./  N[a.+1, y.+1]) .- M[a] # update F estimates after convergence
  Cpred =  F ./ ( F .+ M ) .* N .* (1.0 .- exp.(- (F .- M) ) )
  C ~ arraydist( LogNormal.( Cpred, 0.25 ) )  # # likelihood .. Broadcast not working?
  q ~ arraydist( LogNormal.( log.(q0), 0.25) ) 
  mu = log.(q .* N) 
  muyi = view( mu, :, yi ) # overlapping subset
  A ~ arraydist( LogNormal.( muyi, 0.25 ) )  # # likelihood .. Broadcast not working?
  S ~ arraydist( LogNormal.( mu , 0.25 ) )  # force a sample for all data to be kept
end



Turing.@model function sequential_population_analysis()

# direct transcripion, untested

# Vague priors on first row and column 
  for j in 1:(ny+1)
    LogRecruits[j] ~ Normal(10,.01) 
    LogN[1,j] = LogRecruits[j] 
  end

  for i in 2:(na+1)
    LogInitials[i-1] ~ Normal(10,.01) 
    LogN[i,1] = LogInitials[i-1]
  end

  # Distribution for data, na[i,j]
  for i in 1:na 
  for j in 1:ny 
    MeanLogI[i,j] = Logq[i] + (LogN[i,j] + 11*LogN[i+1,j+1]) / 12;
    Prec[i,j] = PrecLogI / i;
    na[i,j] ~ LogNormal(MeanLogI[i,j], Prec[i,j])
  end
  end

  # Prior on mortalities, hierarchical # NB: log(0.2)=-1.61
  for i in 1:na
    MeanLogAgeM[i] ~ Normal(-1.61, 1)
  end

  for j in 1:ny
    MeanLogDeltaM[j] ~ Normal(0,1)
  end
  
  PrecLogM ~ Gamma(1,1)
  
  for i in 1:na
  for j in 1:ny
    MeanLogM[i,j] = MeanLogAgeM[i] + MeanLogDeltaM[j]
    M[i,j] ~ Truncated( LogNormal( MeanLogM[i,j], PrecLogM), 0, 10.0)
  end
  end
  
  # Cohort equation (starts with young in earliest time period (unlike traditional VPA)
  N = exp(LogN)
  N = cohort_analysis!( N, C, M; Mtype="forward" ) 
  logN = log.(maximum.(N, 10))
  
  # Prior on q's# 
  MeanLogq = -10; PrecLogq = 1;
  for i in 1:A
    Logq[i] ~ Normal(MeanLogq, PrecLogq)
  end

  # Prior on precision of log(na)# 
  PrecLogI~dgamma(10,1);

  # Vector of 1983 numbers for monitoring
  @show exp(LogN[:,6])/1000


#= original code copied from Millar and Meyer 2000 .. bugs version

    # Vague priors on first row and column#######
    for(y in 1:(Y+1)) {LogRecruits[y]~dnorm(10,.01);
      LogN[1,y]<-LogRecruits[y]; }
    for(a in 2:(A+1)){LogInitials[a-1]~dnorm(10,.01);
      LogN[a,1]<-LogInitials[a-1]; }

    ########Distribution for data, A[a,y]######
    for (a in 1:A) {
    for (y in 1:Y) { MeanLogI[a,y]<-Logq[a]+(LogN[a,y]+11*LogN[a+1,y+1])/12;
      Prec[a,y]<-PrecLogI/a;
      A[a,y]~dlnorm(MeanLogI[a,y],Prec[a,y]); } }

    #####Prior on mortalities, hierarchical######NB: log(0.2)=-1.61
    for(a in 1:A){ MeanLogAgeM[a]~dnorm(-1.61,1);}
    for(y in 1:Y){ MeanLogDeltaM[y]~dnorm(0,1);}
    PrecLogM~dgamma(1,1);
    for(a in 1:A) {
    for(y in 1:Y) { MeanLogM[a,y]<-MeanLogAgeM[a]+MeanLogDeltaM[y];
      M[a,y]~dlnorm(MeanLogM[a,y],PrecLogM) T(0,10.0); } }  ## in jags A -> T

    #####Cohort equation#####
    for(a in 2:(A+1)) {
    for(y in 2:(Y+1)) {
    LogN[a,y]<-log(max(exp(LogN[a-1,y-1])*exp(-M[a-1,y-1])-C[a-1,y-1]*exp(-0.5*M[a-1,y-1]),10)); } }


    ######Prior on q's######
    MeanLogq<--10; PrecLogq<-1;
    for (a in 1:A){ Logq[a]~dnorm(MeanLogq,PrecLogq); }

    ######Prior on q's, quadratic######
    #gamma~dnorm(-8,1)A(-10,-6); 
    #alpha~dunif(0.0001,0.02); 
    #beta~dnorm(5,0.25)A(1,11);
    #for (a in 1:A){ MeanLogq[a]<-gamma-alpha*(a-beta)*(a-beta);
    #                Logq[a]~dnorm(MeanLogq[a],16); }

    ######Prior on precision of log(A)######
    PrecLogI~dgamma(10,1);

    ######Vector of 1983 numbers for monitoring#####
    for(a in 1:A) {N1983s[a]<-exp(LogN[a,6])/1000; }
  =#
  
end


function separable_vpa_ices_variation()

# TODO .. not yet converted to Julia/Turing 
# mostly checking on priors and assumptions given lognormal structure  .. Poisson might be better

# Example of a seperable VPA implemented in WinBugs s

# Copied from an ICES CM 2004/D03. Report of the working group on methods of fish stock assessments. Copenhagen, Denmark
#
# Original by:Manica Azevedo: pp 116-127 (Appendix A)
# which was adapted from Nielsen A. 2000. Fish stock assessment using Markov Chain Monte Carlo. MSc thesis.
# Iberian Hake VPA

#=
# likelihoods

# Catch data (ages: 0-8, years: 1982-2002, unit: thousands)
  for( a in 1:A ) {
    for( y in 4:Y ) {
      logC[a,y] <- log( C[a,y] )
      logC[a,y] ~ dnorm( logmu.C[a,y], tau.C )
      logmu.C[a,y] <- log( F[a,y] / Z[a,y] * N[a,y] * ( 1 - exp(-Z[a,y])) ) )
  }}

# Survey data : Spanish and Portugese
# Survey data October (ages: 0-8, years: 1985-2002, unit: thousands)
  for( a in 1:A ) {
    for( y in 4:Y ) {
      logSO[a,y] <- log( SO[a,y-3] )
      logSO[a,y] ~ dnorm( logmu.SO[a,y], tau.s )
      logmu.SO[a,y] <- log( qs1[a,y] * exp(-Z[a,y]*10/12) * N[a,y] ) 
  }}

# Survey data July (ages: 0-8, years: 1989-1993, 1995, 1997-2001, unit: thousands)
  for( a in 1:A ) {
    for( y in 8:12 ) {
      logSJ[a,y] <- log( SJ[a,y-7] )
      logSJ[a,y] ~ dnorm( logmu.SJ[a,y], tau.s )
      logmu.SO[a,y] <- log( qs1[a,y] * exp(-Z[a,y]*7/12) * N[a,y] ) 
  }}

  for( a in 1:A ) {
    y <- 14
      logSJ[a,y] <- log( SJ[a,y-7] )
      logSJ[a,y] ~ dnorm( logmu.SJ[a,y], tau.s )
      logmu.SO[a,y] <- log( qs1[a,y] * exp(-Z[a,y]*7/12) * N[a,y] ) 
  }

# etc ..

# Separable F (age, year)
  for( a in 1:A ) {
    for( y in 1:Y ) {
      Fay[a,y] <- Fa[a] * Fy[y]
  }}

# Natural and total mortality
  for( a in 1:A ) {
    for( y in 1:Y ) {
      Z[a,y] <- F[a,y] + M[a,y]
      F[a,y] <- Fay[a,y]
  }}

  
# Catchability
  for( a in 1:A) {
    for(y in 4:Y) {
      qs1[a,y] <- qs11.sa[a]
  }}

  for( a in 1:A) {
    qs2[a,14] <- qs22.sa[a]
    for( y in 8:12) {
      qs2[a,y] <- qs22.sa[a]
    }
    for( y in 16:20) {
      qs2[a,y] <- qs22.sa[a]
    }
  }}
# etc .,,,

# N-at-age
  for(y in 1:Y) {
    N.year1[y] ~ dunif(0,1)
    N[1,y] <- N.year1[y] * 1.0E6
  }
  N.age[1] <- N.year[1]

  for(a in 1:(A-1)) {
    N.age1[a] ~ dunif(0,1)
    N.age[a+] <- N.age1[a]
    N[a+,1] <- N.age[a+] * 2.0E5
  }

  for(a in 2:A) {
    for(y in 2:(Y+1)) {
      N[a,y] <- N[a-1,y-1] * exp(-Z[a-1, y-1])
  }}


# Priors
  # Fa
    for (a in1:A) {
      Fa1[a] ~ dunif(0,2)  
      Fa[a] <- Fa1[a]
    }

  # Fy
    for(y in 1:Y) {
      Fy1[y] ~ dunif(0,2)
      Fy[y] <- Fy1[y]
    }

  # M
    for(a in 1:A) {
      for(y in 1:Y) {
        M1[a,y] ~ dunif(0.1, 0.3)
        M[a,y] <- M1[a,y]
      }
    }

  # catchability
    for(a in1:A) {
      qs11.sa[a] ~ dunif(6.0E-6,1)
      qs22.sa[a] ~ dunif(6.0E-6,1)
    }

  # tau
    tau.C ~ dgamma(4, 0.4)
    var.C <- 1/tau.C

    tau.s ~ dgamma(100, 100)
    var.s <- 1/tau.s

# Stock characteristics
  # R
    for( y in 1:Y) {
      R[y] <- N[1,y]
    }

  # B, SSB
    for(a in 1:A) {
      for(y in 1:Y) {
        B[a,y] <- N[a,y] * w[a,y]
        SSB[a,y] <- B[a,y] * mat[a,y]
    }}

  # Fbar (age range: mina to maxa)
    for( y in 1:Y) {
      Fsum[y] <- sum(F[mina:maxa,y])
      Fbar[y] <- (Fsum[y])/(maxa-mina+1)
    }

  # Catch and catachability
  # residuals
    for(a in 1:A) {
      for(y in 1:Y) {
        C.est[a,y] <- F[a,y] / Z[a,y] * N[a,y] * (1-exp( -Z[a,y]))
        C.est[a,y] <- C[a,y] - C.est[a,y]
    }}
    C.resT <- sum(C.res[1:A, 1:Y])

  # Catchability
    for(a in 1:A) {
      for(y in 1:Y) {
        q.est[a,y] <- Survey[a,y] / N[a,y]
    }}
=#

end


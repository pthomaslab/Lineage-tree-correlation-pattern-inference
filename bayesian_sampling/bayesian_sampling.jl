## bayesian_sampling.jl
# Please read README.md for full details

## This script does the following:
# 1 - function set up
# 2 - MCMC Simulation with Adaptive Metopolis-within-Gibbs Sampling
# 3 - evaluate relevant output and export

## REQUIREMENTS - see README.txt for details
# Requires a .txt file called [celltype]_ccd_data_corrs.txt
#   which has correlations & standard deviations in interdivision time
#   for mother-daughter, grandmother-granddaughter, sister-sister
#   and cousin-cousin correlations. See README.txt for example format.
# If using the MATLAB pipeline to calculate correlations, this file is
#   automatically created in the correct format.

## files exported - see README.txt for further details on format.
# ..._sampling_params.dat - parameter sampling file
# ..._output.txt - full sampling output file, has all relevant outputs
# ..._output_thin.txt - thinned output file
# ..._inf.txt - produces inferred correlation coefficients and 95% CI
# ..._pattern_dist.txt - % of samples that fall into each pattern
# ..._mlp.txt - maximum likelihood parameter set
# ..._mle.txt - correlations calculated from maximum likelihood parameter set

## LOAD PACKAGES
using MatrixEquations
using Distributions
using Plots
using DelimitedFiles
using LinearAlgebra
using StatsBase
using GaussianMixtures
using Mamba

## CHOOSE INPUTS
println("Part 1: Set up model...")

println("Cell type:")
celltype = readline()
println("Number of samples desired:")
act = parse(Float64,readline())
println("Number of sampled wanted after thinning")
des = parse(Float64,readline())
println("What did you pair (eg IDT, etc)")
choice = readline()

## SET CONDITIONS

# set initial condition
th = [0. 0. 0. 0.]
la = [1. 1.]
ga = [0.]
de = [0. 0. 0.]

pp = [ th la ga de ]

# alphas both = 1, cell cycle phase model
alpha1 = 1;
alpha2 = 1;

# steps for sampling
burnin = Int(act/10) # burnin set to 10%
n = Int(act + burnin) # samples including burnin

# assign name for file labeling
name = string(celltype, "_", (n-burnin))

# obtain file with correlation and distributions
ccd = readdlm( string(pwd(),"\\",celltype,"_",choice,"_corrs.txt"),'\t', Float64, '\n')

## FUNCTION SET UP
# parameters required
# 1 - theta11
# 2 - theta12
# 3 - theta21
# 4 - theta22
# 5 - lambda1
# 6 - lambda2
# 7 - gamma12
# 8 - delta11
# 9 - delta12
# 10 - delta22

# define all functions
function theta(paras)
    return [ [paras[1] paras[2]] ; [paras[3] paras[4]] ]
end

function alphas()
    return [ [alpha1] ; [alpha2] ]*transpose([ [alpha1] ; [alpha2] ])
end

function varmat(paras)
    return [ [paras[5]] ; [paras[6]] ]*transpose([ [paras[5]] ; [paras[6]] ])
end

function gammat(paras)
    return [[ 1 paras[7] ] ;  [ paras[7] 1 ]]
end

function delmat(paras)
    return [[paras[8] paras[9]];
            [paras[9] paras[10]] ]
end

function S1(paras)
    return varmat(paras).*gammat(paras)
end

function S2(paras)
    return varmat(paras).*delmat(paras)
end

function Lambda(paras)
    return [[S1(paras) S2(paras)] ; [S2(paras) S1(paras)]]
end

# define covariance matrices
function omega0(paras)
    return lyapd(theta(paras),S1(paras))
end

function omegaMD(paras)
    return omega0(paras)*transpose(theta(paras))
end

function omegaSS(paras)
    return theta(paras)*omega0(paras)*transpose(theta(paras)) + S2(paras)
end

function omegaCC(paras)
    return theta(paras)*omegaSS(paras)*transpose(theta(paras))
end

function omegaGG(paras)
    return omegaMD(paras)*transpose(theta(paras))
end

function omegaGR(paras)
    return omegaGG(paras)*transpose(theta(paras))
end

function omegaGR2(paras)
    return omegaGR(paras)*transpose(theta(paras))
end

# define covariances
function cov0(paras)
    return sum(omega0(paras).*alphas())
end

function covMD(paras)
    return sum(omegaMD(paras).*alphas())
end

function covSS(paras)
    return sum(omegaSS(paras).*alphas())
end

function covCC(paras)
    return sum(omegaCC(paras).*alphas())
end

function covGG(paras)
    return sum(omegaGG(paras).*alphas())
end

function covGR(paras)
    return sum(omegaGR(paras).*alphas())
end

function covGR2(paras)
    return sum(omegaGR2(paras).*alphas())
end

# define correlations
function omh0(paras)
    return cov0(paras)
end

function corhMD(paras)
    return covMD(paras) / cov0(paras)
end

function corhSS(paras)
    return covSS(paras) / cov0(paras)
end

function corhCC(paras)
    return covCC(paras) / cov0(paras)
end

function corhGG(paras)
    return covGG(paras) / cov0(paras)
end

function corhGR(paras)
    return covGR(paras) / cov0(paras)
end

function corhGR2(paras)
    return covGR2(paras) / cov0(paras)
end

# phase correlation matrices
function cor0(paras)
    D = sqrt(omega0(paras).*Matrix(I,2,2))
    Dinv = inv(D)
    return Dinv*omega0(paras)*Dinv
end

function bigmat(paras)
    b1 = [omega0(paras) omegaMD(paras)]
    b2 = [transpose(omegaMD(paras)) omega0(paras)]
    bigmat = [b1 ; b2]
    return bigmat
end

function MD0(paras)
    D = sqrt(bigmat(paras).*Matrix(I,4,4))
    Dinv = inv(D)
    return Dinv*bigmat(paras)*Dinv
end

## LIKELIHOOD FUNCTION
function like(paras,ccd)

    mu = ccd[1,1]
    om0 = ccd[2,1]

    corMD = ccd[3,1]
    corSS = ccd[4,1]
    corCC = ccd[5,1]
    corGG = ccd[6,1]

    sigom0 = ccd[2,4]

    sigMD = ccd[3,4]
    sigSS = ccd[4,4]
    sigCC = ccd[5,4]
    sigGG = ccd[6,4]

    val =
    ((corMD - corhMD(paras))^2)/(sigMD^2) +
    ((corSS - corhSS(paras))^2)/(sigSS^2) +
    ((corCC - corhCC(paras))^2)/(sigCC^2) +
    ((corGG - corhGG(paras))^2)/(sigGG^2) +
    ((om0 - omh0(paras))^2)/(sigom0^2)

## CONDITIONS ON PARAMETERS, ALWAYS ASSUMING S2 (DELTA MATRIX) = ZEROS
## STANDARD CONDITIONS
    if eltype(theta(paras))<:Real && # THETA VALUES REAL
        all(abs.(eigvals(theta(paras))) .< 1) && # LYAPUNOV STABILITY OF THETA
        ishermitian(Lambda(paras)) && # ENSURE Lambda POSITIVE SEMI-DEFINITE
        all(eigvals(Lambda(paras)) .>= 0 ) && # ENSURE Lambda POSITIVE SEMI-DEFINITE
        all(abs.(gammat(paras)) .<= 1) && # gamma12 BETWEEN -1 AND 1
        all(varmat(paras) .> 0) && # LAMBDA VARIANCES POSITIVE
        all(abs.(delmat(paras)) .<= 1) && # DELTAS BETWEEN -1 AND 1
        all( [ paras[8] paras[9] paras[10] ] .== 0 ) # DELTAS 8-10 FIXED TO ZERO
        ## ABOVE CONDITION SETS S2 = ZEROS

        return val/2

    else

        return Inf

    end

end

# evaluate likelihood function at ccd
logf = function(p)
    return -like(p,ccd)
end

println("Model set up complete")

## MCMC Simulation with Adaptive Metopolis-within-Gibbs Sampling
println("Part 2: Begin sampling")

# get initial condition and evaluate likelihood at IC
pInit = pp
lh1 = like(pInit, ccd)

# file name
fname = string(name,"_sampling.params.dat")

# create file with initial condition (if file not exists)
if !isfile(fname)
  outfile = open(fname,"w")
  for s in pInit
     write(outfile,"$s\t")
  end
  write(outfile,"$lh1\n")
  close(outfile)
end

# read last parameter set from existing (if it exists)
outfile = open(fname)
global lastline=""
global n0=0
for l in eachline(outfile)
  global lastline=l
  global n0+=1
end

pInit=parse.(Float64,split(lastline))[1:end-1]
close(outfile)

println("init:","$pInit")

# start sampling using Mamba library
fits = Mamba.AMWGVariate(pInit, 1.0, logf)
for i in n0:n
  if i>n0
  outfile = open(fname,"a")
  for s in fits[1:length(pInit)]
     write(outfile,"$s\t")
  end

  lh1 = logf(fits)
  write(outfile,"$lh1")
  write(outfile,"\n")
  close(outfile)
  end

  Mamba.sample!(fits, adapt = (i <= burnin)) #sampling
end

println("Sampling complete.")

## CREATE RELEVANT OUTPUT FROM SAMPLING
println("Part 3: Create output file")

# convert .dat file into an array
sparas = readdlm(string(name,"_sampling.params.dat"), '\t', Float64, '\n')
# REMOVE BURN IN
sparas = sparas[burnin+1:end,:]

## FOCUS ON THETA: EIGENVALUES, PERIOD, DAMPING, pattern CLASSIFICATION
prob = zeros(ComplexF64,size(sparas)[1],9)

for i = 1:size(sparas)[1]
    prob[i,1:2] = eigvals(theta(sparas[i,:])) # eigenvalues
    prob[i,3] = 2*pi/angle(prob[i,1]) # period 1
    prob[i,4] = 2*pi/angle(prob[i,2]) # period 2
    prob[i,5] = norm(prob[i,2]) # damping
    if all( [imag(prob[i,1]) imag(prob[i,2])] .==0) # if eigenvalues real
        if all([real(prob[i,1]) real(prob[i,2])].>0) # if both positive
            prob[i,8] = 1 # then aperiodic
        else
            prob[i,7] = 1 # otherwise, alternator
        end
    else
        prob[i,6] = 1 # if complex, oscillator
    end
end

patterns = real(prob[:,6:8])

pattern_dist = (sum(patterns,dims=1)./size(patterns,1)).*100

open(string(name,"_pattern_dist.txt"),"w") do io
    writedlm(io,pattern_dist,' ')
end

# find all alternator
r2 = findall(patterns[:,2].==1)
period = real(prob[:,4])

# for +ve -ve eig val, replace period with 2 if inf period taken
if sum(patterns[:,2]) .> 0
period[r2,:] .= 2
end

damping = real(prob[:,5])

## 'Underlying periods' - shift frequency to obtain alternative periods

perds = zeros(size(sparas)[1],10)
meanccd = ccd[1,1]

for i = 1:size(sparas)[1]
    perds[i,1] = (2*pi/(angle(prob[i,1])))*meanccd
    perds[i,2] = (2*pi/(angle(prob[i,1])+2*pi))*meanccd
    perds[i,3] = (2*pi/(angle(prob[i,1])-2*pi))*meanccd
    perds[i,4] = (2*pi/(angle(prob[i,1])+4*pi))*meanccd
    perds[i,5] = (2*pi/(angle(prob[i,1])-4*pi))*meanccd

    perds[i,6] = (2*pi/(angle(prob[i,2])))*meanccd
    perds[i,7] = (2*pi/(angle(prob[i,2])+2*pi))*meanccd
    perds[i,8] = (2*pi/(angle(prob[i,2])-2*pi))*meanccd
    perds[i,9] = (2*pi/(angle(prob[i,2])+4*pi))*meanccd
    perds[i,10] = (2*pi/(angle(prob[i,2])-4*pi))*meanccd
end

## OUTPUT CORRELATIONS

function corrvec(paras)
    cv = [corhMD(paras) corhGG(paras) corhGR(paras) corhGR2(paras) corhSS(paras) corhCC(paras)]
return cv
end

corrs = zeros(size(sparas)[1],6)
for i = 1:size(sparas)[1]
    corrs[i,:] = corrvec(sparas[i,:])
end

## CREATE MATRIX FOR PLOTTING
# [ median , LQ - median , UQ - median]
inf = zeros(6,3)
for i = 1:6
    inf[i,:] =[ median(corrs[:,i]) quantile(corrs[:,i],0.025)-median(corrs[:,i]) quantile(corrs[:,i],0.975) - median(corrs[:,i]) ]
end
# OUTPUT FILE FOR PLOTTING FIT IN MATHEMATICA
open(string(name,"_inf.txt"),"w") do io
    writedlm(io, inf, ' ')
end

## EIGENVALUES AND PHASE CORRELATIONS
function eigv(ths)
    t = [[ ths[1] ths[2] ] ; [ ths[3] ths[4] ]]
    out = eigvals(t)
    line = [ real(out[1]) imag(out[1]) real(out[2]) imag(out[2])]
    return line
end

eigandcorr = zeros(size(sparas,1),9)

for i = 1:size(sparas,1)
    eigandcorr[i,:] = [eigv(sparas[i,1:4]) cor0(sparas[i,1:10])[1,2] Transpose(MD0(sparas[i,1:10])[3:4,1:2][:])]
end

## OUTPUT
output = [sparas corrs period damping patterns perds eigandcorr]

## THIN OUTPUT
tloc = [1:act./des:act;]
tloc = floor.(Int,tloc)
output_thin = output[tloc,:]

## Export .txt files

mxl = maximum(sparas[2:end,end])
loc = findall(sparas[:,end].==mxl)
mlp = sparas[loc[1],:]

open(string(name,"_mlp.txt"),"w") do io
    writedlm(io,mlp,' ')
end

mlc = corrvec(mlp)
open(string(name,"_mle.txt"),"w") do io
    writedlm(io,mlc,' ')
end

open(string(name,"_output.txt"),"w") do io
    writedlm(io,output,' ')
end

open(string(name,"_output_thin.txt"),"w") do io
    writedlm(io,output_thin,' ')
end

println("Finish")

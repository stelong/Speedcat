## Insert absolute path of "module"
push!(LOAD_PATH, "/home/sl18/julia_scripts/julia_master/Speedcat/module/")

## Import modules
using DelimitedFiles
using Plots
using Revise
using Solver

## Initialize simulation options
rat = "ab" # possible rat phenotype choices are "sham" and "ab"
freq = 6 # possible Hz choices are 1 and 6
n_beats = 1000 # number of cardiac cycles to simulate
n_curves = 1 # number of limit cardiac cycles to observe
soltype = "full" # if set to "fragmentary", the solution will be saved only at stim_period-spaced time points

## Run Gattoni EP model unperturbed
# t, sol = solveEp(rat, freq, n_beats, n_curves, soltype)
	
## Run Gattoni EP model with pertubations given by the k-th row of matrix S
# S = readdlm("data.txt", Float64) # read perturbation matrix (scaling prefactors for conductances)
# k = 10 # choose which row
# perturb = S[k, :] # actual perturbation vector
# t, sol = solveEp(rat, freq, n_beats, n_curves, soltype, perturb)

## Perturbed components are: [g_Na, g_to, g_ss, g_K1, g_f, i_NaK_max, N, g_NCX, g_SERCA, g_pCa, g_SRl]
## NOTE: all ones = unperturbed
perturb = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
t, sol = solveEp(rat, freq, n_beats, n_curves, soltype, perturb)

## Plot action potential and intracellular calcium concentration over time
p1 = plot(t, sol[13, :], label="AP")
p2 = plot(t, 1e3*sol[9, :], label="Ca2+")
plot(p1, p2, layout=(2, 1))

## Save time and calcium transient vectors
# M = [t 1e3*sol[9, :]] # two-column matrix
# mypath = "/home/sl18/Desktop/" # where to save
# filename = rat*"_"*string(freq)*"_"*"Ca_transient.txt"
# open(mypath*filename, "w") do f
# 	writedlm(f, M)
# end
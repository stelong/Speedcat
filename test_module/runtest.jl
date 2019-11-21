# Insert absolute path of "Ep" module
push!(LOAD_PATH, "/home/sl18/julia_scripts/julia_master/Speedcat/module/")

# Import modules
using DelimitedFiles
using Ep
using Plots
using Revise
using Statistics

# Read perturbation matrix (scaling prefactors for conductances)
S = readdlm("data.txt", Float64) 

# Initialize simulation options
rat = "ab" # possible rat phenotypes are "sham" and "ab"
freq = 6 # possible Hz choices are 1 and 6
n_beats = 1000 # number of simulated cycles
n_curves = 10 # number of limit cycles we want to observe
soltype = "full" # if set to "fragmentary", the solution will be saved only at stim_period-spaced time points

# Run Gattoni model (unperturbed)
# t, sol = solve_ep(rat, freq, n_beats, n_curves, soltype)
	
# or run Gattoni model with kth row perturbation
# k = 14 # choose which row
# t, sol = solve_ep(rat, freq, n_beats, n_curves, soltype, S[k, :])
t, sol = solve_ep(rat, freq, n_beats, n_curves, soltype, [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]) # ---> still unperturbed
#														 [g_Na, g_to, g_ss, g_K1, g_f, i_NaK_max, N, g_NCX, g_SERCA, g_pCa, g_SRl]

# Plot action potential and intracellular calcium concentration over time
p1 = plot(t, sol[13, :], label="AP")
p2 = plot(t, 1e3*sol[9, :], label="Ca2+")
plot(p1, p2, layout=(2, 1))
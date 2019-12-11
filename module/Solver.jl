module Solver
export solveEp

using Initializer
using Model
using Revise
using Sundials

function solveEp(rat, freq, n_beats, n_curves, soltype, args...)
	p = initConstants(rat, freq)
	stim_period = p[5]

	if length(args) > 0
		S = args[1]
		l = [11, 13, 17, 18, 19, 23, 34, 54, 55, 57, 60]
		p[l] = S.*p[l]
	end
	
	tspan = (0.0, n_beats*stim_period-1.0)

	if soltype == "fragmentary"
		dt = stim_period
	else
		dt = 1.0
	end

	t = collect((n_beats-n_curves)*stim_period:dt:n_beats*stim_period-1.0)

	u0 = initStates(rat, freq)
	prob = ODEProblem(computeRates, u0, tspan, p)
	
	alg = CVODE_BDF()
	sol = solve(prob, alg, dtmax=1.0, saveat=t)

	return t.-(n_beats-n_curves)*stim_period, sol
end

end
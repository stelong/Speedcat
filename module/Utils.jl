module Utils
export checkSteadystate 

using Statistics

function checkSteadystate(sol)
	Y = hcat(sol[9, :], sol[11, :], sol[13, :])
	vs = abs.((reshape(mapslices(maximum, Y, dims=1), 3)-reshape(mapslices(mean, Y, dims=1), 3))./reshape(mapslices(mean, Y, dims=1), 3))
	vi = abs.((reshape(mapslices(minimum, Y, dims=1), 3)-reshape(mapslices(mean, Y, dims=1), 3))./reshape(mapslices(mean, Y, dims=1), 3))
	vm = mean(0.5*(vi+vs))
	vm < 1e-1
end

end
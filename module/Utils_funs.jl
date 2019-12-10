# Define module
module Utils_funs

export check_steadystate # make this function visible outside

# Import module
using Statistics

# Check if the solution has reached some sort of steady-state
function check_steadystate(sol)
	Y = hcat(sol[9, :], sol[11, :], sol[13, :])
	vs = abs.((reshape(mapslices(maximum, Y, dims=1), 3)-reshape(mapslices(mean, Y, dims=1), 3))./reshape(mapslices(mean, Y, dims=1), 3))
	vi = abs.((reshape(mapslices(minimum, Y, dims=1), 3)-reshape(mapslices(mean, Y, dims=1), 3))./reshape(mapslices(mean, Y, dims=1), 3))
	vm = mean(0.5*(vi+vs))
	vm < 1e-1
end

end
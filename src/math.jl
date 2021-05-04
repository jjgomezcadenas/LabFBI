using GLM
"""
	lfit(x::Array{Float64}, y::Array{Float64})

Linear fit (syntactic sugar over lm function)

#Argument
- `x::Array{Float64}` :: Array represeting x data
- `y::Array{Float64}` :: Array represeting y data
"""
function lfit(x::Array{Float64}, y::Array{Float64})
	fit = lm(hcat(ones(length(x)), x), y)
	return fit
end

"""
	sline(x::Number, x1::Number, x2::Number) = x1 + x2 * x

A straight line
"""
sline(x::Number, x1::Number, x2::Number) = x1 + x2 * x

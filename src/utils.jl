using Printf

"""
	float_to_fstr(val, fmt) (obsolete)

Convert val into a formatted string using @sprintf
In spite of the function name (kept to distract the enemy)
the function will convert any value that can be formatted with @sprintf

# Arguments
- `val::Any`    : value to be formatted
- `fmt::String` : A format
"""
function float_to_fstr(val, fmt)
	ex = quote
	@sprintf $fmt $val
	end
	eval(ex)
end

"""
	to_fstr(val::Any, fmt::String)

Convert val into a formatted string using @sprintf

# Arguments
- `val::Any`    : value to be formatted
- `fmt::String` : A format
"""
function to_fstr(val::Any, fmt::String)
	ex = quote
	@sprintf $fmt $val
	end
	eval(ex)
end

"""
	vect_to_fstr(vect::AbstractVector, fmt::String)

Converts a vector into a formatted string

# Arguments
- `vect::AbstractVector` : vector to be formatted
- `fmt::String`          : A format
"""
function vect_to_fstr(vect::AbstractVector, fmt::String)
	vs = to_fstr.(vect,(fmt,))
	str = ""
	for v in vs[1:end-1]
		str = str * v * ", "
	end
	str = str * vs[end]
	str
end

"""
	vect_to_list(vect::AbstractVector, fmt::String)

Converts a vector into a list
To display in Pluto use: Text(list)

# Arguments
- `vect::AbstractVector` : vector to be formatted
- `fmt::String`          : A format
"""
function vect_to_list(vect::AbstractVector, fmt::String)
	vs = ut.to_fstr.(vect,(fmt,))
	str = ""
	for v in vs[1:end-1]
		str = string(str," - ", v,"\n")
	end
	str = string(str," - ", vs[end])
	str
end


"""
	logrange(x1::Number, x2::Number, n::Number)

returns a logarithmic range

# Arguments
- `x1::Number`     : start of the logrange
- `x2::Number`     : end of the logrange
- `length::Number` : length of the range
"""
function logrange(x1::Number, x2::Number, n::Number)
	return (10^y for y in range(log10(x1), log10(x2), length=n))
end


"""
	unzip(x) = collect(zip(x...))

Syntactic sugar for unzipping a collection

# Arguments
- `x1::any`     : name of collection to unzip
"""
unzip(x) = collect(zip(x...))

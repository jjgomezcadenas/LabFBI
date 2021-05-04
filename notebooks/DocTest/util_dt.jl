### A Pluto.jl notebook ###
# v0.14.3

using Markdown
using InteractiveUtils

# ╔═╡ 9ffc1790-a984-11eb-3277-f3725d671eee
begin
	using CSV
	using DataFrames
	using PlutoUI
	using Shapefile
	using ZipFile
	using LsqFit
	using Plots
	using LaTeXStrings
	using Statistics
	using Dates
	using Unitful
	using UnitfulEquivalences
	using Interpolations
	using QuadGK
	using Test
end

# ╔═╡ a317a082-ea26-4a3c-b2fb-666b84a3a143
using DrWatson

# ╔═╡ 1a1f50f5-27dd-4ea6-9b0e-66cc97a21d73
using Printf

# ╔═╡ d997f4b7-4ad0-4686-9fed-ed166bcf9041
md"# Test and documentation of functions in module `utils.jl`"

# ╔═╡ 6d7c7865-a8d4-4a75-b5ef-d09bbd3c1b5a
@quickactivate "LabFBI"

# ╔═╡ bcc21af7-bfd3-4adf-accf-26edeeacfbfd
projectdir()

# ╔═╡ 8c0d2002-b375-4489-b136-2f80debe0c33
function ingredients(path::String)
	# this is from the Julia source code (evalfile in base/loading.jl)
	# but with the modification that it returns the module instead of the last object
	name = Symbol(basename(path))
	m = Module(name)
	Core.eval(m,
        Expr(:toplevel,
             :(eval(x) = $(Expr(:core, :eval))($name, x)),
             :(include(x) = $(Expr(:top, :include))($name, x)),
             :(include(mapexpr::Function, x) = $(Expr(:top, :include))(mapexpr, $name, x)),
             :(include($path))))
	m
end

# ╔═╡ 67839f51-80c0-4924-8939-e6c194d33f96
ut = ingredients(srcdir("utils.jl"))

# ╔═╡ 62f419a9-fce0-40e3-a300-8faa4005365b
md"## Notebook"

# ╔═╡ 244b025e-ddd8-4969-9a0d-25a12677d619
md"""
	to_fstr(val::Any, fmt::String)

Convert val into a formatted string using @sprintf

\# Arguments
- `val::Any`    : value to be formatted
- `fmt::String` : A format
"""

# ╔═╡ 489348ef-7800-4518-aff3-4c003c9c816a
@sprintf "%7.2f" π

# ╔═╡ 93314cb4-5d95-4444-b4f2-f642721b0c81
@test ut.to_fstr(π, "%7.2f")== @sprintf "%7.2f" π

# ╔═╡ 824e0e46-aa82-48fb-af76-9f2ff5131f70
@test ut.to_fstr(10^π, "%7.3g") =="1.39e+03"

# ╔═╡ 8ee4e9b1-669a-4cb9-96ca-f752224abfba
@test ut.to_fstr(1, "%d") =="1"

# ╔═╡ b4b95535-62b5-4cfc-bff2-9a85d47e2ff9
@test ut.to_fstr("cat", "%s") =="cat"

# ╔═╡ 7daa3461-30ee-4982-9dad-3978b2928200
md"""
	vect_to_fstr(vect::AbstractVector, fmt::String)

Converts a vector into a formatted string

\# Arguments
- `vect::AbstractVector` : vector to be formatted
- `fmt::String`          : A format
"""

# ╔═╡ a4eaf88b-3582-4d5c-8bed-81bdde06d03e
@test ut.vect_to_fstr([1,2,3,4], "%d") == "1, 2, 3, 4"

# ╔═╡ d13d7550-7723-4a3c-a4c3-616b9bbf6f52
@test ut.vect_to_fstr([π,log(π),π^2,sqrt(π)], "%7.4f") ==" 3.1416,  1.1447,  9.8696,  1.7725"

# ╔═╡ 0e64f4ea-33dd-4345-b896-c02b3c4bb9c3
ut.vect_to_fstr([π,log(π),π^2,sqrt(π)], "%7.4f")

# ╔═╡ 5a530357-ccfc-42d4-aaeb-e0ce6a7ba06f
@test ut.vect_to_fstr(["the", "sun", "shines"], "%s") == "the, sun, shines"

# ╔═╡ fc59f914-70b3-4800-b14f-cc72a3b594c5
md"""
	logrange(x1::Number, x2::Number, n::Number)

returns a logarithmic range

\# Arguments
- `x1::Number`     : start of the logrange
- `x2::Number`     : end of the logrange
- `length::Number` : length of the range
"""

# ╔═╡ 5607c259-2c5e-4e8d-a53c-d3250829dc89
c1 = collect(ut.logrange(1,10,10)) 

# ╔═╡ 4014626f-0ae2-47fc-8af0-6d3cc92232fa
c2 = collect((10^y for y in range(log10(1), log10(10), length=10)))

# ╔═╡ 239efc3c-df42-4202-b71b-ae6dd53c4e83
@test c1 ≈ c2

# ╔═╡ Cell order:
# ╠═9ffc1790-a984-11eb-3277-f3725d671eee
# ╠═a317a082-ea26-4a3c-b2fb-666b84a3a143
# ╠═1a1f50f5-27dd-4ea6-9b0e-66cc97a21d73
# ╠═d997f4b7-4ad0-4686-9fed-ed166bcf9041
# ╠═6d7c7865-a8d4-4a75-b5ef-d09bbd3c1b5a
# ╠═bcc21af7-bfd3-4adf-accf-26edeeacfbfd
# ╟─8c0d2002-b375-4489-b136-2f80debe0c33
# ╠═67839f51-80c0-4924-8939-e6c194d33f96
# ╟─62f419a9-fce0-40e3-a300-8faa4005365b
# ╠═244b025e-ddd8-4969-9a0d-25a12677d619
# ╠═489348ef-7800-4518-aff3-4c003c9c816a
# ╠═93314cb4-5d95-4444-b4f2-f642721b0c81
# ╠═824e0e46-aa82-48fb-af76-9f2ff5131f70
# ╠═8ee4e9b1-669a-4cb9-96ca-f752224abfba
# ╠═b4b95535-62b5-4cfc-bff2-9a85d47e2ff9
# ╠═7daa3461-30ee-4982-9dad-3978b2928200
# ╠═a4eaf88b-3582-4d5c-8bed-81bdde06d03e
# ╠═d13d7550-7723-4a3c-a4c3-616b9bbf6f52
# ╠═0e64f4ea-33dd-4345-b896-c02b3c4bb9c3
# ╠═5a530357-ccfc-42d4-aaeb-e0ce6a7ba06f
# ╠═fc59f914-70b3-4800-b14f-cc72a3b594c5
# ╠═5607c259-2c5e-4e8d-a53c-d3250829dc89
# ╠═4014626f-0ae2-47fc-8af0-6d3cc92232fa
# ╠═239efc3c-df42-4202-b71b-ae6dd53c4e83

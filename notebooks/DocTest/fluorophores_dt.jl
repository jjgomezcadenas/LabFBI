### A Pluto.jl notebook ###
# v0.14.3

using Markdown
using InteractiveUtils

# ╔═╡ 59138a9e-5009-4334-9a9d-26da5a2efee2
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

# ╔═╡ fffac8bf-2bad-442b-918b-aa49c0b70901
using DrWatson

# ╔═╡ 968f3d06-a905-11eb-1252-4be9665e5088


# ╔═╡ a341cf04-626b-4cef-81a2-85cc8c059097
md"# Test and documentation of functions in module `fluorophores.jl`"


# ╔═╡ 03fcbd0f-42dd-4e32-a8a5-3b584b4710e6
@quickactivate "LabFBI"

# ╔═╡ 7752d3a0-3121-4d98-83ce-b5bcdc267a76
projectdir()

# ╔═╡ c714d318-9a92-4fbe-8922-a7db4d5cf3a0
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

# ╔═╡ 076677c8-ef19-479c-b80c-5aebb3415736
fluo = ingredients(srcdir("fluorophores.jl"))

# ╔═╡ 0a5d3226-97e5-4140-b951-892c9f2290ce
plt   = ingredients(srcdir("plotDF.jl"))

# ╔═╡ 462e91e3-62f8-45fb-986a-0abb9096d2b6
import Unitful:
    nm, μm, mm, cm, m, km, inch, ft, mi,
    ac,
    mg, g, kg,
    Ra, °F, °C, K,
    rad, °,
    ns, μs, ms, s, minute, hr, d, yr, Hz,
    eV,
    μJ, mJ, J,
	mW, W,
    A, N, mol, mmol, V, L, M

# ╔═╡ 85e2f9ae-fc05-45e3-9d82-d9f688dd2ba6
import PhysicalConstants.CODATA2018: N_A

# ╔═╡ 2d25fdfc-81c6-4102-929e-1bffc47312c0
md"## Notebook"

# ╔═╡ 3a1942f3-4739-405d-bde7-a8ba16b4e2fc
md"""
	struct Fluorophore

Represent a fluorescent molecule

\# Fields
- `name::String`                 : identifies the compound
- `expeak::Units (nm)`           : peak excitation
- `empeak::Units (nm)`           : peak emission
- `ϵ::Units(cm^{-1} M^{-1})` : molar extinction coefficient
- `Q::Float64`                   : Quantum efficiency
- `σ::Units(``cm^2``)`           : attenuation cross section
"""

# ╔═╡ a019b6df-4a3f-4b31-9422-633d0b44e008
fbi = fluo.Fluorophore("FBI", 325nm, 525nm, 13e+4*(M^-1*cm^-1), 0.6)

# ╔═╡ 85211dc2-5804-4424-94fa-9f5e1f61f429
@test unit(fbi.σ) == cm^2

# ╔═╡ 26aff41f-f3f9-46e0-9dec-1f783b36a775
md"""
	struct Solution

Represent a fluorophore in solution

\# Fields

- `name::String`  : identifies the compound
- `c::Units (M)`  : concentration
"""

# ╔═╡ c2172baa-7757-4ebc-a46e-f364bbdb3dd2
fbiStandard = fluo.Solution("FBI Standard", 5e-04M)

# ╔═╡ 12b94252-9f45-44cc-914e-b84b58382547
fbit = fluo.Solution("FBI Standard", 1M)

# ╔═╡ 3dc43b6e-7ab5-4de3-a8d3-f76a20fdc14e
@test uconvert(mol, fbit.c*1L) == 1mol

# ╔═╡ 482d89bf-5cd7-443c-bc00-eb101a6db9c3
md"""
	struct Powder

Represent a fluorophore in powder

\# Fields

- `name::String`         : identifies the compound
- `cs::Units (mmol/mg)`  : concentration relative to substrate
- `rs::Units (mg/``cm^2``)`  : concentration relative to area
"""

# ╔═╡ 76d2ff2e-e0de-4765-ac96-fd36c4c960fa
ps = fluo.Powder("Standard", 1mmol/mg, 1mg/cm^2)

# ╔═╡ b4546474-9f5f-4f18-81a3-11ea6a123662
md"""
	struct Pellet

A pellet of powder

\# Fields

- `area::Units(``cm^2``)`   : area of pellet
- `mass::Units (mg)`        : mass of pellet
"""

# ╔═╡ 9d7adad5-1b52-4df1-8d8b-c8012a2a525a
pe = fluo.Pellet(1μm^2, 1mg)

# ╔═╡ 7bc11b36-c8d7-458b-99d9-8516e2e67fde
md"""
	struct Mlayer

Represents a monolayer

\# Fields

- `name::String`              : name of monolayer (e.g, name fluorophores)
- `σm::Units (``mol/cm^2``)`  : number of molecules per area
"""

# ╔═╡ d7265e5f-9eb9-4679-8eda-965d5469c4e8
ml = fluo.Mlayer("test", 1e+6/μm^2)

# ╔═╡ bc43452a-16ac-4d10-8f49-abac34dbf8c8
@test pe.area * ml.σm ≈ 1e+6

# ╔═╡ 87478ec5-02d3-4e89-8287-882957d1e0ba
@test ps.cs * pe.mass ≈ 1mmol

# ╔═╡ a65c9fb2-da1c-4f15-ad35-cbc0cf127bf0
md"""
	nofv(sol::Solution, vol::Unitful.Volume)

Number of Fluorophores per Volume for a given solution

\# Fields

- `sol::Solution`        : A solution of fluorophores
- `vol::Unitful.Volume`  : The volume of solution considered
"""

# ╔═╡ ca45fe93-bf1a-45b7-be7f-9171d308fba1
fbiS= fluo.Solution("Test", 1M)

# ╔═╡ 4f4d98a2-c44a-4c76-bc2c-ffe065c9e8a2
@test fluo.nofv(fbiS, 1L) ≈ N_A*1.0*mol

# ╔═╡ cf925835-de35-40d4-aff1-47e35f03ee1f
md"""
	nofv(c::Quantity, vol::Unitful.Volume)

Number of Fluorophores per Volume for a given solution

\# Fields

- `c::(M)`               : A concentration of fluorophores
- `vol::Unitful.Volume`  : The volume of solution considered
"""

# ╔═╡ 879ea89d-56f6-4e31-8753-93eee3646421
@test fluo.nofv(1M, 1L) ≈ N_A*1.0*mol

# ╔═╡ 79c6ca52-0964-475e-828d-4502303893a1
md"""
	nofa(P::Powder, p::Pellet, unitarea::Unitful.Area)

Number of Fluorophores per area for a pellet p made of powder P

\# Fields

- `P::Powder`     : Define Powder
- `p::Pellet`     : Define Pellet
-`unitarea`       : unit area
"""

# ╔═╡ 0214487d-f44e-4cf4-98a1-446dca43ad43
pt = fluo.Pellet(1mm^2, 1mg)

# ╔═╡ 047bd1e4-e8df-4ccc-ade3-fba515ac7106
Pt2 = fluo.Powder("Powder test", 1mol/mg, 1mg/mm^2)

# ╔═╡ a6033c08-4d07-4c77-a2ea-5d35106eb3bb
@test fluo.nofa(Pt2, pt, 1mm^2) ≈ N_A*1.0*mol

# ╔═╡ 0dfce8ff-1862-4a7e-9e77-f108c12d58be
md"""
	nofa(P::Powder, area::Unitful.Area)

Number of Fluorophores per area for  powder P

\# Fields

- `P::Powder`         : Define Powder
- `area::(``cm^2``)`  : area
"""

# ╔═╡ ade5f991-6421-44f9-be48-bd7d5b0f07d8
@test fluo.nofa(Pt2,1mm^2) ≈ N_A*1.0*mol

# ╔═╡ 57c6ce77-2430-47d1-8f9d-5163750e9748
md"""
	fluorescence(f::Fluorophore, I::Quantity)

Number of  photons emitted per unit time when fluorosphore F
is illuminated with a laser of photon density I

\# Fields

- `f::Fluorophore`    : Define fluorophore
- `I::(``nγ/cm^2``)`  : Photon density
"""

# ╔═╡ e093a574-10e9-4c26-99bf-898a44a3e39a
fluo.fluorescence(fbi, 1Hz/cm^2)

# ╔═╡ Cell order:
# ╠═968f3d06-a905-11eb-1252-4be9665e5088
# ╠═59138a9e-5009-4334-9a9d-26da5a2efee2
# ╠═fffac8bf-2bad-442b-918b-aa49c0b70901
# ╠═a341cf04-626b-4cef-81a2-85cc8c059097
# ╠═03fcbd0f-42dd-4e32-a8a5-3b584b4710e6
# ╠═7752d3a0-3121-4d98-83ce-b5bcdc267a76
# ╟─c714d318-9a92-4fbe-8922-a7db4d5cf3a0
# ╠═076677c8-ef19-479c-b80c-5aebb3415736
# ╠═0a5d3226-97e5-4140-b951-892c9f2290ce
# ╠═462e91e3-62f8-45fb-986a-0abb9096d2b6
# ╠═85e2f9ae-fc05-45e3-9d82-d9f688dd2ba6
# ╠═2d25fdfc-81c6-4102-929e-1bffc47312c0
# ╟─3a1942f3-4739-405d-bde7-a8ba16b4e2fc
# ╠═a019b6df-4a3f-4b31-9422-633d0b44e008
# ╠═85211dc2-5804-4424-94fa-9f5e1f61f429
# ╟─26aff41f-f3f9-46e0-9dec-1f783b36a775
# ╠═c2172baa-7757-4ebc-a46e-f364bbdb3dd2
# ╠═12b94252-9f45-44cc-914e-b84b58382547
# ╠═3dc43b6e-7ab5-4de3-a8d3-f76a20fdc14e
# ╟─482d89bf-5cd7-443c-bc00-eb101a6db9c3
# ╠═76d2ff2e-e0de-4765-ac96-fd36c4c960fa
# ╟─b4546474-9f5f-4f18-81a3-11ea6a123662
# ╠═9d7adad5-1b52-4df1-8d8b-c8012a2a525a
# ╟─7bc11b36-c8d7-458b-99d9-8516e2e67fde
# ╠═d7265e5f-9eb9-4679-8eda-965d5469c4e8
# ╠═bc43452a-16ac-4d10-8f49-abac34dbf8c8
# ╠═87478ec5-02d3-4e89-8287-882957d1e0ba
# ╟─a65c9fb2-da1c-4f15-ad35-cbc0cf127bf0
# ╠═ca45fe93-bf1a-45b7-be7f-9171d308fba1
# ╠═4f4d98a2-c44a-4c76-bc2c-ffe065c9e8a2
# ╟─cf925835-de35-40d4-aff1-47e35f03ee1f
# ╠═879ea89d-56f6-4e31-8753-93eee3646421
# ╟─79c6ca52-0964-475e-828d-4502303893a1
# ╠═0214487d-f44e-4cf4-98a1-446dca43ad43
# ╠═047bd1e4-e8df-4ccc-ade3-fba515ac7106
# ╠═a6033c08-4d07-4c77-a2ea-5d35106eb3bb
# ╟─0dfce8ff-1862-4a7e-9e77-f108c12d58be
# ╠═ade5f991-6421-44f9-be48-bd7d5b0f07d8
# ╠═57c6ce77-2430-47d1-8f9d-5163750e9748
# ╠═e093a574-10e9-4c26-99bf-898a44a3e39a

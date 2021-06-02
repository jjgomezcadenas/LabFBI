#include("dffunctions.jl")
using Unitful
using UnitfulEquivalences
using PhysicalConstants.CODATA2018
using DataFrames

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

import PhysicalConstants.CODATA2018: N_A

#TYPES
"""
	struct Fluorophore

Represent a fluorescent molecule

# Fields
- `name::String`                 : identifies the compound
- `expeak::Units (nm)`           : peak excitation
- `empeak::Units (nm)`           : peak emission
- `ϵ::Units(``cm^{-1} M^{-1}``)` : molar extinction coefficient
- `Q::Float64`                   : Quantum efficiency
- `σ::Units(``cm^2``)`           : attenuation cross section
"""
struct Fluorophore
    name::String
    expeak::typeof(1.0nm)
    empeak::typeof(1.0nm)
    ϵ::typeof(1.0/(cm*M))
    Q::Float64
	σ::typeof(1.0cm^2)
	function Fluorophore(name, ex, en, ϵ, Q)
		σ = log(10) * uconvert(cm^2/mol, ϵ) / N_A
		new(name, ex, en, ϵ, Q, σ)
	end
end

"""
	 FbiFluorophores
Defines Fbi and FbiBa response to different lasers

#Fields
- `fbi   ::Dict{String, Fluorophore}`
- `fbiba ::Dict{String, Fluorophore}`

"""
struct FbiFluorophores
	fbi   ::Dict{String, Fluorophore}
	fbiba ::Dict{String, Fluorophore}
end


"""
	struct Solution

Represent a fluorophore in solution

# Fields

- `name::String`  : identifies the compound
- `c::Units (M)`  : concentration
"""
struct Solution
    name::String
    c::typeof(1.0M)
end


"""
	struct Powder

Represent a fluorophore in powder

# Fields

- `name::String`         : identifies the compound
- `cs::Units (mmol/mg)`  : concentration relative to substrate
- `rs::Units (mg/``cm^2``)`  : concentration relative to area
"""
struct Powder
    name::String
    cs::typeof(1.0mmol/mg)
    rs::typeof(1.0mg/cm^2)
end


"""
	struct Pellet

A pellet of powder

# Fields

- `area::Units(``cm^2``)`   : area of pellet
- `mass::Units (mg)`        : mass of pellet
"""
struct Pellet
    area::typeof(1.0mm^2)
    mass::typeof(1.0mg)
end


"""
	struct Mlayer

Represents a monolayer

# Fields

- `name::String`              : name of monolayer (e.g, name fluorophores)
- `σm::Units (``mol/cm^2``)`  : number of molecules per area
"""
struct Mlayer
    name::String
    σm::typeof(1.0/mm^2)

end


#FUNCTIONS

"""
fbi_fluorophores(adf::DataFrame, gs::Vector{String}=["g1", "g2"],
						  ls::Vector{Int64} = [325,405],
						  qs::Vector{Float64} = [0.67,0.67])

Return an object of type FbiFluorophore, defining the fluorophores associated
with the versions of the molecule and the excitation wavelengths.

# Fields
- `adf::DataFrame`    :  A DF that contains data on extinction coefficients
- `gs::Vector{String}`:  Version of molecule (g1, g2...)
- `ls::Vect{Int64}`   :  Wavelengths (in nm)
- `qs::Vect{Float64}` :  Quantum efficiencies

"""
function fbi_fluorophores(adf::DataFrame, gs::Vector{String}=["g1", "g2"],
	                      ls::Vector{Int64} = [325,405],
						  qs::Vector{Float64} = [0.67,0.67])
	ffbi = Dict()
	ffbiba = Dict()

	for gn in gs
		efbi = string("ϵFbi",titlecase(gn))
		efbiba = string("ϵFbiBa",titlecase(gn))
		for l in ls
			lfbi = string("l", l, gn)
			tfbi = string("Fbi",l,gn)
			tfbiba = string("FbiBa",l,gn)

			ffbi[lfbi]   = Fluorophore(tfbi, l*nm, 489nm,
						   select_element(adf, "λ", l, efbi)/(M*cm), qs[1])


			ffbiba[lfbi] = Fluorophore(tfbiba, l*nm, 405nm,
						   select_element(adf, "λ", l, efbiba)/(M*cm), qs[2])
		end
	end

	return FbiFluorophores(ffbi, ffbiba)
end


"""
	nofv(sol::Solution, vol::Unitful.Volume)

Number of Fluorophores per Volume for a given solution

# Fields

- `sol::Solution`        : A solution of fluorophores
- `vol::Unitful.Volume`  : The volume of solution considered
"""
function nofv(sol::Solution, vol::Unitful.Volume)
	return N_A * uconvert(mol, sol.c * vol)
end


"""
	nofv(c::Quantity, vol::Unitful.Volume)

Number of Fluorophores per Volume for a given solution

# Fields

- `c::(M)`               : A concentration of fluorophores
- `vol::Unitful.Volume`  : The volume of solution considered
"""
function nofv(c::Quantity, vol::Unitful.Volume)
	return N_A * uconvert(mol, c * vol)
end


"""
	nofa(P::Powder, p::Pellet, unitarea::Unitful.Area)

Number of Fluorophores per area for a pellet p made of powder P

# Fields

- `P::Powder`     : Define Powder
- `p::Pellet`     : Define Pellet
-`unitarea`       : unit area
"""
function nofa(P::Powder, p::Pellet, unitarea::Unitful.Area)
	return N_A * uconvert(mol, unitarea* p.mass * P.cs/p.area)
end


"""
	nofa(P::Powder, area::Unitful.Area)

Number of Fluorophores per area for  powder P

# Fields

- `P::Powder`         : Define Powder
- `area::(``cm^2``)`  : area
"""
function nofa(P::Powder, area::Unitful.Area)
	return N_A* uconvert(mol, area * P.cs * P.rs)
end


"""
	nofa(M::Mlayer, area::Unitful.Area)

Number of Fluorophores per area for monolayer M

# Fields

- `M::Mlayer`         : Define monolayer
- `area::(``cm^2``)`  : area
"""
function nofa(M::Mlayer, area::Unitful.Area)
	return  uconvert(NoUnits, area * M.σm)
end

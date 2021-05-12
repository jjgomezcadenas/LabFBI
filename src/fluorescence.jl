include("fluorophores.jl")
include("setup.jl")
using Unitful
using UnitfulEquivalences

"""
	Fluorescence
Defines the fluorescence reponse of Fluorophores to different lasers
# Fields
- `fbi::Dict{String, typeof(1.0Hz)}`
- `fbiba::Dict{String, typeof(1.0Hz)}`

"""
struct Fluorescence
	fbi::Dict{String, typeof(1.0Hz)}
	fbiba::Dict{String, typeof(1.0Hz)}
end

"""
	fluorescence(f::Fluorophore, I::Quantity)

Number of  photons emitted per unit time when fluorosphore F
is illuminated with a laser of photon density I

# Fields

- `f::Fluorophore`    : Define fluorophore
- `I::(``nγ/cm^2``)`  : Photon density
"""
function fluorescence(f::Fluorophore, I::Quantity)
	return uconvert(Hz, f.σ * f.Q * I)
end


"""
	emitted_fluorescence(fbif::FbiFluorophores, lT::LaserSetup,
						 gs::Vector{String}= ["g1", "g2"],
						 ls::Vector{Int64} = [325,405])

Return emitted fluorescence, given a set of fluorophores and a laser setup.

# Fields

- `fbif::FbiFluorophores`: Fluorophores setup
- `lT::LaserSetup` : Laser setup
- `gs::Vector{String}´: molecule versions
- `ls::Vector{Int64}´ : laser wavelengths
"""
function emitted_fluorescence(fbif::FbiFluorophores, lT::LaserSetup,
	                          gs::Vector{String}= ["g1", "g2"],
						      ls::Vector{Int64} = [325,405])

	dfbi = Dict{String, typeof(1.0Hz)}()
	dfbiba = Dict{String, typeof(1.0Hz)}()
	for gn in gs
		for l in ls
			klt = string("l",l)
			kff = string(klt,gn)
			dfbi[kff] =   fluorescence(fbif.fbi[kff], lT.Is[klt])
			dfbiba[kff] = fluorescence(fbif.fbiba[kff], lT.Is[klt])
		end
	end

	return Fluorescence(dfbi, dfbiba)
end

"""
	molecules_in_fov(lT::LaserSetup, solfbi::Solution)

Return the number of molecules in the field of view (fov) defined by the
laser setup.

# Fields

- `lT::LaserSetup`   : Laser setup (may include several lasers)
- `solfbi::Solution` : Solution of molecules (defined by a given concentration)

"""

function molecules_in_fov(lT::LaserSetup, solfbi::Solution)
	mfov = Dict{String, Float64}()
	for lname in keys(lT.lasers)
		mfov[lname] = nofv(solfbi, lT.fovs[lname].v)
	end
	return mfov
end

"""
	emitted_fluorescence_fov(efpm::lfi.LabFbi.Fluorescence,
		                     mfov::Dict{String, Float64},
		                     gs::Vector{String}= ["g1", "g2"],
						     ls::Vector{Int64} = [325,405])

Return the emitted fluorescence in the FoV defined by the laser.
This equals the product of the number of molecules in that FoV and
the fluorescence per molecule.

# Fields

- `efpm::Fluorescence` : emitted fluorescence per molecule
- `mfov::Dict{String, Float64}`: number of molecules in the fov
- `gs::Vector{String}`: Versions of FBI molecule
- `ls::Vector{Int64}`: Wavelengths of lasers

"""
function emitted_fluorescence_fov(efpm::Fluorescence,
		                          mfov::Dict{String, Float64},
		                          gs::Vector{String}= ["g1", "g2"],
						          ls::Vector{Int64} = [325,405])

	dfbi   = Dict{String, typeof(1.0Hz)}()
	dfbiba = Dict{String, typeof(1.0Hz)}()
	for gn in gs
		for l in ls
			klt = string("l",l)
			kff = string(klt,gn)
			dfbi[kff]   =   efpm.fbi[kff] * mfov[klt]
			dfbiba[kff] =   efpm.fbiba[kff] * mfov[klt]
		end
	end

	return Fluorescence(dfbi, dfbiba)
end

"""
	filtered_fluorescence(fl::Fluorescence, filter::Float64,
		                  gs::Vector{String}= ["g1", "g2"],
						  ls::Vector{Int64} = [325,405])

Return the fluorescence that passes through filter f.

# Fields

- `fl::Fluorescence` : Fluorescence before filter
- `filter::Float64`  : Filtering value (a fraction of 1.)
- `gs::Vector{String}`: Versions of FBI molecule
- `ls::Vector{Int64}`: Wavelengths of lasers

"""
function filtered_fluorescence(fl::Fluorescence, filter::Float64,
		                       gs::Vector{String}= ["g1", "g2"],
						       ls::Vector{Int64} = [325,405])

	dfbi   = Dict{String, typeof(1.0Hz)}()
	dfbiba = Dict{String, typeof(1.0Hz)}()
	for gn in gs
		for l in ls
			kff = string("l",l,gn)
			dfbi[kff]   =   fl.fbi[kff] * filter
			dfbiba[kff] =   fl.fbiba[kff] * filter
		end
	end

	return Fluorescence(dfbi, dfbiba)
end

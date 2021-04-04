### A Pluto.jl notebook ###
# v0.14.0

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 4f35dcfc-8644-11eb-1c86-75608b38bf5f
begin
	using Pkg
	Pkg.add.(["Unitful", "UnitfulEquivalences", "Plots","Interpolations"])
end

# ╔═╡ 2e3c7fda-904f-11eb-2988-25604b1caad0
begin
	Pkg.add("PhysicalConstants")
end

# ╔═╡ 4ea80a97-44a6-4205-ae30-eaf239309313
Pkg.add.(["Images", "ImageIO", "ImageView"])

# ╔═╡ 3207f446-8643-11eb-37ba-c9aec47fcb8f
begin
	using Markdown
	using InteractiveUtils
	using PlutoUI
	using Test
	using Plots
	using LsqFit
	using Interpolations
	using Images
	using ImageIO
	using ImageView
	using CSV
	using DataFrames
end

# ╔═╡ 5115917a-8644-11eb-19fc-0528741ca75d
begin
	using Unitful
	using UnitfulEquivalences
end

# ╔═╡ 8ba78a98-904f-11eb-1977-2750643d2f9f
using PhysicalConstants.CODATA2018

# ╔═╡ f2faa594-9316-11eb-2451-33e4b4d3a4f5
using StatsPlots

# ╔═╡ 5356699e-93ca-11eb-0afd-e7d0edc4231e
using QuadGK

# ╔═╡ fc8b79a2-8728-11eb-2da7-e3ffa3ceef08
using DrWatson

# ╔═╡ 79cfd2fc-9046-11eb-2b13-1b877d57645d
md"# LabFBI"

# ╔═╡ 0f2f4c78-8729-11eb-2bab-27812ce8c47e
@quickactivate "LabFBI"

# ╔═╡ 5ee27d52-86fd-11eb-365e-9f2b1a095575
;pwd()

# ╔═╡ 621ec96c-86fd-11eb-1c41-379cc17180dc
datadir()

# ╔═╡ 0b1cedc7-8ada-45a5-adef-fbae794dee3e
markercolors = [:green :orange :black :purple :red  :yellow :brown :white]

# ╔═╡ 70525ca4-8643-11eb-121e-196d80a60fcf
md"## Laboratory Setup"

# ╔═╡ e5e02216-86fb-11eb-1cf6-6f6962313a09
load("laserFBI.png")

# ╔═╡ 3c01d7c8-8645-11eb-372d-4d66ae46ae72
md"## Units"

# ╔═╡ 46b5465a-8645-11eb-0291-612455795518
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

# ╔═╡ 443e90e2-8664-11eb-19bc-51f8cfc282d5
md"- To convert a number to the desired units simply multiply by the unit"

# ╔═╡ 33ad8312-8661-11eb-1909-87eb9ed385bc
begin
	x  = 3.0         # Float65
	ux = x*mm        # Units of length
	y  = ux/mm       # no units again
end

# ╔═╡ e6e6be0a-9046-11eb-01ae-c572f38fca9e
@test typeof(x) == Float64

# ╔═╡ ebd713a6-9046-11eb-28f3-c14460ab4b0b
typeof(ux)

# ╔═╡ f7a2b46a-9046-11eb-04b7-2971adb83918
@test typeof(y) == Float64

# ╔═╡ 4c8ab940-8661-11eb-30a9-9301eb73ea33
md"- We can convert between equivalent energy systems (e.g, SI to natural units). For example, the wavelength of a phothon of 1 eV and viceversa"

# ╔═╡ c1615b46-8663-11eb-21dc-5b746b6ae780
uconvert(nm, 1.0eV, Spectral())

# ╔═╡ a73a6004-8664-11eb-0cd2-6bbf19682226
uconvert(eV, 1.0nm, Spectral())

# ╔═╡ 0f9a7e86-8665-11eb-0810-43b5d2227a40
md"- We can use ustrip to strip the units (also divide by the unit)"

# ╔═╡ 0374c94a-8665-11eb-28b9-8127833558fc
ustrip(1.0nm)

# ╔═╡ 341791a4-665f-4bc1-bc3d-35799fd4afec
md"- To compute the relation between units of the same dimension (e.g, nm and cm) one uses uconvert with NoUnits"

# ╔═╡ 39d064d4-93cd-11eb-0530-7d4dd6b3c2d5
xxt = 1cm

# ╔═╡ 414f254c-93cd-11eb-07ca-f323b0b1b36c
uconvert(Unitful.NoUnits, xxt/nm)

# ╔═╡ 754ab06a-9048-11eb-3953-bd88a367d21e
md"- We can use non SI units such as liter"

# ╔═╡ 4dc66f48-9048-11eb-2334-6d9115a0a741
typeof(L)

# ╔═╡ 658252da-904a-11eb-15f0-036d2412b386
md"- mol"

# ╔═╡ 58bfee7c-904a-11eb-07f8-cdafd513ac04
typeof(mol)

# ╔═╡ a607fc3a-9048-11eb-2074-f188307b35a9
md" - Molarity is also defined"

# ╔═╡ 537c6c58-9048-11eb-0176-97f2bb8e06f7
typeof(M)

# ╔═╡ 6dd9d1b0-904a-11eb-1fec-9fe82cf1e5ee
typeof(mol/L)

# ╔═╡ 8dd2c6a4-904a-11eb-28f9-8f9ede41f1dd
@test dimension(mol/L) == dimension(M)

# ╔═╡ 05ad4620-904b-11eb-3f68-13a12b8b5f27
md"- The molar attenuation coefficient is a measurement of how strongly a chemical species attenuates light at a given wavelength. It is an intrinsic property of the species. Attenuation coefficient is tusually expressed in terms of M−1⋅cm−1. The molar attenuation coefficient is also known as the molar extinction coefficient and molar absorptivity."

# ╔═╡ a04ab7ba-904a-11eb-21b9-8b09059f1ca5
dimension(1.0/(cm*M))

# ╔═╡ ee7f63de-904f-11eb-0e8f-1f8ec9d50811
md"## Physical constants"

# ╔═╡ 631833b8-0649-4f45-b1dd-e4a53091d321
md"#### import Avogadro Number"

# ╔═╡ f8d3f48a-904f-11eb-2659-51f7508b434d
import PhysicalConstants.CODATA2018: N_A

# ╔═╡ 9fdae9e2-6d2e-4901-8f48-14764b6800c2
N_A

# ╔═╡ 42d84ad0-904c-11eb-08d6-7b4b117bf2a6
md"## Types and functions"

# ╔═╡ 5f56cf78-5b2d-471b-b21f-f242b7d78175
md" - A (Fluorescent) Molecule is defined by the excitation and emission peaks, its molar attenuation coefficient and its quantum efficieny. The attenuation cross section can be directly computed from the molar attenuation coefficient"

# ╔═╡ 3f63bf34-9046-11eb-07a2-6b5f90371cb6
struct Molecule
        name::String
        ex::typeof(1.0nm)      # peak excitation
        en::typeof(1.0nm)      # peak emission
        ϵ::typeof(1.0/(cm*M))  # molar attenuation coefficient
        Q::Float64             # Quantum efficiency
	    σ::typeof(1.0cm^2)     # attenuation cross section
	function Molecule(name, ex, en, ϵ, Q)
		σ = log(10) * uconvert(cm^2/mol, ϵ) / N_A
		new(name, ex, en, ϵ, Q, σ)
		end
end

# ╔═╡ a10b23d9-e366-4e38-80b6-a1532e30302e
md"- A solution is defined by the concentration of solute"

# ╔═╡ a90d8eea-904b-11eb-2c8d-b9fbc016aadb
struct Solution
    name::String
    c::typeof(1.0M) # concentration
end

# ╔═╡ 7de08c7a-7c93-45af-b46a-247cfc30c859
md" - Powder is defined by the concentration of solute with respect to support and the concentration per area of support"

# ╔═╡ 1694d934-904c-11eb-0b79-5f88092d69ce
struct Powder
    name::String
    cs::typeof(1.0mmol/mg) # concentration to substrate
    rs::typeof(1.0mg/cm^2) # concentration per area
end

# ╔═╡ ea7a474e-0d12-436a-b8a6-e0737507d921
md" - A pellet is defined by its area and its mass"

# ╔═╡ d6f6bfda-904c-11eb-31d1-3550103982f9
struct Pellet
    area::typeof(1.0mm^2)
    mass::typeof(1.0mg)
end

# ╔═╡ 15c9d8e0-904e-11eb-35e4-3d920e199e97
function nof_molecules_area(P::Powder, p::Pellet, unitarea)
	return uconvert(unitarea,
		            N_A * uconvert(mol, p.mass * P.cs)/p.area)
end

# ╔═╡ bf51632d-4543-42b6-b9d1-f54c240412d4
md"- Given a powder p, one can compute the number of molecules per unit area"

# ╔═╡ dece508c-9052-11eb-3d03-431affa529dc
function nof_molecules_area(P::Powder, unitarea)
	return N_A* uconvert(mol*unitarea, P.cs * P.rs)
end

# ╔═╡ 84e3a8e1-665b-4302-9f9d-a1fb82c13aa2
md"- The simplest laser is defined by its wavelength and its power"

# ╔═╡ a61dcc86-8643-11eb-274a-db5011271345
struct Laser
	λ::typeof(1.0nm)
	P::typeof(1.0mW) # power
end

# ╔═╡ 52d85aa0-feaa-40aa-91a2-93b1b3534da6
md" - A field of view (FOV) is defined by its diameter and thickness and has computable attributes such as area and volume"

# ╔═╡ 2287af6e-9128-11eb-3070-f7af71689367
struct FoV
    d::Unitful.Length
    z::Unitful.Length
	a::Unitful.Area
    v::Unitful.Volume

	function FoV(d,z)
		a = π * (d/2.)^2
		v = a * z
		new(d,z,a,v)
	end
end

# ╔═╡ 792cfd0e-89e0-478e-abc1-ee8fcbba79e8
md"- We can define a convenience function that related the wavelength of a photon and its energy"

# ╔═╡ 3adf7eba-865c-11eb-2cb2-01884cac9b61
function photon_energy(λ::Unitful.Length)
	uconvert(eV, λ, Spectral())
end

# ╔═╡ ae7ed83e-562a-4e13-a24e-a8e4c793a27b
md"- The number of photons (per unit time) is a characteristic of a laser"

# ╔═╡ 54c84cf6-908b-11eb-26e7-3d833153a5bc
n_photons(l::Laser) = uconvert(Hz, l.P / photon_energy(l.λ))

# ╔═╡ dc537c9f-3764-4055-a9bc-9eacdaf81fbb
md"- This can also be computed directly from wavelength and power"

# ╔═╡ 7b4761e4-908b-11eb-2840-293f5286d1d5
n_photons(λ::Unitful.Length, p::Unitful.Power) = uconvert(Hz, p / photon_energy(λ))

# ╔═╡ 1810e902-9ea4-41da-a7fc-60e03197bc7e
md" - The number of photons in a given time unit is also a function of laser"

# ╔═╡ 89528344-0c5c-41d9-a1d5-62b315b90ef7
md"- Given Length, power and area is possible to compute photon density"

# ╔═╡ 58f59d86-5a11-4dfd-952d-d3d5655999f7
md"- Also given a laser and a fov"

# ╔═╡ 0630250e-8602-4556-8948-1f6858e0b587
md"- Delievered energy is a function of laser and time"

# ╔═╡ 9f25699e-8665-11eb-13af-0be83c293556
function delivered_energy(laser::Laser, t::Unitful.Time)
	laser.P * t
end

# ╔═╡ 6a28f556-8667-11eb-3614-d59ecab119d6
function n_photons(laser::Laser, t::Unitful.Time)
	uconvert(eV,delivered_energy(laser, t)) / photon_energy(laser.λ)
end

# ╔═╡ 811b7904-90bd-11eb-3182-4969edbde50e
photon_density(λ::Unitful.Length, p::Unitful.Power, a::Unitful.Area) = n_photons(λ, p)/ a

# ╔═╡ 28d38d4c-9129-11eb-310e-1996bd60a771
photon_density(l::Laser, fov::FoV) = n_photons(l) / fov.a

# ╔═╡ 815a3ba6-4a61-4d96-9aca-d32318c7a99f
md" - One can compute the fluorescence (number of photons emitted per unit time) of a molecule that is illuminated with a laser of photon density I"

# ╔═╡ 9596d471-bf22-47f3-8f2d-577b8d875382
function fluorescence_per_molecule(m::Molecule, I::Quantity)
	return uconvert(Hz, m.σ * m.Q * I)
end


# ╔═╡ cb4d9e52-9134-11eb-2557-d94ef2d4e541
md"function `geometrical_acceptance(d,D)` computes the fraction of photons that make it through a lens or hole of diameter D located at a distance d from the emission point (no focusing)"

# ╔═╡ acba0a22-9134-11eb-3bd8-a7e99e8a710c
geometrical_acceptance(d::Float64, D::Float64) = 0.5(1. - d/sqrt(d^2 + (D/2.)^2))

# ╔═╡ 45bdd766-9131-11eb-31dc-f9ee1ced0db9
md"# Notebook"

# ╔═╡ 87199dd6-9047-11eb-312d-1fea3493dd46
md"### Molecule"

# ╔═╡ c3b5e2a0-90bc-11eb-243e-eb6908eeeb13
md"#### FBI and FBI-Ba2+ in solution"

# ╔═╡ 304c3140-90bc-11eb-1bbd-632bdcce5ced
mfbisol250 = Molecule("FBI-solution", 250nm, 489nm, 11260/(M*cm), 0.67)

# ╔═╡ 6ff7fe82-90bc-11eb-3b33-271f008fe4b4
mfbibasol250 = Molecule("FBI-Ba2-solution", 250nm, 428nm, 8060/(M*cm), 0.45)

# ╔═╡ 538eccf0-904c-11eb-0a29-491192a77e31
md"#### Solution and silica powder used commonly in LabFBI"

# ╔═╡ b7de26de-904b-11eb-223d-b9d5f33b1c5f
solfbi = Solution("FBI-Standard-Solution", 5e-5M)

# ╔═╡ e82978a2-904b-11eb-15af-8f1fedc428a5
solfbiba = Solution("FBIBa-Standard-Solution", 5e-5mol/L)

# ╔═╡ 075148ae-904c-11eb-1969-1be58cdfb088
@test solfbi.c == solfbiba.c

# ╔═╡ 7a93e858-904c-11eb-0a45-17e5bb3d405a
sifbi = Powder("FBI-Standard", 2.25e-5mmol/mg, 50mg/cm^2)

# ╔═╡ 91d4f7c8-904c-11eb-0486-23d11f6936df
sifbiba = Powder("FBIBa-Standard", 7.38e-8mmol/mg, 50mg/cm^2)

# ╔═╡ 934c7fcc-904c-11eb-1339-d5a1e82ac071
hsifbi = Powder("FBI-Homeopatic", 7.4e-15mmol/mg, 50mg/cm^2)

# ╔═╡ 0773cad6-904d-11eb-1ba0-9b46e5d1046c
hsipellet = Pellet(π*5*mm^2, 36.74mg)

# ╔═╡ c9017918-9050-11eb-2b1b-a9d8243e00a7
md" #### Number of molecules in a homeophatic pellet of FBI powder"

# ╔═╡ f2d90250-9054-11eb-3008-01b74530d71e
md"- Computing number of molecules per area in terms of pellet parameters for homeopatic pellet"

# ╔═╡ 9b8f1f65-8476-4a2c-adc2-b1973776110e
@test nof_molecules_area(hsifbi, hsipellet, μm^-2) ≈ uconvert(μm^-2, N_A * uconvert(mol, hsipellet.mass * hsifbi.cs)/hsipellet.area)

# ╔═╡ 50c68b42-9052-11eb-0dd3-41e31d7a1947
nof_molecules_area(hsifbi, hsipellet, μm^-2)

# ╔═╡ 07fdc33c-9055-11eb-2f90-b9077aa5fe3d
md"- Number of molecules per area from powder parameters for homeopatic solution"

# ╔═╡ f63079f2-9053-11eb-0e1b-cd8542581fc2
nof_molecules_area(hsifbi, μm^-2)

# ╔═╡ 2f084308-9055-11eb-2ce2-bd48996c0f5e
md"- Number of molecules per mm^2 for standard FBI powder"

# ╔═╡ 3dbd2044-9055-11eb-3ca4-75af1ab90b72
nof_molecules_area(sifbi, mm^-2)

# ╔═╡ 66a99aaa-9055-11eb-1126-57af32a69b8e
md"- Number of molecules per area for standard FBIBa2+ powder"

# ╔═╡ 73885928-9055-11eb-20a0-41f136464b4f
nof_molecules_area(sifbiba, mm^-2)

# ╔═╡ cadd0620-9315-11eb-228e-135bd066038d
md"### Molecule spectra"

# ╔═╡ d65c78aa-9315-11eb-094c-01765ce6cc2f
begin
	path=string(datadir(),"/fluorimeter")
	selFbiBa =["EDI044_FBI_Ba_round4_em325_375_405.csv",
             "EDI044_FBI_Ba_round5_em325_375_405.csv"]
	selFbi =["EDI044_FBI_round4_em325_375_405.csv",
               "EDI044_FBI_round5_em325_375_405.csv"]

	fsfbi =   [string(path,"/", fbi) for fbi in selFbi]
    fsfbiba = [string(path,"/", fbi) for fbi in selFbiBa]

	csvSelFbi = CSV.File(fsfbi[1]; delim=';', decimal=',')
	fbidf = DataFrame(csvSelFbi);
	csvSelFbiBa = CSV.File(fsfbiba[1]; delim=';', decimal=',')
	fbibadf = DataFrame(csvSelFbiBa);
end

# ╔═╡ 8ad489e2-9318-11eb-1531-7515b778e6b8
begin
    wf = 335.
    we = 701.
    ws = 2.
    wl = wf:ws:we
end

# ╔═╡ 0eb3c064-9316-11eb-3c97-afba4053eb94
function plot_fbi(fbidf, fbibadf, fbi375, fbiba375, w)
	L =collect(wl)
	FBI = fbi375.(L)
	FBIBA = fbiba375.(L);
	@df fbidf plot(:λ, :E375,
               colour = :green,
               shape = :circle,
               label  ="FBI",
               legend=:topright)
	plot!(L, FBI,
    colour = :green,
    label  ="FBI Interpol",
    legend=:topright)

	@df fbibadf plot!(:λ, :E375,
               colour = :blue,
               shape = :circle,
               label  ="FBI+Ba2+",
               legend=:topright)

	plot!(L, FBIBA,
    colour = :blue,
    label  ="FBI+Ba2+ Interpol",
    legend=:topright)


	xlabel!("λ (nm)")
	ylabel!("Intensity(AU)")
	title!("FBI/BA response, 375 nm")
end

# ╔═╡ 322216ee-93cb-11eb-272c-b7f5e5d87ee5
function plot_fbin(fbin, fbiban,w)
	L =collect(wl)
	FBI = fbin.(L)
	FBIBA = fbiban.(L);
	plot(L, FBI,
    colour = :green,
    label  ="FBI Normalised",
    legend=:topright)

	plot!(L, FBIBA,
    colour = :blue,
    label  ="FBI+Ba2+ Normalised",
    legend=:topright)


	xlabel!("λ (nm)")
	ylabel!("Intensity(normalised)")
	title!("FBI/BA response, 375 nm")
end

# ╔═╡ 9989f79a-9318-11eb-07a3-a9986579342d
begin
	fbi375 = LinearInterpolation(wl, fbidf[!,"E375"]);
	fbiba375 = LinearInterpolation(wl, fbibadf[!,"E375"]);
end

# ╔═╡ 5361bf8e-9318-11eb-1acd-e9a10c2b30bf
plot_fbi(fbidf, fbibadf,fbi375,fbiba375,wl)

# ╔═╡ 2ecbd182-93ca-11eb-08f6-3d0c1710f9f5
md"#### Normalize to the area of fbi and fbiba"

# ╔═╡ eba136e2-93c9-11eb-21f8-8b36a5158624
fbi375N, eps = quadgk(fbi375, 375., 700.)

# ╔═╡ 6deec574-93ca-11eb-0b5f-1b0d7ca04b61
fbiba375N, eps2 = quadgk(fbiba375, 375., 700.)

# ╔═╡ 46411c0e-93cd-11eb-1039-8fe43a7679d5
md"`fpdf(fi, xmin, xmax, N)` takes an interpolated function fi, lower (xmin) and upper (xmax) values of interpolation and a normalisation constant and returns a function equal to fi(x) / N in the range xmin-xmax and zero otherwise."

# ╔═╡ 353b7b37-c037-49f9-b149-56541db6cc28
function gfpdf(fi, xmin, xmax, N)
	function fn(x)
		if x < xmin || x > xmax
			return 0
		else
			return fi(x) / N
		end
	end

	return fn
end

# ╔═╡ 566a984a-f0fc-482f-9b04-029cba608aef
fpdf(fi, xmin::Float64, xmax::Float64, N::Float64) = gfpdf(fi, xmin, xmax, N)

# ╔═╡ fee9eace-93cd-11eb-114c-c174cc6637ce
fbipdf = fpdf(fbi375, 375.0, 701.0, fbi375N)

# ╔═╡ 4d17046e-ce2e-4863-8231-2ded5722ff44
fbibapdf = fpdf(fbiba375, 375.0, 701.0, fbiba375N)

# ╔═╡ fc40d4e4-93ce-11eb-02ed-57cb27c0db82
LL =collect(wl)

# ╔═╡ f06bb3a3-247c-4985-91cf-8ce79744945e
length(LL)

# ╔═╡ 2d0e3c95-fcd6-4d96-9ca0-bd8c8be01edb
Ifbipdf, eps3 = quadgk(fbipdf, 375., 700.)

# ╔═╡ 65dbf890-3d5b-4f1c-aa83-4d72556f9173
@test Ifbipdf ≈ 1.0

# ╔═╡ c1f8e6d0-52b2-4ca7-8aaf-6ed25bbfce08
Ifbibapdf, eps4 = quadgk(fbibapdf, 375., 700.)

# ╔═╡ 0e678e6b-dcbd-4cd2-9be9-5f9decc2dcc8
@test Ifbibapdf ≈ 1.0

# ╔═╡ e4e97bd1-28c2-4340-ade7-50dab201fe1a
fbi375N * fbipdf(500.0)

# ╔═╡ b3287b7b-3a7f-4ef9-aaf1-4ed2754b7e72
fbi375(500.0)

# ╔═╡ 56e1f815-7337-499d-a4bc-ed65121984c5
F1 = fbi375N * fbipdf.(LL)

# ╔═╡ b55a5691-97b0-470e-8fd3-b143c26f724d
F2 = fbi375.(LL)

# ╔═╡ 054700b9-a1a6-4c04-94b7-037df2a890c4
@test F1 ≈ F2

# ╔═╡ 76a569c4-93cb-11eb-3bb4-55e0c9221e53
plot_fbin(fbipdf, fbibapdf, wl)

# ╔═╡ 9bf0b7e6-8643-11eb-142b-7fc904562cc8
md"### Laser"

# ╔═╡ e9076d72-8648-11eb-0c02-9b00b0605e72
l400P1mw = Laser(405.0nm, 1.0mW)

# ╔═╡ 2aef669a-864b-11eb-3290-359ff40cece4
@test l400P1mw.λ == 405.0nm

# ╔═╡ 4bb9cb3e-864b-11eb-3b51-eb59b6da2c91
@test l400P1mw.P == 1.0mW

# ╔═╡ 1947e7a4-908b-11eb-0477-99f5152c6526
begin
	Λ =collect(170:1:800) * nm
	PE = photon_energy.(Λ)
	plot(Λ./nm, PE./eV, leg=false,lw=2)
	xlabel!("λ (nm)")
	ylabel!("Photon energy (eV)")
	title!("Photon energy as a function of λ")
end

# ╔═╡ 2aff4f8c-908b-11eb-29af-93923b73c692
md"#### Number of photons per unit time as a function of power"

# ╔═╡ 6ce764b6-908b-11eb-3013-0f904fb0a0d1
n_photons(l400P1mw)

# ╔═╡ 9ab0e624-908b-11eb-2d84-d93e1317e3bd
n_photons(405nm, 1mW)

# ╔═╡ a773efa0-908b-11eb-12a4-7bce3bb0b3d6
begin
	P =collect(1:1:10^3) * mW
	NP = n_photons.((405nm,), P )
	plot(P./mW, NP./Hz, leg=false,lw=2)
	xlabel!("P (mW)")
	ylabel!("number of photons (Hz)")
	title!("Number of photons as a function of P")
end

# ╔═╡ bb4fa686-908b-11eb-0aa0-27ab50bb0a4e
md"#### delivered energy"

# ╔═╡ 1c31637a-8666-11eb-2ea6-17b0d4507e10
@bind tt Slider(1:100)

# ╔═╡ 3bf19b4c-8668-11eb-1025-534324f5d70e
np = n_photons(l400P1mw, tt * 1.0s)

# ╔═╡ 28c0a218-8666-11eb-291a-71fd93c60da5
de = uconvert(mJ,delivered_energy(l400P1mw, tt * 1.0s))

# ╔═╡ 013420a8-9129-11eb-0d7f-ab8fb7c70786
md"#### FOV and photon density: no focus (FOV of 2 mm diameter)"

# ╔═╡ 7ff6fe7c-9129-11eb-1805-0dd052ec01a8
fovnf = FoV(2mm, 1μm)

# ╔═╡ bc92a90a-9129-11eb-136b-efe4bb6659ec
I_nf = photon_density(l400P1mw, fovnf)

# ╔═╡ f7e62bbc-9129-11eb-3366-fb91dfe68531
photon_density(405nm, 1mW, fovnf.a)

# ╔═╡ 5b1cb8ce-912a-11eb-3320-ffe3dfc7be67
@test photon_density(405nm, 1mW, fovnf.a) ≈ photon_density(l400P1mw, fovnf)

# ╔═╡ a80e6642-908c-11eb-3ac3-6737f11f4eff
md"### Number of photons in Silica powder (a first estimate):
- Standard concentration
- Assume a standard beam spot of radius 1 mm
- Use the cross section for FBI measured in solution"

# ╔═╡ 1d515f5e-908d-11eb-153b-eda702ae7084
mfbi = mfbisol250

# ╔═╡ b982d466-908d-11eb-3efb-5d537abfad23
mfbi.σ

# ╔═╡ ecdbdf64-9138-11eb-2f8c-8588c4177ba5
mfbiba = mfbibasol250

# ╔═╡ 0a6ef384-9139-11eb-161b-c3598134f47c
mfbiba.σ

# ╔═╡ 635f9390-9133-11eb-3eec-0b15c5c40ba7
md"- The nuber of photons **produced** per unit time (in Hz) is given by the fluorescence pero molecule (mfbi) corresponding to the non-focused spot (I_nf) times the number of molecules in the spot"

# ╔═╡ ebd52068-908d-11eb-2c2b-67c2571da02d
@test fluorescence_per_molecule(mfbi, I_nf) ≈ I_nf * uconvert(mm^2, mfbi.σ) * mfbi.Q

# ╔═╡ d409e8ba-9137-11eb-00f3-75db3de2e8eb
nsibimm2 = nof_molecules_area(sifbi, mm^-2)

# ╔═╡ 1641a013-e977-498b-b0d8-019d5e7ae8f9
nsifbibamm2 = nof_molecules_area(sifbiba, mm^-2)

# ╔═╡ 36ce81c2-9138-11eb-1353-81a30079575e
md" #### Number of photons produced in FBI/Ba powder, laser of 1 mm2 radius"

# ╔═╡ d773c5ac-9137-11eb-2832-edbc3b7f7328
γ_fbinf = fluorescence_per_molecule(mfbi, I_nf) * nof_molecules_area(sifbi, mm^-2) * fovnf.a

# ╔═╡ 1a52dd68-9138-11eb-1a6b-e705887868bc
md" - FBIBa2+ non-focused"

# ╔═╡ 4ee0b9d8-9138-11eb-113c-17b29b79628f
γ_fbibanf = fluorescence_per_molecule(mfbiba, I_nf) * nof_molecules_area(sifbiba, mm^-2) * fovnf.a

# ╔═╡ afe6ab7a-9133-11eb-1f9f-dd7c0b27e4d0
md"#### Number of photons detected:
 - non-focused setup
 - only geometrical acceptance taken into account
 - Laser light is collected by a lens of 25.4 mm located at 350 mm from target"

# ╔═╡ 55dc84a0-9139-11eb-1694-bba96f7bad9e
geometrical_acceptance(350.0, 25.4)

# ╔═╡ 67b2b186-9139-11eb-2f29-f9f36be38b6f
begin
	d_target_ccd = 350.0mm
	D_ccd        = 25.4mm
end

# ╔═╡ 437b6cee-9135-11eb-1e08-a70d334d271b
n_fbinf = γ_fbinf * geometrical_acceptance(d_target_ccd/mm, D_ccd/mm)

# ╔═╡ e112a0fe-9139-11eb-3e41-872c8d481b43
n_fbibanf = γ_fbibanf * geometrical_acceptance(d_target_ccd/mm, D_ccd/mm)

# ╔═╡ f96e9b12-9139-11eb-08cb-5319bca7a01d
md"## Setup
- Setup includes an objective (if present), a CCD to collect light and a collection of filters"

# ╔═╡ cfbd347a-866e-11eb-28f1-138f1ead212f
md"### Objective"

# ╔═╡ db362c80-866e-11eb-0fdd-0f7c8cfbfae1
struct Objective
	name              ::String
    numerical_aperture::Float64
    magnification     ::Float64
end

# ╔═╡ d3aaf1ce-8670-11eb-3f29-05b28699efcc
function transmission(objective::Objective)
	A = objective.numerical_aperture
	A <=1 ? (1 - sqrt(1 - A^2)) /2 : 0.5
end

# ╔═╡ b4badd12-8671-11eb-1f1b-69b36be8f7e4
begin
	naL = collect(0.1:0.05:1)
	obL = [Objective("Nikon-SLVD", na, 100) for na in naL]
	trL = [transmission(ob) for ob in obL]

	plot(naL, trL, leg=false, lw=2)

	xlabel!("numerical aperture")
	ylabel!("Transmission")
	title!("Tranmission as a function of NA")
end

# ╔═╡ 6188da98-867f-11eb-11b7-e9c11f87914f
md"### CCD"

# ╔═╡ 6873ba8a-867f-11eb-1a10-1fa51d56b5be
struct CCD
    name ::String
	eff  ::Any
	function CCD(name)
		function eff(ll::typeof(1.0nm))
			l = ustrip(ll)
			if l < 350.
				return 0.
			elseif l > 1000.
				return 0.
			else
				return e(l)
			end
		end
		l = 350.:50.:1000.
		ϵ = [0.3, 0.4,0.65,0.78,0.82,0.82,0.8,0.72,0.62,0.5,0.37,
                      0.24,0.12,0.07]
		e = CubicSplineInterpolation(l, ϵ)

		new(name, eff)
	end
end

# ╔═╡ 2bd60598-1716-4b43-8163-50b7f9471e92
function ccd(lmin::Float64=350.0, lmax::Float64=1000.0)
	function eff(l::Float64)
		if l < lmin || l > lmax
			return 0.
		else
			wl = 350.:50.:1000.
			ϵ = [0.3, 0.4,0.65,0.78,0.82,0.82,0.8,0.72,0.62,0.5,0.37,
			  0.24,0.12,0.07]
			e = CubicSplineInterpolation(wl, ϵ)
			return e(l)
		end
	end
	return eff
end

# ╔═╡ ebcccab3-6f16-436c-a2dc-c2909cab8ca2


# ╔═╡ 79af9066-8680-11eb-0d91-15759e995f6c
orca = ccd()

# ╔═╡ 2f259d3a-8688-11eb-1733-13978ae756ce
begin
	xL = collect(350.:5.:1000.)
	xEff =orca.(xL)
	plot(xL, xEff, leg=false, lw=2)

	xlabel!("λ (nm)")
	ylabel!("Efficiency")
	title!("Efficiency of CCD as a function of λ")
end

# ╔═╡ 56d72e3c-86fa-11eb-2928-798e913a5976
md"### Filters"

# ╔═╡ 11547512-86fe-11eb-2de2-19e00dbb9cf5
filter_names = Dict("FF01405-10" => "/Users/jj/JuliaProjects/LabFBI/data/filters/BandPass405.csv", "FB430-10" => "/Users/jj/JuliaProjects/LabFBI/data/filters/BandPass430.csv", "DoubleNotchFilter_405_532" => "/Users/jj/JuliaProjects/LabFBI/data/filters/DoubleNotchFilter_405_532.csv", "DMLP425"=> "/Users/jj/JuliaProjects/LabFBI/data/filters/LongPass425.csv","FELH0450"=>"/Users/jj/JuliaProjects/LabFBI/data/filters/LongPass450.csv")

# ╔═╡ 262f9496-8701-11eb-18fb-b5648cb8058a
function load_filter(filter_names::Dict, filtername::String)
	fdf = DataFrame(CSV.File(filter_names[filtername]))
	function scaledf(df)
		df[!, :T] = map(t -> t/100.0, df[:, :T])
		return df
	end
	function fixdf(df)
		df = df[!,3:4]
		rename!(df, [Symbol("λ"),  Symbol("T")])
		return scaledf(df)
	end

	if filtername == "FB430-10"
		fdf = fixdf(fdf)
		sort!(fdf, rev = false)
	elseif filtername == "DoubleNotchFilter_405_532"
		fdf = scaledf(fdf)
	elseif filtername == "DMLP425" || filtername == "FELH0450"
		fdf = fixdf(fdf)
	end

	fitp = interpolate((fdf.λ,), fdf.T, Gridded(Linear()))
	return fdf, fitp
end

# ╔═╡ e4855e4c-8703-11eb-0873-9f064139b2fa
f405df, f405fi = load_filter(filter_names, "FF01405-10");

# ╔═╡ 9b1aac58-8703-11eb-3aa8-c1c7b909f02f
function plot_filter(filtername, fdf, fint, wl)
	function labels()
		xlabel!("λ (nm)")
		ylabel!("T")
		title!(filtername)
	end

	ifbi = [fint(x*1.0) for x in wl]
	p1 = plot(fdf.λ, fdf.T,
		label = "Filter data",
    	shape = :circle,
    	color = :black,
    	markersize = 3, leg=:topright)
	p2 = plot!(p1,wl, ifbi,
		label = "Filter data interpolation",
    	color = :blue,
    	markersize = 3, leg=:topright)
	labels()
	p3 = plot(wl, ifbi,
		label = "Filter data interpolation",
    	color = :blue,
    	markersize = 3, leg=:topright)
	labels()

	plot(p2, p3, layout = (1, 2), legend = false, fmt = :png)
end

# ╔═╡ 86a72216-92f6-11eb-179a-d7a40a7c9cd2
plot_filter("FF01405-10",f405df, f405fi, 345:1:445)

# ╔═╡ b4ff7a5a-92f6-11eb-3aaa-4fb7c21f4b12
f430df, f430fi = load_filter(filter_names, "FB430-10");

# ╔═╡ f830f3d0-92f6-11eb-2ec6-dd2cc6983e70
plot_filter("FB430-10",f430df, f430fi, 410.0:1.0:450.0)

# ╔═╡ 132eb3da-930a-11eb-292b-d7abf1e8ab2b
fDN, fDNfi = load_filter(filter_names, "DoubleNotchFilter_405_532");

# ╔═╡ 6088f582-930a-11eb-09be-fd288fcd1be9
plot_filter("DoubleNotchFilter_405_532",fDN, fDNfi, 400:1:800)

# ╔═╡ 978fe7d4-930a-11eb-2792-15d899027d32
flp425df, flp425fi = load_filter(filter_names, "DMLP425");

# ╔═╡ ccd696ea-930a-11eb-211a-3da303ff1eea
plot_filter("DMLP425",flp425df, flp425fi, 400:1:800)

# ╔═╡ ed122ef6-930a-11eb-04ab-457fe2e1fcc9
flh450df, flh450fi = load_filter(filter_names, "FELH0450");

# ╔═╡ b0a2b0de-930b-11eb-3247-a3cab75a9de4
plot_filter("FELH0450",flh450df, flh450fi, 400:1:600)

# ╔═╡ f89d8589-80a9-4d06-80d6-3b9f1a1f78af
f405 = fpdf(f405fi, 375.0, 701.0, 1.0)

# ╔═╡ ca88a081-8a7a-4d83-a727-38b264266c18
f430bp = fpdf(f430fi, 375.0, 701.0, 1.0)

# ╔═╡ f0524443-8b01-4662-b241-341f829d457e
fdn = fpdf(fDNfi, 375.0, 701.0, 1.0)

# ╔═╡ 6283cc6e-38d3-4ecb-a1b8-c678413a3a64
f425lp = fpdf(flp425fi, 375.0, 701.0, 1.0)

# ╔═╡ 164ffc2c-52fc-482e-a8fd-e356f5865263
f450lp = fpdf(flh450fi, 375.0, 701.0, 1.0)

# ╔═╡ 80f3d9eb-2c40-4761-ab73-0119ff185cc0
function plot_filters(filterlist, filternames, wl)
	function labels(filtername)
		xlabel!("λ (nm)")
		ylabel!("T")
		title!(filtername)
	end
	PP = []
	for (i, ff) in enumerate(filterlist)
		ifbi = [ff(x*1.0) for x in wl]

		p = plot(wl, ifbi,
				 label = filternames[i],
    			 color = markercolors[i],
    			 lw=2, leg=:topright)
		labels(filternames[i])
		push!(PP,p)
	end

	return PP
end

# ╔═╡ 382aea80-665c-476c-9e68-63841b4472c2
PP = plot_filters([f405,f430bp,fdn,f425lp,f450lp],["Band 405", "Band 430", "Double Notch", "LP 425", "LP 450"],wl);

# ╔═╡ 37fc5ad1-f184-4384-bc2d-08beb88e1494
plot(PP[1], PP[2], layout = (1, 2), legend = false, fmt = :png)

# ╔═╡ 1ec7f29c-0a2f-4276-b1eb-2d03f08da2a5
plot(PP[3], PP[4], PP[5], layout = (1, 3), legend = false, fmt = :png)

# ╔═╡ b0b10062-93c9-11eb-1b87-0f76a1735fa1
md"## Laser transport
- Start with shape of FBI or FBIBa2+ (number of counts)
- Compute the shape after all filters
- Multiply by number of events
"

# ╔═╡ c126f2d9-0ccb-4b5c-9600-e14ad96e9591
filterSetBand430(λ) = map(x->f425lp(x)^2 * fdn(x) * f430bp(x) * orca(x), λ)

# ╔═╡ 7ccc9401-2fee-4d5e-a6ca-9a863c60d5a6
filterSetLP450(λ) = map(x->f425lp(x)^2 * fdn(x) * f450lp(x) * orca(x), λ)

# ╔═╡ 9bd87f70-24e6-4476-acf7-7647face1b3e
@test f425lp(430.0)^2 * fdn(430.0) * f430bp(430.0) * orca(430.) ≈ filterSetBand430(430.0)

# ╔═╡ a40e5ab4-ff6c-4cbc-a377-db60b76d1963
@test f425lp(460.0)^2 * fdn(460.0) * f450lp(460.0) * orca(460.0) ≈ filterSetLP450(460.0)

# ╔═╡ 84deb385-b6be-4b89-a43d-bd1b1730473a
function plot_filterset(fset, fname, wl)
	function labels()
		xlabel!("λ (nm)")
		ylabel!("T")
		title!(fname)
	end

	ifbi = [fset(x) for x in wl]
	p    = plot(wl, ifbi,
		label = fname,
    	color = :blue,
    	lw=2, leg=:false)
	labels()
end

# ╔═╡ e4abed36-f47b-4d77-9d10-e7455edd9325
plot_filterset(filterSetBand430, "filterSetBand430", 400.0:500.0)

# ╔═╡ 202f831f-8e54-4d17-9f27-aea1a07ab919
plot_filterset(filterSetLP450, "filterSetLP450", 400.0:700.0)

# ╔═╡ 241a2416-6bc7-4804-b6ae-874d9336a02a
md"### Filtered pdfs for FBI and FBIBa2+"

# ╔═╡ 85dddc12-0a69-4745-a895-2508ea3ffc5c
filteredFbi430BP(λ) = map(x->filterSetBand430(x) * fbipdf(x), λ)

# ╔═╡ e72d34ff-1a0a-4dd4-a088-49e671d01408
filteredFbiBa430BP(λ) = map(x->filterSetBand430(x) * fbibapdf(x), λ)

# ╔═╡ 60badadd-4df8-426a-b367-86e7bdc0c25f
filteredFbi450LP(λ) = map(x->filterSetLP450(x) * fbipdf(x), λ)

# ╔═╡ 6a5f6e6f-06cb-46cd-bab4-f3a0778af2b6
filteredFbiBa450LP(λ) = map(x->filterSetLP450(x) * fbibapdf(x), λ)

# ╔═╡ e10d9e49-fe4a-4304-80f8-87da9bfe6ac3
plot_filterset(filteredFbi430BP, "filteredFbi430BP", 400.0:700.0)

# ╔═╡ 03935fbd-a05b-47b3-afe8-43ce8fd532c6
plot_filterset(filteredFbiBa430BP, "filteredFbiBa430BP", 400.0:700.0)

# ╔═╡ b9453c4f-055a-4b8e-a14b-d7e116374dd0
plot_filterset(filteredFbi450LP, "filteredFbi450LP", 400.0:700.0)

# ╔═╡ c5e21983-0c21-44c1-b16f-3578e3c762ae
plot_filterset(filteredFbiBa450LP, "filteredFbiBa450LP", 400.0:700.0)

# ╔═╡ 95ce45d1-d000-46c4-9996-2aa6db2876cf
md"### Ratios and double ratio for filtered spectra
- Ratio of events in BP430/LP450 for both FBI and FBIBa
- Double ratio"

# ╔═╡ 12fdb81e-e1f2-44c0-a051-41b5be990342
fi430bp, eps5 = quadgk(filteredFbi430BP, 375.0,  700.0)

# ╔═╡ edfc6f80-f0dc-444c-bd74-ae4f226f9583
fi450lp, eps6 = quadgk(filteredFbi450LP, 375.0,  700.0)

# ╔═╡ cfef5324-ba28-412a-a927-f4f91480194b
rfbi = fi430bp / fi450lp

# ╔═╡ 0f1df1d8-2e18-4656-b44e-0bc6505a5bbd
fiba430bp, eps7 = quadgk(filteredFbiBa430BP, 375.0,  700.0)

# ╔═╡ a75af200-eede-48c9-9b3a-173e1b023bdc
fiba450lp, eps8 = quadgk(filteredFbiBa450LP, 375.0,  700.0)

# ╔═╡ 3325a9f2-004f-4b5b-84bc-bb2618a51a43
rfbiba = fiba430bp / fiba450lp

# ╔═╡ 6162ef9a-1337-4152-a688-61e3c32ed947
r2 = rfbiba / rfbi

# ╔═╡ 1a1b4a30-3966-40ab-8469-1509539defa2
md"### Ratios and double ratio for unfiltered specta"

# ╔═╡ e0b1c1d0-1d65-48c0-b09f-844bd27187fb
ufi430bp, ups1 = quadgk(fbipdf, 425.0,  435.0)

# ╔═╡ 2eb4afcd-a2c4-4962-83aa-2207640a861a
ufi450lp, ups2 = quadgk(fbipdf, 450.0,  700.0)

# ╔═╡ 5c3b7dd9-686d-42c2-8b41-8d559e52cf31
urfbi = ufi430bp / ufi450lp

# ╔═╡ abe6cb26-73c9-4e99-8672-9f8186548ff5
ufiba430bp, ups3 = quadgk(fbibapdf, 425.0,  435.0)

# ╔═╡ dea8bc12-e165-471f-8689-2d329a92b364
ufiba450lp, ups4 = quadgk(fbibapdf, 450.0,  700.0)

# ╔═╡ 16e956b2-8c1e-4a4f-90ac-13174fcab66c
urfbiba = ufiba430bp / ufiba450lp

# ╔═╡ 8e5ea441-fe39-4e32-ae80-986331f786d9
ur2 = urfbiba / urfbi

# ╔═╡ 1e325c2e-893a-4136-9302-45a61aa469df
md"### Number of events reaching CCD for FBI and FBIBa"

# ╔═╡ 5a48f61e-0de6-4b36-9f3e-6e539269ec03
md"`nccd(n::Float64, fpdf::Function)` takes the number of photons that would reach the CCD for a flat FBI spectrum and unit filter transfer function (n) and multiplies by fpdf which is the product of the FBI pdf and the filters transfer function"

# ╔═╡ 2197fe18-811a-4c2e-a1fc-daac14ff1fa0
function nccd(n::Float64, fpdf::Function, wmin::Float64, wmax::Float64)
	fi, eps = quadgk(fpdf, wmin, wmax)
	return n * fi
end

# ╔═╡ 14d074ec-1880-449f-abfa-eab2f431fd76
function nccd(n::Quantity, fpdf::Function, wmin::Float64, wmax::Float64)
	fi, eps = quadgk(fpdf, wmin, wmax)
	return n * fi
end

# ╔═╡ 0bb2ec9b-1741-451b-8dd2-31fde86098ad
md"`nccd_nf_430bp_fbi` is the number of photons that reach the CCD in a non-focused (nf) setup, passing the 430BP filter transfer function, for FBI"

# ╔═╡ 02e26731-8762-4680-8d15-3da475b3f3ce
n_fbinf

# ╔═╡ 5dc7bbb5-d23a-4a19-b353-b05f758406f0
nccd_nf_430bp_fbi = nccd(n_fbinf, filteredFbi430BP, 375.0, 700.0)

# ╔═╡ 144f6082-b4c7-4b9c-abf4-579516503966
@test n_fbinf * fi430bp ≈ nccd_nf_430bp_fbi

# ╔═╡ ac46c2cb-59e0-48a0-83c7-23bd2bf71f68
md"`nccd_nf_450lp_fbi` is the number of photons that reach the CCD in a non-focused (nf) setup, passing the 450LP filter transfer function, for FBI"

# ╔═╡ bbcba001-a680-439f-ade6-234613c97a93
nccd_nf_450lp_fbi = nccd(n_fbinf, filteredFbi450LP, 375.0, 700.0)

# ╔═╡ d62b1ebb-3e39-43ff-b05f-a357b98ff5f5
md"- The ratio between the number of events reaching the CCD must be the same than the ratio between sepectra"

# ╔═╡ 939c02c9-f488-4da2-8916-4a8089486f17
r_nccd_fbi = nccd_nf_430bp_fbi / nccd_nf_450lp_fbi

# ╔═╡ 8f3b2a4e-bbd9-45f6-8912-0686970c32ca
@test rfbi ≈ r_nccd_fbi

# ╔═╡ f48d5ff1-768c-48b1-9f05-da49c5457b17
md"`nccd_nf_430bp_fbiba` is the number of photons that reach the CCD in a non-focused (nf) setup, passing the 430BP filter transfer function, for FBIBa"

# ╔═╡ b69fbabd-9a06-45cf-900d-88038f3e804e
nccd_nf_430bp_fbiba = nccd(n_fbibanf, filteredFbiBa430BP, 375.0, 700.0)

# ╔═╡ 5526d5ec-a080-4832-b9e1-f6dcda1af954
md"`nccd_nf_450lp_fbiba` is the number of photons that reach the CCD in a non-focused (nf) setup, passing the 450LP filter transfer function, for FBIBa"

# ╔═╡ daeb3af4-ee9b-4d70-8a5c-9d37550bb199
nccd_nf_450lp_fbiba = nccd(n_fbibanf, filteredFbiBa450LP, 375.0, 700.0)

# ╔═╡ 985aedab-8704-4ff6-a925-3be3e3636e99
r_nccd_fbiba = nccd_nf_430bp_fbiba / nccd_nf_450lp_fbiba

# ╔═╡ 7178e04b-a890-44a7-8b72-ee56ccdef9af
@test rfbiba ≈ r_nccd_fbiba

# ╔═╡ Cell order:
# ╠═79cfd2fc-9046-11eb-2b13-1b877d57645d
# ╠═4f35dcfc-8644-11eb-1c86-75608b38bf5f
# ╠═2e3c7fda-904f-11eb-2988-25604b1caad0
# ╠═4ea80a97-44a6-4205-ae30-eaf239309313
# ╠═3207f446-8643-11eb-37ba-c9aec47fcb8f
# ╠═5115917a-8644-11eb-19fc-0528741ca75d
# ╠═8ba78a98-904f-11eb-1977-2750643d2f9f
# ╠═f2faa594-9316-11eb-2451-33e4b4d3a4f5
# ╠═5356699e-93ca-11eb-0afd-e7d0edc4231e
# ╠═fc8b79a2-8728-11eb-2da7-e3ffa3ceef08
# ╠═0f2f4c78-8729-11eb-2bab-27812ce8c47e
# ╠═5ee27d52-86fd-11eb-365e-9f2b1a095575
# ╠═621ec96c-86fd-11eb-1c41-379cc17180dc
# ╠═0b1cedc7-8ada-45a5-adef-fbae794dee3e
# ╟─70525ca4-8643-11eb-121e-196d80a60fcf
# ╠═e5e02216-86fb-11eb-1cf6-6f6962313a09
# ╟─3c01d7c8-8645-11eb-372d-4d66ae46ae72
# ╠═46b5465a-8645-11eb-0291-612455795518
# ╟─443e90e2-8664-11eb-19bc-51f8cfc282d5
# ╠═33ad8312-8661-11eb-1909-87eb9ed385bc
# ╠═e6e6be0a-9046-11eb-01ae-c572f38fca9e
# ╠═ebd713a6-9046-11eb-28f3-c14460ab4b0b
# ╠═f7a2b46a-9046-11eb-04b7-2971adb83918
# ╟─4c8ab940-8661-11eb-30a9-9301eb73ea33
# ╠═c1615b46-8663-11eb-21dc-5b746b6ae780
# ╠═a73a6004-8664-11eb-0cd2-6bbf19682226
# ╟─0f9a7e86-8665-11eb-0810-43b5d2227a40
# ╠═0374c94a-8665-11eb-28b9-8127833558fc
# ╟─341791a4-665f-4bc1-bc3d-35799fd4afec
# ╠═39d064d4-93cd-11eb-0530-7d4dd6b3c2d5
# ╠═414f254c-93cd-11eb-07ca-f323b0b1b36c
# ╟─754ab06a-9048-11eb-3953-bd88a367d21e
# ╠═4dc66f48-9048-11eb-2334-6d9115a0a741
# ╠═658252da-904a-11eb-15f0-036d2412b386
# ╠═58bfee7c-904a-11eb-07f8-cdafd513ac04
# ╟─a607fc3a-9048-11eb-2074-f188307b35a9
# ╠═537c6c58-9048-11eb-0176-97f2bb8e06f7
# ╠═6dd9d1b0-904a-11eb-1fec-9fe82cf1e5ee
# ╠═8dd2c6a4-904a-11eb-28f9-8f9ede41f1dd
# ╟─05ad4620-904b-11eb-3f68-13a12b8b5f27
# ╠═a04ab7ba-904a-11eb-21b9-8b09059f1ca5
# ╠═ee7f63de-904f-11eb-0e8f-1f8ec9d50811
# ╠═631833b8-0649-4f45-b1dd-e4a53091d321
# ╠═f8d3f48a-904f-11eb-2659-51f7508b434d
# ╠═9fdae9e2-6d2e-4901-8f48-14764b6800c2
# ╠═42d84ad0-904c-11eb-08d6-7b4b117bf2a6
# ╠═5f56cf78-5b2d-471b-b21f-f242b7d78175
# ╠═3f63bf34-9046-11eb-07a2-6b5f90371cb6
# ╠═a10b23d9-e366-4e38-80b6-a1532e30302e
# ╠═a90d8eea-904b-11eb-2c8d-b9fbc016aadb
# ╠═7de08c7a-7c93-45af-b46a-247cfc30c859
# ╠═1694d934-904c-11eb-0b79-5f88092d69ce
# ╠═ea7a474e-0d12-436a-b8a6-e0737507d921
# ╠═d6f6bfda-904c-11eb-31d1-3550103982f9
# ╠═15c9d8e0-904e-11eb-35e4-3d920e199e97
# ╟─bf51632d-4543-42b6-b9d1-f54c240412d4
# ╠═dece508c-9052-11eb-3d03-431affa529dc
# ╟─84e3a8e1-665b-4302-9f9d-a1fb82c13aa2
# ╠═a61dcc86-8643-11eb-274a-db5011271345
# ╟─52d85aa0-feaa-40aa-91a2-93b1b3534da6
# ╠═2287af6e-9128-11eb-3070-f7af71689367
# ╟─792cfd0e-89e0-478e-abc1-ee8fcbba79e8
# ╠═3adf7eba-865c-11eb-2cb2-01884cac9b61
# ╟─ae7ed83e-562a-4e13-a24e-a8e4c793a27b
# ╠═54c84cf6-908b-11eb-26e7-3d833153a5bc
# ╟─dc537c9f-3764-4055-a9bc-9eacdaf81fbb
# ╠═7b4761e4-908b-11eb-2840-293f5286d1d5
# ╟─1810e902-9ea4-41da-a7fc-60e03197bc7e
# ╠═6a28f556-8667-11eb-3614-d59ecab119d6
# ╟─89528344-0c5c-41d9-a1d5-62b315b90ef7
# ╠═811b7904-90bd-11eb-3182-4969edbde50e
# ╟─58f59d86-5a11-4dfd-952d-d3d5655999f7
# ╠═28d38d4c-9129-11eb-310e-1996bd60a771
# ╟─0630250e-8602-4556-8948-1f6858e0b587
# ╠═9f25699e-8665-11eb-13af-0be83c293556
# ╟─815a3ba6-4a61-4d96-9aca-d32318c7a99f
# ╠═9596d471-bf22-47f3-8f2d-577b8d875382
# ╟─cb4d9e52-9134-11eb-2557-d94ef2d4e541
# ╠═acba0a22-9134-11eb-3bd8-a7e99e8a710c
# ╟─45bdd766-9131-11eb-31dc-f9ee1ced0db9
# ╟─87199dd6-9047-11eb-312d-1fea3493dd46
# ╟─c3b5e2a0-90bc-11eb-243e-eb6908eeeb13
# ╠═304c3140-90bc-11eb-1bbd-632bdcce5ced
# ╠═6ff7fe82-90bc-11eb-3b33-271f008fe4b4
# ╟─538eccf0-904c-11eb-0a29-491192a77e31
# ╠═b7de26de-904b-11eb-223d-b9d5f33b1c5f
# ╠═e82978a2-904b-11eb-15af-8f1fedc428a5
# ╠═075148ae-904c-11eb-1969-1be58cdfb088
# ╠═7a93e858-904c-11eb-0a45-17e5bb3d405a
# ╠═91d4f7c8-904c-11eb-0486-23d11f6936df
# ╠═934c7fcc-904c-11eb-1339-d5a1e82ac071
# ╠═0773cad6-904d-11eb-1ba0-9b46e5d1046c
# ╟─c9017918-9050-11eb-2b1b-a9d8243e00a7
# ╟─f2d90250-9054-11eb-3008-01b74530d71e
# ╠═9b8f1f65-8476-4a2c-adc2-b1973776110e
# ╠═50c68b42-9052-11eb-0dd3-41e31d7a1947
# ╟─07fdc33c-9055-11eb-2f90-b9077aa5fe3d
# ╠═f63079f2-9053-11eb-0e1b-cd8542581fc2
# ╟─2f084308-9055-11eb-2ce2-bd48996c0f5e
# ╠═3dbd2044-9055-11eb-3ca4-75af1ab90b72
# ╟─66a99aaa-9055-11eb-1126-57af32a69b8e
# ╠═73885928-9055-11eb-20a0-41f136464b4f
# ╠═cadd0620-9315-11eb-228e-135bd066038d
# ╠═d65c78aa-9315-11eb-094c-01765ce6cc2f
# ╠═0eb3c064-9316-11eb-3c97-afba4053eb94
# ╠═322216ee-93cb-11eb-272c-b7f5e5d87ee5
# ╠═8ad489e2-9318-11eb-1531-7515b778e6b8
# ╠═9989f79a-9318-11eb-07a3-a9986579342d
# ╠═5361bf8e-9318-11eb-1acd-e9a10c2b30bf
# ╠═2ecbd182-93ca-11eb-08f6-3d0c1710f9f5
# ╠═eba136e2-93c9-11eb-21f8-8b36a5158624
# ╠═6deec574-93ca-11eb-0b5f-1b0d7ca04b61
# ╠═46411c0e-93cd-11eb-1039-8fe43a7679d5
# ╠═353b7b37-c037-49f9-b149-56541db6cc28
# ╠═566a984a-f0fc-482f-9b04-029cba608aef
# ╠═fee9eace-93cd-11eb-114c-c174cc6637ce
# ╠═4d17046e-ce2e-4863-8231-2ded5722ff44
# ╠═fc40d4e4-93ce-11eb-02ed-57cb27c0db82
# ╠═f06bb3a3-247c-4985-91cf-8ce79744945e
# ╠═2d0e3c95-fcd6-4d96-9ca0-bd8c8be01edb
# ╠═65dbf890-3d5b-4f1c-aa83-4d72556f9173
# ╠═c1f8e6d0-52b2-4ca7-8aaf-6ed25bbfce08
# ╠═0e678e6b-dcbd-4cd2-9be9-5f9decc2dcc8
# ╠═e4e97bd1-28c2-4340-ade7-50dab201fe1a
# ╠═b3287b7b-3a7f-4ef9-aaf1-4ed2754b7e72
# ╠═56e1f815-7337-499d-a4bc-ed65121984c5
# ╠═b55a5691-97b0-470e-8fd3-b143c26f724d
# ╠═054700b9-a1a6-4c04-94b7-037df2a890c4
# ╠═76a569c4-93cb-11eb-3bb4-55e0c9221e53
# ╟─9bf0b7e6-8643-11eb-142b-7fc904562cc8
# ╠═e9076d72-8648-11eb-0c02-9b00b0605e72
# ╠═2aef669a-864b-11eb-3290-359ff40cece4
# ╠═4bb9cb3e-864b-11eb-3b51-eb59b6da2c91
# ╠═1947e7a4-908b-11eb-0477-99f5152c6526
# ╟─2aff4f8c-908b-11eb-29af-93923b73c692
# ╠═6ce764b6-908b-11eb-3013-0f904fb0a0d1
# ╠═9ab0e624-908b-11eb-2d84-d93e1317e3bd
# ╠═3bf19b4c-8668-11eb-1025-534324f5d70e
# ╠═a773efa0-908b-11eb-12a4-7bce3bb0b3d6
# ╠═bb4fa686-908b-11eb-0aa0-27ab50bb0a4e
# ╠═1c31637a-8666-11eb-2ea6-17b0d4507e10
# ╠═28c0a218-8666-11eb-291a-71fd93c60da5
# ╟─013420a8-9129-11eb-0d7f-ab8fb7c70786
# ╠═7ff6fe7c-9129-11eb-1805-0dd052ec01a8
# ╠═bc92a90a-9129-11eb-136b-efe4bb6659ec
# ╠═f7e62bbc-9129-11eb-3366-fb91dfe68531
# ╠═5b1cb8ce-912a-11eb-3320-ffe3dfc7be67
# ╟─a80e6642-908c-11eb-3ac3-6737f11f4eff
# ╠═1d515f5e-908d-11eb-153b-eda702ae7084
# ╠═b982d466-908d-11eb-3efb-5d537abfad23
# ╠═ecdbdf64-9138-11eb-2f8c-8588c4177ba5
# ╠═0a6ef384-9139-11eb-161b-c3598134f47c
# ╟─635f9390-9133-11eb-3eec-0b15c5c40ba7
# ╠═ebd52068-908d-11eb-2c2b-67c2571da02d
# ╠═d409e8ba-9137-11eb-00f3-75db3de2e8eb
# ╠═1641a013-e977-498b-b0d8-019d5e7ae8f9
# ╟─36ce81c2-9138-11eb-1353-81a30079575e
# ╠═d773c5ac-9137-11eb-2832-edbc3b7f7328
# ╟─1a52dd68-9138-11eb-1a6b-e705887868bc
# ╠═4ee0b9d8-9138-11eb-113c-17b29b79628f
# ╠═afe6ab7a-9133-11eb-1f9f-dd7c0b27e4d0
# ╠═55dc84a0-9139-11eb-1694-bba96f7bad9e
# ╠═67b2b186-9139-11eb-2f29-f9f36be38b6f
# ╠═437b6cee-9135-11eb-1e08-a70d334d271b
# ╠═e112a0fe-9139-11eb-3e41-872c8d481b43
# ╠═f96e9b12-9139-11eb-08cb-5319bca7a01d
# ╠═cfbd347a-866e-11eb-28f1-138f1ead212f
# ╠═db362c80-866e-11eb-0fdd-0f7c8cfbfae1
# ╠═d3aaf1ce-8670-11eb-3f29-05b28699efcc
# ╠═b4badd12-8671-11eb-1f1b-69b36be8f7e4
# ╠═6188da98-867f-11eb-11b7-e9c11f87914f
# ╠═6873ba8a-867f-11eb-1a10-1fa51d56b5be
# ╠═2bd60598-1716-4b43-8163-50b7f9471e92
# ╠═ebcccab3-6f16-436c-a2dc-c2909cab8ca2
# ╠═79af9066-8680-11eb-0d91-15759e995f6c
# ╠═2f259d3a-8688-11eb-1733-13978ae756ce
# ╟─56d72e3c-86fa-11eb-2928-798e913a5976
# ╠═11547512-86fe-11eb-2de2-19e00dbb9cf5
# ╠═262f9496-8701-11eb-18fb-b5648cb8058a
# ╠═e4855e4c-8703-11eb-0873-9f064139b2fa
# ╠═9b1aac58-8703-11eb-3aa8-c1c7b909f02f
# ╠═86a72216-92f6-11eb-179a-d7a40a7c9cd2
# ╠═b4ff7a5a-92f6-11eb-3aaa-4fb7c21f4b12
# ╠═f830f3d0-92f6-11eb-2ec6-dd2cc6983e70
# ╠═132eb3da-930a-11eb-292b-d7abf1e8ab2b
# ╠═6088f582-930a-11eb-09be-fd288fcd1be9
# ╠═978fe7d4-930a-11eb-2792-15d899027d32
# ╠═ccd696ea-930a-11eb-211a-3da303ff1eea
# ╠═ed122ef6-930a-11eb-04ab-457fe2e1fcc9
# ╠═b0a2b0de-930b-11eb-3247-a3cab75a9de4
# ╠═f89d8589-80a9-4d06-80d6-3b9f1a1f78af
# ╠═ca88a081-8a7a-4d83-a727-38b264266c18
# ╠═f0524443-8b01-4662-b241-341f829d457e
# ╠═6283cc6e-38d3-4ecb-a1b8-c678413a3a64
# ╠═164ffc2c-52fc-482e-a8fd-e356f5865263
# ╠═80f3d9eb-2c40-4761-ab73-0119ff185cc0
# ╠═382aea80-665c-476c-9e68-63841b4472c2
# ╠═37fc5ad1-f184-4384-bc2d-08beb88e1494
# ╠═1ec7f29c-0a2f-4276-b1eb-2d03f08da2a5
# ╠═b0b10062-93c9-11eb-1b87-0f76a1735fa1
# ╠═c126f2d9-0ccb-4b5c-9600-e14ad96e9591
# ╠═7ccc9401-2fee-4d5e-a6ca-9a863c60d5a6
# ╠═9bd87f70-24e6-4476-acf7-7647face1b3e
# ╠═a40e5ab4-ff6c-4cbc-a377-db60b76d1963
# ╠═84deb385-b6be-4b89-a43d-bd1b1730473a
# ╠═e4abed36-f47b-4d77-9d10-e7455edd9325
# ╠═202f831f-8e54-4d17-9f27-aea1a07ab919
# ╠═241a2416-6bc7-4804-b6ae-874d9336a02a
# ╠═85dddc12-0a69-4745-a895-2508ea3ffc5c
# ╠═e72d34ff-1a0a-4dd4-a088-49e671d01408
# ╠═60badadd-4df8-426a-b367-86e7bdc0c25f
# ╠═6a5f6e6f-06cb-46cd-bab4-f3a0778af2b6
# ╠═e10d9e49-fe4a-4304-80f8-87da9bfe6ac3
# ╠═03935fbd-a05b-47b3-afe8-43ce8fd532c6
# ╠═b9453c4f-055a-4b8e-a14b-d7e116374dd0
# ╠═c5e21983-0c21-44c1-b16f-3578e3c762ae
# ╠═95ce45d1-d000-46c4-9996-2aa6db2876cf
# ╠═12fdb81e-e1f2-44c0-a051-41b5be990342
# ╠═edfc6f80-f0dc-444c-bd74-ae4f226f9583
# ╠═cfef5324-ba28-412a-a927-f4f91480194b
# ╠═0f1df1d8-2e18-4656-b44e-0bc6505a5bbd
# ╠═a75af200-eede-48c9-9b3a-173e1b023bdc
# ╠═3325a9f2-004f-4b5b-84bc-bb2618a51a43
# ╠═6162ef9a-1337-4152-a688-61e3c32ed947
# ╠═1a1b4a30-3966-40ab-8469-1509539defa2
# ╠═e0b1c1d0-1d65-48c0-b09f-844bd27187fb
# ╠═2eb4afcd-a2c4-4962-83aa-2207640a861a
# ╠═5c3b7dd9-686d-42c2-8b41-8d559e52cf31
# ╠═abe6cb26-73c9-4e99-8672-9f8186548ff5
# ╠═dea8bc12-e165-471f-8689-2d329a92b364
# ╠═16e956b2-8c1e-4a4f-90ac-13174fcab66c
# ╠═8e5ea441-fe39-4e32-ae80-986331f786d9
# ╠═1e325c2e-893a-4136-9302-45a61aa469df
# ╠═5a48f61e-0de6-4b36-9f3e-6e539269ec03
# ╠═2197fe18-811a-4c2e-a1fc-daac14ff1fa0
# ╠═14d074ec-1880-449f-abfa-eab2f431fd76
# ╟─0bb2ec9b-1741-451b-8dd2-31fde86098ad
# ╠═02e26731-8762-4680-8d15-3da475b3f3ce
# ╠═5dc7bbb5-d23a-4a19-b353-b05f758406f0
# ╠═144f6082-b4c7-4b9c-abf4-579516503966
# ╟─ac46c2cb-59e0-48a0-83c7-23bd2bf71f68
# ╠═bbcba001-a680-439f-ade6-234613c97a93
# ╠═d62b1ebb-3e39-43ff-b05f-a357b98ff5f5
# ╠═939c02c9-f488-4da2-8916-4a8089486f17
# ╠═8f3b2a4e-bbd9-45f6-8912-0686970c32ca
# ╟─f48d5ff1-768c-48b1-9f05-da49c5457b17
# ╠═b69fbabd-9a06-45cf-900d-88038f3e804e
# ╟─5526d5ec-a080-4832-b9e1-f6dcda1af954
# ╠═daeb3af4-ee9b-4d70-8a5c-9d37550bb199
# ╠═985aedab-8704-4ff6-a925-3be3e3636e99
# ╠═7178e04b-a890-44a7-8b72-ee56ccdef9af

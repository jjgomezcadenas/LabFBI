### A Pluto.jl notebook ###
# v0.14.1

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
	Pkg.add.(["Unitful", "UnitfulEquivalences", "Plots","Interpolations",
			"PhysicalConstants","Images", "ImageIO", "ImageView","Peaks","Formatting",             "StrLiterals", "StrFormat", "Format"])
end

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
	using PhysicalConstants.CODATA2018
	using StatsPlots
	using QuadGK
	using Peaks
	using Formatting
	using Printf
	using StrLiterals
	using StrFormat
	using Format
end

# ╔═╡ fc8b79a2-8728-11eb-2da7-e3ffa3ceef08
using DrWatson

# ╔═╡ 79cfd2fc-9046-11eb-2b13-1b877d57645d
md"# TOPATU setup

- Studies the TOPATU setup, described below
"

# ╔═╡ 68e738e7-88bd-41c2-89e4-594f07d64ddc
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

# ╔═╡ 0f2f4c78-8729-11eb-2bab-27812ce8c47e
@quickactivate "LabFBI"

# ╔═╡ 5ee27d52-86fd-11eb-365e-9f2b1a095575
;pwd()

# ╔═╡ 621ec96c-86fd-11eb-1c41-379cc17180dc
datadir()

# ╔═╡ 9b853f27-4288-42a6-8f12-ca004e1773b7
srcdir()

# ╔═╡ cf89b973-7b0f-483d-8ce9-ba426f1df2a6
fbi = ingredients(srcdir("fbi.jl"))

# ╔═╡ a574c593-c72a-486b-bac4-91827581ee2e
lfbi = ingredients(srcdir("labFbi.jl"))

# ╔═╡ a1821b92-8b46-494b-a063-ae8146caf294
filters = ingredients(srcdir("filters.jl"))

# ╔═╡ 0b1cedc7-8ada-45a5-adef-fbae794dee3e
markercolors = [:green :orange :black :purple :red  :yellow :brown :white]

# ╔═╡ 70525ca4-8643-11eb-121e-196d80a60fcf
md"## TOPATU Setup"

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
    A, N, mol, mmol, V, L, mL, μL, M

# ╔═╡ 631833b8-0649-4f45-b1dd-e4a53091d321
md"#### Avogadro Number"

# ╔═╡ f8d3f48a-904f-11eb-2659-51f7508b434d
import PhysicalConstants.CODATA2018: N_A

# ╔═╡ 9fdae9e2-6d2e-4901-8f48-14764b6800c2
N_A

# ╔═╡ cb4d9e52-9134-11eb-2557-d94ef2d4e541
md"function `geometrical_acceptance(d,D)` computes the fraction of photons that make it through a lens or hole of diameter D located at a distance d from the emission point (no focusing)"

# ╔═╡ 45bdd766-9131-11eb-31dc-f9ee1ced0db9
md"# Notebook"

# ╔═╡ 9bf0b7e6-8643-11eb-142b-7fc904562cc8
md"### Laser

- Topatu has a Blue laser, with wavelength 405 nm and power up to 100 mW
"

# ╔═╡ e9076d72-8648-11eb-0c02-9b00b0605e72
l400P1mw = lfbi.Laser(405.0nm, 1.0mW);

# ╔═╡ 2aef669a-864b-11eb-3290-359ff40cece4
@test l400P1mw.λ == 405.0nm

# ╔═╡ 4bb9cb3e-864b-11eb-3b51-eb59b6da2c91
@test l400P1mw.P == 1.0mW

# ╔═╡ 1947e7a4-908b-11eb-0477-99f5152c6526
begin
	Λ =collect(170:1:800) * nm
	PE = lfbi.photon_energy.(Λ)
	plot(Λ./nm, PE./eV, leg=false,lw=2)
	xlabel!("λ (nm)")
	ylabel!("Photon energy (eV)")
	title!("Photon energy as a function of λ")
end

# ╔═╡ 2aff4f8c-908b-11eb-29af-93923b73c692
md"#### Number of photons per unit time as a function of power"

# ╔═╡ 6ce764b6-908b-11eb-3013-0f904fb0a0d1
lfbi.n_photons(l400P1mw)

# ╔═╡ 9ab0e624-908b-11eb-2d84-d93e1317e3bd
@test lfbi.n_photons(l400P1mw)≈lfbi.n_photons(405nm, 1mW)

# ╔═╡ a773efa0-908b-11eb-12a4-7bce3bb0b3d6
begin
	P =collect(1:1:10^3) * mW
	NP = lfbi.n_photons.((405nm,), P )
	plot(P./mW, NP./Hz, leg=false,lw=2)
	xlabel!("P (mW)")
	ylabel!("number of photons (Hz)")
	title!("Number of photons as a function of P")
end

# ╔═╡ bb4fa686-908b-11eb-0aa0-27ab50bb0a4e
md"#### delivered energy"

# ╔═╡ 1c31637a-8666-11eb-2ea6-17b0d4507e10
@bind tt Slider(1:100)

# ╔═╡ 28c0a218-8666-11eb-291a-71fd93c60da5
de = uconvert(mJ, lfbi.delivered_energy(l400P1mw, tt * 1.0s))

# ╔═╡ 013420a8-9129-11eb-0d7f-ab8fb7c70786
md"#### FOV and photon density: no focus (FOV of 2 mm diameter)"

# ╔═╡ 7ff6fe7c-9129-11eb-1805-0dd052ec01a8
fovnf = lfbi.FoV(2mm, 1μm);

# ╔═╡ bc92a90a-9129-11eb-136b-efe4bb6659ec
Inf = lfbi.photon_density(l400P1mw, fovnf)

# ╔═╡ f7e62bbc-9129-11eb-3366-fb91dfe68531
@test Inf ≈ lfbi.photon_density(405nm, 1mW, fovnf.a)

# ╔═╡ a80e6642-908c-11eb-3ac3-6737f11f4eff
md"#### Number of photons produced in Silica powder (a first estimate):
- Standard concentration
- Assume a standard beam spot of radius 1 mm (2 mm FoV)
- Use the cross section for FBI measured in solution"

# ╔═╡ 1d515f5e-908d-11eb-153b-eda702ae7084
mfbi = lfbi.Molecule("FBI-solution", 250nm, 489nm, 11260/(M*cm), 0.67);

# ╔═╡ b982d466-908d-11eb-3efb-5d537abfad23
mfbi.σ

# ╔═╡ fbb47644-74e8-423f-b867-5d51c31bceac
@test lfbi.fluorescence_per_molecule(mfbi, Inf) ≈uconvert(Hz, mfbi.σ * mfbi.Q * Inf)

# ╔═╡ 4ed47097-0d51-4639-a19f-38bbaa6c508d
γm = f"\%7.1f(lfbi.fluorescence_per_molecule(mfbi, Inf)/Hz)"

# ╔═╡ 5493a8d1-a31a-4338-91cc-688a8df4b048
md" - The number of photons per molecule is : $\gamma:m = \sigma Q I$ where I is the photon density. We find γm = $γm Hz"

# ╔═╡ 991eb30d-c5eb-4ed6-978e-1556961df89a
fbiSipRef = lfbi.Powder("FBI-Standard", 2.27e-5mmol/mg, 30mg/cm^2);


# ╔═╡ dd61fcbc-f60d-4426-b90a-81023b135eae
@test lfbi.nof_molecules_area(fbiSipRef, fovnf)≈  lfbi.nof_molecules_area(fbiSipRef, fovnf.a)

# ╔═╡ 0973ae7d-cb35-4e0d-95cb-fe0258114e27
nmL = f"\%7.2g(lfbi.nof_molecules_area(fbiSipRef, fovnf))"

# ╔═╡ 6b16643a-e656-453c-9fce-84c3d3d76963
md"- In the case of a Power we know the concentration of molecules per area (the assumption here is that of a thin film). For a FoV of 2 mm in the Reference FBI poweder concentration we find that the number of molecules is nmL = $nmL"

# ╔═╡ c620fc65-fa0e-47fe-a8dc-b04e2d7064b7
nγfov =  f"\%7.2g(lfbi.fluorescence_per_molecule(mfbi, Inf) * lfbi.nof_molecules_area(fbiSipRef, fovnf)/Hz)"

# ╔═╡ 506fde2b-935a-414d-a51e-afb75959e03e
md"- The number of photons produced in the beam FoV equals the number of molecules in the FoV times the fluorescence per molecule Thus, the rate of photons produced in the FOV is nγ_fov = $nγfov Hz "

# ╔═╡ afe6ab7a-9133-11eb-1f9f-dd7c0b27e4d0
md"#### Number of photons reaching CCDs (in the absence of filters):
 - non-focused setup
 - only geometrical acceptance taken into account"

# ╔═╡ 67b2b186-9139-11eb-2f29-f9f36be38b6f
begin
	d_target_ccd = 350.0mm
	D_ccd        = 25.4mm
end

# ╔═╡ 55dc84a0-9139-11eb-1694-bba96f7bad9e
ga = f"\%7.2g(lfbi.geometrical_acceptance(350.0, 25.4))"

# ╔═╡ 3d52a298-7ee2-4237-9f1f-39448d624867
md"- The geometrical acceptance for the unfocused setup (Laser light is collected by a lens of $(D_ccd/mm) mm located at $(d_target_ccd/mm) mm from the target) is ga = $ga "

# ╔═╡ 3756d51b-bdd0-4128-b345-b9cae6e48aef
nγga =  f"\%7.2g(lfbi.fluorescence_per_molecule(mfbi, Inf) * lfbi.nof_molecules_area(fbiSipRef, fovnf) * lfbi.geometrical_acceptance(d_target_ccd/mm, D_ccd/mm)/Hz)"

# ╔═╡ 88b4313c-c645-40b5-8894-8b0eb2ced5cf
md" - The number o photons reaching the CCD in the absence of filter (only geometrical acceptance) is nγga = $nγga"

# ╔═╡ f96e9b12-9139-11eb-08cb-5319bca7a01d
md"## Setup
- Setup includes an objective (if present), a CCD to collect light and a collection of filters"

# ╔═╡ b4badd12-8671-11eb-1f1b-69b36be8f7e4
begin
	naL = collect(0.1:0.05:1)
	obL = [lfbi.Objective("Nikon-SLVD", na, 100) for na in naL]
	trL = [lfbi.transmission(ob) for ob in obL]

	plot(naL, trL, leg=false, lw=2)

	xlabel!("numerical aperture")
	ylabel!("Transmission")
	title!("Tranmission as a function of NA")
end

# ╔═╡ 79af9066-8680-11eb-0d91-15759e995f6c
orca = lfbi.ccd();

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

# ╔═╡ 920bf9b9-7422-46cf-a772-1ee05e814474
begin
	λmin = 350.0nm
	λmax = 800.0nm
end

# ╔═╡ 52fb0e6a-5dba-4c09-be76-ac9665c09376
#fnames = ["BandPass405", "LongPass425", "BandPass430", "LongPass450", "DoubleNotch405_522"]

# ╔═╡ 58c093c7-5146-47c6-a01c-fb74a18185ef
fnames = filters.fnames

# ╔═╡ d55da805-3a0d-40f7-a514-aa8dd775cac5
filter_names = filters.filter_names

# ╔═╡ e4855e4c-8703-11eb-0873-9f064139b2fa
begin
	fbp405df, fbp405 = filters.load_filter(filter_names, "BandPass405")
	fbp405 = fbi.gfpdf(fbp405, λmin/nm, λmax/nm)
end

# ╔═╡ 86a72216-92f6-11eb-179a-d7a40a7c9cd2
plot_filter("BandPass405",fbp405df, fbp405, 345:1:445)

# ╔═╡ b4ff7a5a-92f6-11eb-3aaa-4fb7c21f4b12
begin
	fbp430df, fbp430 = filters.load_filter(filter_names, "BandPass430")
	fbp430 = fbi.gfpdf(fbp430, λmin/nm, λmax/nm);
end

# ╔═╡ f830f3d0-92f6-11eb-2ec6-dd2cc6983e70
plot_filter("BandPass430",fbp430df, fbp430, 410.0:1.0:450.0)

# ╔═╡ 132eb3da-930a-11eb-292b-d7abf1e8ab2b
begin
	fdn405_522df, fdn405_522 = filters.load_filter(filter_names, "DoubleNotch405_522")
	fdn405_522 = fbi.gfpdf(fdn405_522, λmin/nm, λmax/nm);
end

# ╔═╡ 6088f582-930a-11eb-09be-fd288fcd1be9
plot_filter("DoubleNotch405_522",fdn405_522df, fdn405_522, 400:1:800)

# ╔═╡ 978fe7d4-930a-11eb-2792-15d899027d32
begin
	flp425df, flp425 = filters.load_filter(filter_names, "LongPass425")
	flp425 = fbi.gfpdf(flp425, λmin/nm, λmax/nm)
end

# ╔═╡ ccd696ea-930a-11eb-211a-3da303ff1eea
plot_filter("LongPass425",flp425df, flp425, 400:1:800)

# ╔═╡ 56ac8fe0-cd1b-477d-93e0-0e0798a22f25
begin
	flp450df, flp450 = filters.load_filter(filter_names, "LongPass450")
	flp450 = fbi.gfpdf(flp425, λmin/nm, λmax/nm)
end

# ╔═╡ b0a2b0de-930b-11eb-3247-a3cab75a9de4
plot_filter("LongPass450",flp450df, flp450, 400:1:600)

# ╔═╡ 382aea80-665c-476c-9e68-63841b4472c2
PP = plot_filters([fbp405,flp425,fbp430,flp450,fdn405_522],fnames,350:800);

# ╔═╡ 37fc5ad1-f184-4384-bc2d-08beb88e1494
plot(PP[1], PP[2], layout = (1, 2), legend = false, fmt = :png)

# ╔═╡ 1ec7f29c-0a2f-4276-b1eb-2d03f08da2a5
plot(PP[3], PP[4], PP[5], layout = (1, 3), legend = false, fmt = :png)

# ╔═╡ b0b10062-93c9-11eb-1b87-0f76a1735fa1
md"## Laser transport

The TOPATU setup has two configurations: 
- Band-430: Select light in the band 430 +- 10. This is done combining:
	- The 405 nm BP filter (not needed for this calculation, since the beam is filtered before reaching the optical path, the filter contributes, together with other factors to a normalisation factor regarding power that needs to be computed)
	-The 425 nm LP filter, which enters twice in the optical path and eliminates high frequencies (below 425 nm).
	-The double noth which selects a region between 405 and 522 nm 
	-The 430 nm BP filter, which selects the signal (FbIBa2+ peak) region

- LP-450: Select light in the long pass Band 450 nm. This is done using the same combination than before except the last filter which is a 450 nm LP filter. This selects primarily the region of background (FBI).

- In addition one needs to add the efficiency of the CCD to have the complete transport function for both cases
"

# ╔═╡ c126f2d9-0ccb-4b5c-9600-e14ad96e9591
filterSetBand430(λ) = map(x->flp425(x)^2 * fdn405_522(x) * fbp430(x) * orca(x), λ)

# ╔═╡ 7ccc9401-2fee-4d5e-a6ca-9a863c60d5a6
filterSetLP450(λ) = map(x->flp425(x)^2 * fdn405_522(x) * flp450(x) * orca(x), λ)

# ╔═╡ 9bd87f70-24e6-4476-acf7-7647face1b3e
@test flp425(430.0)^2 * fdn405_522(430.0) * fbp430(430.0) * orca(430.) ≈ filterSetBand430(430.0)

# ╔═╡ a40e5ab4-ff6c-4cbc-a377-db60b76d1963
@test flp425(460.0)^2 * fdn405_522(460.0) * flp450(460.0) * orca(460.0) ≈ filterSetLP450(460.0)

# ╔═╡ e4abed36-f47b-4d77-9d10-e7455edd9325
plot_filterset(filterSetBand430, "filterSetBand430", 400.0:500.0)

# ╔═╡ 202f831f-8e54-4d17-9f27-aea1a07ab919
plot_filterset(filterSetLP450, "filterSetLP450", 400.0:700.0)

# ╔═╡ 241a2416-6bc7-4804-b6ae-874d9336a02a
md"### Filtered pdfs for FBI and FBIBa2+"

# ╔═╡ 23fabeb9-a947-4d44-9da2-a929cf9db6a4
md"- Get the spectrum from FBI and FBIBa powder"

# ╔═╡ 263448a0-a803-46b9-b7b5-d162e578f91f
fbidfname = datadir("fluorimeter/EDI044_FBI_round4_em325_375_405.csv")

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
# ╟─79cfd2fc-9046-11eb-2b13-1b877d57645d
# ╠═4f35dcfc-8644-11eb-1c86-75608b38bf5f
# ╠═3207f446-8643-11eb-37ba-c9aec47fcb8f
# ╠═5115917a-8644-11eb-19fc-0528741ca75d
# ╠═fc8b79a2-8728-11eb-2da7-e3ffa3ceef08
# ╠═68e738e7-88bd-41c2-89e4-594f07d64ddc
# ╠═0f2f4c78-8729-11eb-2bab-27812ce8c47e
# ╠═5ee27d52-86fd-11eb-365e-9f2b1a095575
# ╠═621ec96c-86fd-11eb-1c41-379cc17180dc
# ╠═9b853f27-4288-42a6-8f12-ca004e1773b7
# ╠═cf89b973-7b0f-483d-8ce9-ba426f1df2a6
# ╠═a574c593-c72a-486b-bac4-91827581ee2e
# ╠═a1821b92-8b46-494b-a063-ae8146caf294
# ╠═0b1cedc7-8ada-45a5-adef-fbae794dee3e
# ╟─70525ca4-8643-11eb-121e-196d80a60fcf
# ╠═e5e02216-86fb-11eb-1cf6-6f6962313a09
# ╟─3c01d7c8-8645-11eb-372d-4d66ae46ae72
# ╠═46b5465a-8645-11eb-0291-612455795518
# ╠═631833b8-0649-4f45-b1dd-e4a53091d321
# ╠═f8d3f48a-904f-11eb-2659-51f7508b434d
# ╠═9fdae9e2-6d2e-4901-8f48-14764b6800c2
# ╟─cb4d9e52-9134-11eb-2557-d94ef2d4e541
# ╟─45bdd766-9131-11eb-31dc-f9ee1ced0db9
# ╠═9bf0b7e6-8643-11eb-142b-7fc904562cc8
# ╠═e9076d72-8648-11eb-0c02-9b00b0605e72
# ╠═2aef669a-864b-11eb-3290-359ff40cece4
# ╠═4bb9cb3e-864b-11eb-3b51-eb59b6da2c91
# ╠═1947e7a4-908b-11eb-0477-99f5152c6526
# ╟─2aff4f8c-908b-11eb-29af-93923b73c692
# ╠═6ce764b6-908b-11eb-3013-0f904fb0a0d1
# ╠═9ab0e624-908b-11eb-2d84-d93e1317e3bd
# ╠═a773efa0-908b-11eb-12a4-7bce3bb0b3d6
# ╠═bb4fa686-908b-11eb-0aa0-27ab50bb0a4e
# ╠═1c31637a-8666-11eb-2ea6-17b0d4507e10
# ╠═28c0a218-8666-11eb-291a-71fd93c60da5
# ╟─013420a8-9129-11eb-0d7f-ab8fb7c70786
# ╠═7ff6fe7c-9129-11eb-1805-0dd052ec01a8
# ╠═bc92a90a-9129-11eb-136b-efe4bb6659ec
# ╠═f7e62bbc-9129-11eb-3366-fb91dfe68531
# ╠═a80e6642-908c-11eb-3ac3-6737f11f4eff
# ╠═1d515f5e-908d-11eb-153b-eda702ae7084
# ╠═b982d466-908d-11eb-3efb-5d537abfad23
# ╠═fbb47644-74e8-423f-b867-5d51c31bceac
# ╠═4ed47097-0d51-4639-a19f-38bbaa6c508d
# ╠═5493a8d1-a31a-4338-91cc-688a8df4b048
# ╠═991eb30d-c5eb-4ed6-978e-1556961df89a
# ╠═dd61fcbc-f60d-4426-b90a-81023b135eae
# ╠═0973ae7d-cb35-4e0d-95cb-fe0258114e27
# ╠═6b16643a-e656-453c-9fce-84c3d3d76963
# ╠═c620fc65-fa0e-47fe-a8dc-b04e2d7064b7
# ╠═506fde2b-935a-414d-a51e-afb75959e03e
# ╠═afe6ab7a-9133-11eb-1f9f-dd7c0b27e4d0
# ╠═67b2b186-9139-11eb-2f29-f9f36be38b6f
# ╠═55dc84a0-9139-11eb-1694-bba96f7bad9e
# ╠═3d52a298-7ee2-4237-9f1f-39448d624867
# ╠═3756d51b-bdd0-4128-b345-b9cae6e48aef
# ╠═88b4313c-c645-40b5-8894-8b0eb2ced5cf
# ╠═f96e9b12-9139-11eb-08cb-5319bca7a01d
# ╠═b4badd12-8671-11eb-1f1b-69b36be8f7e4
# ╠═79af9066-8680-11eb-0d91-15759e995f6c
# ╠═2f259d3a-8688-11eb-1733-13978ae756ce
# ╟─56d72e3c-86fa-11eb-2928-798e913a5976
# ╠═9b1aac58-8703-11eb-3aa8-c1c7b909f02f
# ╠═80f3d9eb-2c40-4761-ab73-0119ff185cc0
# ╠═84deb385-b6be-4b89-a43d-bd1b1730473a
# ╠═920bf9b9-7422-46cf-a772-1ee05e814474
# ╠═52fb0e6a-5dba-4c09-be76-ac9665c09376
# ╠═58c093c7-5146-47c6-a01c-fb74a18185ef
# ╠═d55da805-3a0d-40f7-a514-aa8dd775cac5
# ╠═e4855e4c-8703-11eb-0873-9f064139b2fa
# ╠═86a72216-92f6-11eb-179a-d7a40a7c9cd2
# ╠═b4ff7a5a-92f6-11eb-3aaa-4fb7c21f4b12
# ╠═f830f3d0-92f6-11eb-2ec6-dd2cc6983e70
# ╠═132eb3da-930a-11eb-292b-d7abf1e8ab2b
# ╠═6088f582-930a-11eb-09be-fd288fcd1be9
# ╠═978fe7d4-930a-11eb-2792-15d899027d32
# ╠═ccd696ea-930a-11eb-211a-3da303ff1eea
# ╠═56ac8fe0-cd1b-477d-93e0-0e0798a22f25
# ╠═b0a2b0de-930b-11eb-3247-a3cab75a9de4
# ╠═382aea80-665c-476c-9e68-63841b4472c2
# ╠═37fc5ad1-f184-4384-bc2d-08beb88e1494
# ╠═1ec7f29c-0a2f-4276-b1eb-2d03f08da2a5
# ╠═b0b10062-93c9-11eb-1b87-0f76a1735fa1
# ╠═c126f2d9-0ccb-4b5c-9600-e14ad96e9591
# ╠═7ccc9401-2fee-4d5e-a6ca-9a863c60d5a6
# ╠═9bd87f70-24e6-4476-acf7-7647face1b3e
# ╠═a40e5ab4-ff6c-4cbc-a377-db60b76d1963
# ╠═e4abed36-f47b-4d77-9d10-e7455edd9325
# ╠═202f831f-8e54-4d17-9f27-aea1a07ab919
# ╠═241a2416-6bc7-4804-b6ae-874d9336a02a
# ╠═23fabeb9-a947-4d44-9da2-a929cf9db6a4
# ╠═263448a0-a803-46b9-b7b5-d162e578f91f
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

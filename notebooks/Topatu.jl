### A Pluto.jl notebook ###
# v0.14.4

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

- Description of the TOPATU setup
- Calculations of the expected response in solution & silica
"

# ╔═╡ 4f35dcfc-8644-11eb-1c86-75608b38bf5f
#begin
#	using Pkg
#	Pkg.add.(["Unitful", "UnitfulEquivalences", "Plots","Interpolations",
#			"PhysicalConstants","Images", "ImageIO", "ImageView","Peaks","Formatting", #            "StrLiterals", "StrFormat", "Format","LsqFit"])
#end

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
projectdir()

# ╔═╡ 621ec96c-86fd-11eb-1c41-379cc17180dc
datadir()

# ╔═╡ 9b853f27-4288-42a6-8f12-ca004e1773b7
srcdir()

# ╔═╡ 04c59c63-ec43-4472-8947-ab277435932e
begin
	lfi = ingredients(srcdir("LabFbi.jl"))
	lti = ingredients(srcdir("LabTools.jl"))
	lpi = ingredients(srcdir("LabPlots.jl"))
end

# ╔═╡ 7c1ea236-c262-4dcc-9e2f-96c1ddf4f0b9
#begin
#	include(srcdir("LabFbi.jl"))
#	include(srcdir("LabTools.jl"))
#	include(srcdir("LabPlots.jl"))
#end

# ╔═╡ 0b1cedc7-8ada-45a5-adef-fbae794dee3e
markercolors = [:green :orange :black :purple :red  :yellow :brown :white]

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
	mW, μW, W,
    A, N, mol, mmol, V, L, mL, μL, M

# ╔═╡ 631833b8-0649-4f45-b1dd-e4a53091d321
md"#### Avogadro Number"

# ╔═╡ f8d3f48a-904f-11eb-2659-51f7508b434d
import PhysicalConstants.CODATA2018: N_A

# ╔═╡ 9fdae9e2-6d2e-4901-8f48-14764b6800c2
N_A

# ╔═╡ 70525ca4-8643-11eb-121e-196d80a60fcf
md"## TOPATU Setup"

# ╔═╡ e5e02216-86fb-11eb-1cf6-6f6962313a09
load("laserFBI.png")

# ╔═╡ 9bf0b7e6-8643-11eb-142b-7fc904562cc8
md"### Laser

- Topatu has a Blue laser, with wavelength 405 nm and power up to 100 mW
"

# ╔═╡ e9076d72-8648-11eb-0c02-9b00b0605e72
l400P1mw = lfi.LabFbi.Laser(405.0nm, 1.0mW);

# ╔═╡ 052941b7-6a8e-4930-bce7-7deec4c13ff0
begin
	Lx = collect(250:800) * nm
	Ex = lfi.LabFbi.photon_energy.(Lx)
end

# ╔═╡ 8ec6de87-0f1e-49a4-be7d-c3d9b6233718
lpi.LabPlots.plot_xy(Lx/nm, Ex/eV, 
	        "λ (nm)", "E (eV)", "Photon energy as a function of wavelength")


# ╔═╡ bbdce5e5-7fa5-4b06-a047-9346a7c23bb3
begin
	Px = collect(lti.LabTools.logrange(10^-2, 10^3, 100)) *mW
	Np = lfi.LabFbi.n_photons.((405nm,), Px)
end

# ╔═╡ 6b4e68f6-1e92-42bc-8113-a5797cf294bd
lpi.LabPlots.loglog_xy(Px/mW, Np/Hz, 
	        "P (mW)", "N", "Number of photons as a function of Power")


# ╔═╡ bb4fa686-908b-11eb-0aa0-27ab50bb0a4e
md"#### delivered energy: Set the time below"

# ╔═╡ 1c31637a-8666-11eb-2ea6-17b0d4507e10
@bind tt NumberField(1:10; default=1)

# ╔═╡ e0bbf162-fea4-427b-9133-569b51632746
md"- time = $tt second"

# ╔═╡ 28c0a218-8666-11eb-291a-71fd93c60da5
de = uconvert(mJ, lfi.LabFbi.delivered_energy(l400P1mw, tt * 1.0s))

# ╔═╡ 943550b1-7a39-4e99-9054-4b8e1c584b6f
md"### Objective"

# ╔═╡ 0bddff38-99ad-475e-9515-579de07e8ab3
obj = lfi.LabFbi.Objective("Topatu", 0.6, 75.0)

# ╔═╡ 77a68d3b-8fc5-4ef6-bb20-3140779e6535
gl400P1mw = lfi.LabFbi.GaussianLaser(l400P1mw, obj)

# ╔═╡ bb2e22d4-a54e-4519-9a65-89dd1f82a77b
begin
	w0s = lti.LabTools.to_fstr(gl400P1mw.w0/nm,"%5.2g")
	zrs = lti.LabTools.to_fstr(gl400P1mw.zr/nm,"%5.2g")
end

# ╔═╡ 2cd0e291-da66-489a-90af-38ce273c950d
md"### Diffractive limit"

# ╔═╡ 7f081083-3b93-4100-ac1f-4866708fad77
dl = lfi.LabFbi.diffraction_limit(gl400P1mw)

# ╔═╡ 90fea6e7-9803-4833-8f66-f339ed9dd9b2
md"### Gaussian Laser 
- The TOPATU laser is approximately gaussian.
- A gaussian laser is defined with a laser and its objective.
- The waist of the laser is $w0s nm
- The waist in z (zr) is $zrs nm
- The diffraction limit is $(lti.LabTools.to_fstr(dl/nm,\"%5.2g\")) nm
"

# ╔═╡ 586958c7-c69a-483a-b7c8-eb065db96f3b
md"### Field of View (FOV)
- Corresponds to the region (area, volume) iluminated by the laser
- It is defined by a diameter and a length (or thickness). 
- For example, for a gaussian laser beam, the FOV can be defined in terms of w0 and zr, or in terms of the diffraction limit
"

# ╔═╡ 4589e4cc-dd61-4239-b027-b135c3525f42
fovdl = lfi.LabFbi.Fov(dl, 2*gl400P1mw.zr)

# ╔═╡ 30371eba-26ac-429a-bb2e-30de65091288
md"### Spatial distribution of a guassian beam"

# ╔═╡ b78e625e-4c71-4f1c-8711-cfe90c27519b
begin
	R = collect(lti.LabTools.logrange(1, 500, 100))
	gir0(z) = lfi.LabFbi.gI(gl400P1mw, z*nm, 0*nm)
	giz0(r) = lfi.LabFbi.gI(gl400P1mw, 0*nm, r*nm)
	piz0 = lpi.LabPlots.plot_xy(R, giz0.(R), 
		   "R(nm)", "IG(z0)", "Beam profile in R at Z= 0");
	pir0 = lpi.LabPlots.plot_xy(R, gir0.(R), 
		   "Z(nm)", "IG(r0)", "Beam profile in Z at R= 0");
	plot(piz0,pir0, layout = (1, 2), legend=false, fmt = :png)
end

# ╔═╡ f0b9c27a-9623-48ee-b096-6a3fb56bc35c
md" - Notice that essentially 100 % of the beam is contained in the FoV"

# ╔═╡ 013420a8-9129-11eb-0d7f-ab8fb7c70786
md"#### Photon density"

# ╔═╡ fe18d00d-3fe9-4bf4-b926-2dcc12d513ab
Ip = lfi.LabFbi.photon_density.((405nm,), Px, (fovdl.a,))

# ╔═╡ 55046ba7-f5a6-4d9a-a767-3dcdb85bb1dc
lpi.LabPlots.loglog_xy(Px/mW, Ip * μm^2/Hz, 
	        "P (mW)", "ρ", "Photon density in FoV as a function of Power")


# ╔═╡ 53f112f6-ba9f-4089-ae96-83baad162fcc
md"- Define the power of the laser (in mW)"

# ╔═╡ 990240b5-b30d-4c50-911c-42f01d7c98c3
@bind lp NumberField(10^-1:10^3; default=0.1)

# ╔═╡ d0c7c8d5-e94b-434b-8920-c81688dfa5de
begin
	Ixp = lfi.LabFbi.photon_density.(405nm, lp*mW, fovdl.a)
	fovas = lti.LabTools.to_fstr(fovdl.a/nm^2,"%5.2g")
	Is    = lti.LabTools.to_fstr(Ixp/(Hz*μm^-2), "%5.2g")
	lps   = lti.LabTools.to_fstr(lp, "%5.2g")
	end

# ╔═╡ 877170f1-4f5c-438d-8949-854d3e892a2e
md" 
- FoV diameter        = $(fovdl.d)  
- Fov area            = $fovas nm2
- Laser power         = $lps mW
- Photon density       = $Is (Hz/``\mu m^{2}``)
"

# ╔═╡ a80e6642-908c-11eb-3ac3-6737f11f4eff
md"## Fluorescence in solution"

# ╔═╡ 9cb2c012-608a-4cff-b312-88bd2071da1f
md"#### Molar extinction coefficient
- λ in nm
- ϵ in M``^{-1}``cm``^{-1}``
"

# ╔═╡ cd707e24-df49-4f74-9890-0cdfcc1b9296
adf = lfi.LabFbi.load_df_from_csv(datadir("fbi"),
	                          "molar_extinction_coefficient_G1G2.csv", 
	                          lfi.LabFbi.enG) 


# ╔═╡ 482b76cf-b3fb-4a30-8678-56ebe539fc5e
md"##### Define fluorophores"

# ╔═╡ a35dd220-5ab1-4f70-901c-6a55f8e50123
begin
	fbi325G1 = lfi.LabFbi.Fluorophore("FBI-solution", 325nm, 489nm, 
		lfi.LabFbi.select_element(adf, "λ", "325","ϵFbiG1")/(M*cm), 0.67)
	fbiba325G1 = lfi.LabFbi.Fluorophore("FBIBa-solution", 325nm, 489nm, 
		lfi.LabFbi.select_element(adf, "λ", "325","ϵFbiBaG1")/(M*cm), 0.67)
	fbi405G1 = lfi.LabFbi.Fluorophore("FBI-solution", 405nm, 489nm, 
		lfi.LabFbi.select_element(adf, "λ", "405(500)","ϵFbiG1")/(M*cm), 0.67)
	fbiba405G1 = lfi.LabFbi.Fluorophore("FBIBa-solution", 325nm, 489nm, 
		lfi.LabFbi.select_element(adf, "λ", "405(500)","ϵFbiBaG1")/(M*cm), 0.67)
end

# ╔═╡ aaf97203-1412-493f-94cb-75591ff981f7
md"#### Emitted fluorescence

- Laser power    = $lp mW
- Photon density = $Is (Hz/``\mu m^{2}``)
"

# ╔═╡ 4ed47097-0d51-4639-a19f-38bbaa6c508d
begin
	γfbi325G1 = lfi.LabFbi.fluorescence(fbi325G1, Ixp)
	γfbiba325G1 = lfi.LabFbi.fluorescence(fbiba325G1, Ixp)
	γfbi405G1 = lfi.LabFbi.fluorescence(fbi405G1, Ixp)
	γfbiba405G1 = lfi.LabFbi.fluorescence(fbiba405G1, Ixp)
end

# ╔═╡ 45fd4819-12bd-449e-9f37-5885cf3604ca
begin
	γfbi325G1s = lti.LabTools.to_fstr(γfbi325G1/Hz,"%7.3g")
	γfbiba325G1s = lti.LabTools.to_fstr(γfbiba325G1/Hz,"%7.3g")
	γfbi405G1s = lti.LabTools.to_fstr(γfbi405G1/Hz,"%7.3g")
	γfbiba405G1s = lti.LabTools.to_fstr(γfbiba405G1/Hz,"%7.3g")
end

# ╔═╡ 7f595bfa-71d6-42f3-913c-b7574c0d4fc2
md""" #### Fluorescence per molecule in FoV

| λ (nm) | Fbi G1 (Hz) | FbiBa G1 (Hz)
|:---------- | ---------- |:------------:|
| 325    | $γfbi325G1s| $γfbiba325G1s | 
| 405    | $γfbi405G1s| $γfbiba405G1s   | 
"""

# ╔═╡ ccb649f8-0516-4a64-8ccc-5af13a7fe6e1
md""" #### Concentration in solution"""

# ╔═╡ a6db60d8-371d-4b93-b164-38d44f00d17a
md"- Set units (in M)"

# ╔═╡ 96639a5e-81c2-4eac-837c-fff5bc7ea5aa
@bind csu NumberField(1.0:10.0; default=5.0)

# ╔═╡ 639922e2-c0e5-48ee-a2b2-9e3afd543e1e
md"- Set power (in M)"

# ╔═╡ 205870db-bf67-4ece-8de4-239aaa5edc47
@bind csp NumberField(1e-5:1e-9; default=1e-5)

# ╔═╡ 21c6d2f8-519e-41ff-8895-ae3c77eb7a02
cs = csu * csp ;

# ╔═╡ 5347f077-66b4-45e6-a2f1-b0e6f52b0768
solfbi = lfi.LabFbi.Solution("FBI solution", cs*M)

# ╔═╡ 087da75f-f6a3-42f7-a5ee-c36d25326dcd
begin
	ρfbi = lfi.LabFbi.nofv(solfbi, fovdl.v)
	ρfbis = lti.LabTools.to_fstr(ρfbi, "%5.2g") 
end

# ╔═╡ 8b918c37-dcc2-4f23-b0e5-bced29d3c654
md"
- The concentration of FBI in solution is $(lti.LabTools.to_fstr(cs, \"%7.1g\")) M
- This corresponds to $ρfbis in the diffractive spot volume
"

# ╔═╡ ea543fe2-fc1d-42fe-b563-8fca369b1974
begin
	γfbi325G1Vol   = γfbi325G1 * ρfbi
	γfbiba325G1Vol = γfbiba325G1 * ρfbi
	γfbi405G1Vol   = γfbi405G1 * ρfbi
	γfbiba405G1Vol = γfbiba405G1 * ρfbi
end

# ╔═╡ eda71706-2771-47a6-8435-0274ef4b77cd
begin
	γfbi325G1Vs = lti.LabTools.to_fstr(γfbi325G1Vol/Hz,"%7.3g")
	γfbiba325G1Vs = lti.LabTools.to_fstr(γfbiba325G1Vol/Hz,"%7.3g")
	γfbi405G1Vs = lti.LabTools.to_fstr(γfbi405G1Vol/Hz,"%7.3g")
	γfbiba405G1Vs = lti.LabTools.to_fstr(γfbiba405G1Vol/Hz,"%7.3g")
end

# ╔═╡ a7f980e2-4910-464f-8b83-155fcfba5501
md""" #### Fluorescence emitted by the solution

| λ (nm) | Fbi G1 (Hz) | FbiBa G1 (Hz)
|:---------- | ---------- |:------------:|
| 325    | $γfbi325G1Vs| $γfbiba325G1Vs | 
| 405    | $γfbi405G1Vs| $γfbiba405G1Vs   | 
"""

# ╔═╡ afe6ab7a-9133-11eb-1f9f-dd7c0b27e4d0
md"### Number of photons reaching CCDs"

# ╔═╡ e2525056-d023-49c8-b086-f9a6debecb9e
begin
	naL = collect(0.1:0.05:1)
	obL = [lfi.LabFbi.Objective("Nikon-SLVD", na, 100) for na in naL]
	trL = [lfi.LabFbi.transmission(ob) for ob in obL]

	plot(naL, trL, leg=false, lw=2)

	xlabel!("numerical aperture")
	ylabel!("Transmission")
	title!("Tranmission of the objective as a function of NA")
end

# ╔═╡ 67b2b186-9139-11eb-2f29-f9f36be38b6f
begin
	to = lfi.LabFbi.transmission(obj)
	tos = lti.LabTools.to_fstr(to, "%5.2g") 
end

# ╔═╡ 38979d58-a88f-4e3c-ba0e-e0ed7675054d
md" 1. TOPATU objective:
- NA           = $(obj.NA) 
- Transmision  = $tos"

# ╔═╡ a7bd7368-6bf3-4bac-9483-3b61528f20c7
begin
	nfbi325G1Vol   = γfbi325G1Vol * to
	nfbiba325G1Vol = γfbiba325G1Vol * to
	nfbi405G1Vol   = γfbi405G1Vol * to
	nfbiba405G1Vol = γfbiba405G1Vol * to
	
	nVs= lti.LabTools.to_fstr.([nfbi325G1Vol,
							nfbiba325G1Vol,
							nfbi405G1Vol,
							nfbiba405G1Vol]./Hz,"%7.3g")
end

# ╔═╡ 47e4f919-3b72-47b9-8389-4a62587bb76b
md""" #### Fluorescence transmitted by the objective

| λ (nm) | Fbi G1 (Hz) | FbiBa G1 (Hz)
|:---------- | ---------- |:------------:|
| 325    | $(nVs[1])| $(nVs[2]) | 
| 405    | $(nVs[3])| $(nVs[4])   | 
"""

# ╔═╡ f96e9b12-9139-11eb-08cb-5319bca7a01d
md"## Topatu setup
- Laser
- Objective
- Filters
- CCD
"

# ╔═╡ 4624886f-f61b-4eaa-8a08-f275d74e9ac5
md"### CCD"

# ╔═╡ 79af9066-8680-11eb-0d91-15759e995f6c
orca = lfi.LabFbi.ccd()

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

# ╔═╡ 920bf9b9-7422-46cf-a772-1ee05e814474
begin
	λmin = 350.0nm
	λmax = 800.0nm
end

# ╔═╡ cd53b04b-749a-47b3-93fe-0ad194c4f053
Filters = lfi.LabFbi.load_filter.(lfi.LabFbi.fnames, (datadir("filters"),));

# ╔═╡ e2d92deb-8ac2-4991-a842-d9906316d6a8
fdfs, fints = collect(zip(Filters...));

# ╔═╡ ae0c0e93-94eb-4249-a1a5-0e06071d2620
lfi.LabFbi.fnames

# ╔═╡ fc629ff0-471e-46eb-8c5e-888f02801d8e
fbp405, flp425, fbp430, flp450, fdn405_522, fn405 = fints

# ╔═╡ 48753733-a844-4bbf-9ad2-8ffe7a2086f9
begin
	bp405p = lpi.LabPlots.plot_filterset(fbp405, "BandPass405", 390:420)
	bp430p = lpi.LabPlots.plot_filterset(fbp430, "BandPass430", 420:440);
	lp425p = lpi.LabPlots.plot_filterset(flp425, "LongPass425", 420:520);
	lp450p = lpi.LabPlots.plot_filterset(flp450, "LongPass450", 440:540);
	n405p = lpi.LabPlots.plot_filterset(fn405, "NF405", 390:420);
	dn405_522p = lpi.LabPlots.plot_filterset(fdn405_522, "DoubleNotch405_522", 390:450);
	true 
end

# ╔═╡ ce01c6b0-ad7f-49da-86a5-92ce67a7ee9a
plot(bp405p, bp430p, lp425p, lp450p, n405p, dn405_522p, layout = (3, 2), legend=false, fmt = :png)

# ╔═╡ b0b10062-93c9-11eb-1b87-0f76a1735fa1
md"## Laser transport

The TOPATU setup has two configurations: 
- Band-430: Select light in the band 430 +- 10. This is done combining:

  - The 405 nm BP filter (not needed for this calculation, since the beam is filtered before reaching the optical path, the filter contributes, together with other factors to a normalisation factor regarding power that needs to be computed)
  - The 425 nm LP filter, which enters twice in the optical path and eliminates high frequencies (below 425 nm).
  - The double noth which selects a region between 405 and 522 nm 
  - The 430 nm BP filter, which selects the signal (FbIBa2+ peak) region

- LP-450: Select light in the long pass Band 450 nm. This is done using the same combination than before except the last filter which is a 450 nm LP filter. This selects primarily the region of background (FBI).

- In addition one needs to add the efficiency of the CCD to have the complete transport function for both cases
"

# ╔═╡ 5eaa4bc1-25d4-4a0e-837d-4c1f24f868c1
fdn405 = fdn405_522

# ╔═╡ c126f2d9-0ccb-4b5c-9600-e14ad96e9591
filterSetBand430(λ) = map(x->flp425(x)^2 * fdn405(x) * fbp430(x) * orca(x), λ)

# ╔═╡ 7ccc9401-2fee-4d5e-a6ca-9a863c60d5a6
filterSetLP450(λ) = map(x->flp425(x)^2 * fdn405(x) * flp450(x) * orca(x), λ)

# ╔═╡ a40e5ab4-ff6c-4cbc-a377-db60b76d1963
@test flp425(460.0)^2 * fdn405(460.0) * flp450(460.0) * orca(460.0) ≈ filterSetLP450(460.0)

# ╔═╡ e4abed36-f47b-4d77-9d10-e7455edd9325
pfilterSetBand430 =lpi.LabPlots.plot_filterset(filterSetBand430, "filterSetBand430", 400.0:500.0);

# ╔═╡ 202f831f-8e54-4d17-9f27-aea1a07ab919
pfilterSetLP450 =lpi.LabPlots.plot_filterset(filterSetLP450, "filterSetLP450", 400.0:700.0);

# ╔═╡ 33b05bd8-12a6-47bf-b229-f4b97ea279d0
plot(pfilterSetBand430, pfilterSetLP450, layout = (1, 2), legend=false, fmt = :png)

# ╔═╡ 92196df1-1448-41aa-afd1-0c4de5dc0dc8
md"### Beam reflection

- Set below the reflectivity of the beam in the substrate at 405 nm. 
"

# ╔═╡ 78f98759-f580-4d42-85d1-fde4662b666c
@bind Rs NumberField(0.01:0.1:1; default=0.07)

# ╔═╡ b8fe9096-adfd-46f9-84ce-9be60ec49eff
md"- The reflectivity of the beam is set to Rs =$Rs"

# ╔═╡ f273fd6a-e1eb-4d80-86ab-ecdb43be3e53
t405bp430 = Rs * filterSetBand430(405.0)

# ╔═╡ df379f4e-c27b-422d-a1f3-1373cd1f7e76
t405lp430 = Rs * filterSetLP450(405.0)

# ╔═╡ 3d982cfb-7c5e-422a-8120-b6048b3cd658
mdt405lp430 = lti.LabTools.to_fstr(t405bp430, "%7.1g")

# ╔═╡ da582c1b-0f85-48cf-afd1-61e8a83e6c48
md" - The tranmission in the 430 band for 405 nm is $mdt405lp430. It follows that the beam reflection must be negligible."

# ╔═╡ 241a2416-6bc7-4804-b6ae-874d9336a02a
md"### Filtered pdfs for FBI and FBIBa2+"

# ╔═╡ f45d93f2-cc4d-4395-b9f3-f9a49983b8e2
begin
    wf = 350.
    we = 800.
    ws = 2.
    wl = wf:ws:we
end

# ╔═╡ 12d79258-842a-45dd-bc09-fcad9f9ad7be
dffbi = lfi.LabFbi.load_df_from_csv(datadir("fluorimeter/325nm"), "FBI_G1_vs_FBI_G2_em325nm.csv", lfi.LabFbi.spG)

# ╔═╡ 8f52fa68-7f07-4d06-8a5e-9c89708f29cd
pfbi = lpi.LabPlots.plotdf_xy(dffbi,"λ", "FBIG1","λ (nm)", "I (au)", color=:green);

# ╔═╡ cbbc1b49-967f-4ea1-8a0f-7afed8e8858c
pfbiba = lpi.LabPlots.plotdf_xy(dffbi,"λ", "FBIBaG1","λ (nm)", "I (au)",  
	color=:blue);

# ╔═╡ e66c1dab-c4eb-4d49-9785-5f9a72a3b16f
lpi.LabPlots.merge_plots!(pfbi, pfbiba)

# ╔═╡ 24800b8e-34f4-4169-848d-1e2c52e0600b
gfbi = lfi.LabFbi.dftogf(wl, dffbi, "FBIG1")

# ╔═╡ 23a4ab0e-36ce-4c89-81b7-81aa7da74bff
gfbiba = lfi.LabFbi.dftogf(wl, dffbi, "FBIBaG1")

# ╔═╡ d7c27f64-fccb-4ca1-8e25-a8d4a91ff21c
pfbif =lpi.LabPlots.plotdf_gfs([gfbi,gfbiba], wl, ["FBI", "FBIBa"], [:green, :blue],
	                    "λ (nm)", "I (au)", "FBIG1");

# ╔═╡ dd73c8a9-a9fd-491d-9b6d-2ee926588d7e
pfbipdf = lpi.LabPlots.plotdf_gfs([gfbi,gfbiba], wl, ["FBI", "FBIBa"], [:green, :blue],
	                    "λ (nm)", "I (au)", "FBIG1", pdf=false);

# ╔═╡ 28a7d41b-159e-4753-ae99-f4d1a72a7e54
plot(pfbif, pfbipdf, layout = (1, 2), legend=false, fmt = :png)

# ╔═╡ 85dddc12-0a69-4745-a895-2508ea3ffc5c
filteredFbi430BP(λ) = map(x->filterSetBand430(x) * gfbi.pdf(x), λ)

# ╔═╡ e72d34ff-1a0a-4dd4-a088-49e671d01408
filteredFbiBa430BP(λ) = map(x->filterSetBand430(x) * gfbiba.pdf(x), λ)

# ╔═╡ 60badadd-4df8-426a-b367-86e7bdc0c25f
filteredFbi450LP(λ) = map(x->filterSetLP450(x) * gfbi.pdf(x), λ)

# ╔═╡ 6a5f6e6f-06cb-46cd-bab4-f3a0778af2b6
filteredFbiBa450LP(λ) = map(x->filterSetLP450(x) * gfbiba.pdf(x), λ)

# ╔═╡ 29ad555a-932b-4ee2-96c0-1e5550f7a1fd
begin
	n325FbiG1F430BP(λ) = filteredFbi430BP(λ) * nfbi325G1Vol/Hz
	n325FbiBaG1F430BP(λ) = filteredFbiBa430BP(λ) * nfbiba325G1Vol/Hz
	n325FbiG1F450LP(λ) = filteredFbi450LP(λ) * nfbi325G1Vol/Hz
	n325FbiBaG1F450LP(λ) = filteredFbiBa450LP(λ) * nfbiba325G1Vol/Hz

	n405FbiG1F430BP(λ) = filteredFbi430BP(λ) * nfbi405G1Vol/Hz
	n405FbiBaG1F430BP(λ) = filteredFbiBa430BP(λ) * nfbiba405G1Vol/Hz
	n405FbiG1F450LP(λ) = filteredFbi450LP(λ) * nfbi405G1Vol/Hz
	n405FbiBaG1F450LP(λ) = filteredFbiBa450LP(λ) * nfbiba405G1Vol/Hz
end

# ╔═╡ 67709bf1-298a-4371-8f91-1b6fab98ef6b
begin
	pfilteredFbi430BP =lpi.LabPlots.plot_filterset(filteredFbi430BP, "filteredFbi430BP", 400.0:700.0);
	pfilteredFbiBa430BP =lpi.LabPlots.plot_filterset(filteredFbiBa430BP, "filteredFbiBa430BP", 400.0:700.0);
	pfilteredFbi450LP = lpi.LabPlots.plot_filterset(filteredFbi450LP, "filteredFbi450LP", 400.0:700.0);
	pfilteredFbiBa450LP =lpi.LabPlots.plot_filterset(filteredFbiBa450LP, "filteredFbiBa450LP", 400.0:700.0);
	true
end

# ╔═╡ f613b185-0f04-4926-8e7d-f1d9d19046b5
plot(pfilteredFbi430BP, pfilteredFbiBa430BP, pfilteredFbi450LP, pfilteredFbiBa450LP, layout = (2, 2), legend=false, fmt = :png)

# ╔═╡ 3ab61c40-38a1-4654-99c5-6d5f6f36ea85
begin
	pn325FbiG1F430BP =lpi.LabPlots.plot_filterset(n325FbiG1F430BP, 
		"n325FbiG1F430BP", 400.0:700.0)
	pn325FbiBaG1F430BP =lpi.LabPlots.plot_filterset(n325FbiBaG1F430BP, "n325FbiBaG1F430BP", 400.0:700.0)
	pn325FbiG1F450LP = lpi.LabPlots.plot_filterset(n325FbiG1F450LP, "n325FbiG1F450LP", 400.0:700.0)
	pn325FbiBaG1F450LP =lpi.LabPlots.plot_filterset(n325FbiBaG1F450LP, "n325FbiBaG1F450LP", 400.0:700.0)
	true
end

# ╔═╡ 74ad1f8d-1a38-48de-a158-7f4c09bc7031
begin
	pn405FbiG1F430BP =lpi.LabPlots.plot_filterset(n405FbiG1F430BP, 
		"n405FbiG1F430BP", 400.0:700.0)
	pn405FbiBaG1F430BP =lpi.LabPlots.plot_filterset(n405FbiBaG1F430BP, "n405FbiBaG1F430BP", 400.0:700.0)
	pn405FbiG1F450LP = lpi.LabPlots.plot_filterset(n405FbiG1F450LP, "n405FbiG1F450LP", 400.0:700.0)
	pn405FbiBaG1F450LP =lpi.LabPlots.plot_filterset(n405FbiBaG1F450LP, "n405FbiBaG1F450LP", 400.0:700.0)
	true
end

# ╔═╡ 6b95832a-7c78-4ec9-891b-aaae5e382411
plot(pn325FbiG1F430BP, pn325FbiBaG1F430BP, pn325FbiG1F450LP, pn325FbiBaG1F450LP, layout = (2, 2), legend=false, fmt = :png)

# ╔═╡ 53dd2eaf-a0ac-41aa-a1c8-eb8c6daa9b17
plot(pn405FbiG1F430BP, pn405FbiBaG1F430BP, pn405FbiG1F450LP, pn405FbiBaG1F450LP, layout = (2, 2), legend=false, fmt = :png)

# ╔═╡ 6f9c9bfd-b533-4e0a-9a55-8cebaf773660
begin
	ccd325FbiG1F430BP = lfi.LabFbi.qpdf(n325FbiG1F430BP, 420.0, 440.0)
	ccdn325FbiBaG1F430BP = lfi.LabFbi.qpdf(n325FbiBaG1F430BP, 420.0, 440.0)
	ccdn325FbiG1F450LP = lfi.LabFbi.qpdf(n325FbiG1F450LP, 450.0, 700.0)
	ccdn325FbiBaG1F450LP = lfi.LabFbi.qpdf(n325FbiBaG1F450LP, 450.0, 700.0)

	ccd405FbiG1F430BP = lfi.LabFbi.qpdf(n405FbiG1F430BP, 420.0, 440.0)
	ccdn405FbiBaG1F430BP = lfi.LabFbi.qpdf(n405FbiBaG1F430BP, 420.0, 440.0)
	ccdn405FbiG1F450LP = lfi.LabFbi.qpdf(n405FbiG1F450LP, 450.0, 700.0)
	ccdn405FbiBaG1F450LP = lfi.LabFbi.qpdf(n405FbiBaG1F450LP, 450.0, 700.0)

end

# ╔═╡ d5b81d54-ba97-43d8-852b-e511f65949ac
nXs325= lti.LabTools.to_fstr.([ccd325FbiG1F430BP,
							ccdn325FbiBaG1F430BP,
							ccdn325FbiG1F450LP,
							ccdn325FbiBaG1F450LP],"%7.3g")



# ╔═╡ 8925c81a-6c09-4ff4-bd05-ea4137a5752a
nXs405= lti.LabTools.to_fstr.([ccd405FbiG1F430BP,
							ccdn405FbiBaG1F430BP,
							ccdn405FbiG1F450LP,
							ccdn405FbiBaG1F450LP],"%7.3g")

# ╔═╡ 31c984a7-2a4f-4941-8824-3042d8bdfcb0
md""" #### Fluorescence recorded by the CCD (in pes) λ = 325 nm 

| Filter | Fbi G1 (Hz) | FbiBa G1 (Hz)
|:---------- | ---------- |:------------:|
| BP430    | $(nXs325[1])| $(nXs325[2]) | 
| LP450    | $(nXs325[3])| $(nXs325[4])   | 
"""

# ╔═╡ a5e0febb-1bb2-4f01-acfa-927c2dcb433b
md""" #### Fluorescence recorded by the CCD (in pes) λ = 405 nm 

| Filter | Fbi G1 (Hz) | FbiBa G1 (Hz)
|:---------- | ---------- |:------------:|
| BP430    | $(nXs405[1])| $(nXs405[2]) | 
| LP450    | $(nXs405[3])| $(nXs405[4])   | 
"""

# ╔═╡ 30156129-845e-4261-b549-e6c576c749ef
begin
	fi430bp = quadgk(filteredFbi430BP, 375.0,  700.0)[1]
	fi450lp = quadgk(filteredFbi450LP, 375.0,  700.0)[1]
	rfbibplp = fi430bp / fi450lp
end

# ╔═╡ bc05ba67-4fdc-48bb-8985-6a9f62718f2e
rfbibplpN = ccd325FbiG1F430BP / ccdn325FbiG1F450LP

# ╔═╡ dc9773bd-b1f0-4262-9a6a-d8636fb5fe2b
rfbibpFbiBabp = ccdn405FbiBaG1F430BP / ccd405FbiG1F430BP

# ╔═╡ ea5bd4b0-9a3a-461c-9389-463c375b7c98
rFbiBabplpN = ccdn405FbiBaG1F430BP / ccdn405FbiBaG1F450LP

# ╔═╡ 31b1020e-a8e2-44fd-a08b-e49dca478d22
rx405= lti.LabTools.to_fstr.([rfbibplpN,
							rfbibpFbiBabp,
							rFbiBabplpN],"%7.3g")


# ╔═╡ b71c2f66-52ae-466c-8241-36fdfdeacb88
md"
- Ratio BP/LP for Fbi       = $(rx405[1])
- Ratio BP/BP for FbiBa/Fbi = $(rx405[3])
- Ratio BP/LP for FbiBa     = $(rx405[2])
"

# ╔═╡ Cell order:
# ╠═79cfd2fc-9046-11eb-2b13-1b877d57645d
# ╠═4f35dcfc-8644-11eb-1c86-75608b38bf5f
# ╠═3207f446-8643-11eb-37ba-c9aec47fcb8f
# ╠═5115917a-8644-11eb-19fc-0528741ca75d
# ╠═fc8b79a2-8728-11eb-2da7-e3ffa3ceef08
# ╠═68e738e7-88bd-41c2-89e4-594f07d64ddc
# ╠═0f2f4c78-8729-11eb-2bab-27812ce8c47e
# ╠═5ee27d52-86fd-11eb-365e-9f2b1a095575
# ╠═621ec96c-86fd-11eb-1c41-379cc17180dc
# ╠═9b853f27-4288-42a6-8f12-ca004e1773b7
# ╠═04c59c63-ec43-4472-8947-ab277435932e
# ╠═7c1ea236-c262-4dcc-9e2f-96c1ddf4f0b9
# ╠═0b1cedc7-8ada-45a5-adef-fbae794dee3e
# ╟─3c01d7c8-8645-11eb-372d-4d66ae46ae72
# ╠═46b5465a-8645-11eb-0291-612455795518
# ╟─631833b8-0649-4f45-b1dd-e4a53091d321
# ╠═f8d3f48a-904f-11eb-2659-51f7508b434d
# ╠═9fdae9e2-6d2e-4901-8f48-14764b6800c2
# ╟─70525ca4-8643-11eb-121e-196d80a60fcf
# ╠═e5e02216-86fb-11eb-1cf6-6f6962313a09
# ╟─9bf0b7e6-8643-11eb-142b-7fc904562cc8
# ╠═e9076d72-8648-11eb-0c02-9b00b0605e72
# ╠═052941b7-6a8e-4930-bce7-7deec4c13ff0
# ╠═8ec6de87-0f1e-49a4-be7d-c3d9b6233718
# ╠═bbdce5e5-7fa5-4b06-a047-9346a7c23bb3
# ╠═6b4e68f6-1e92-42bc-8113-a5797cf294bd
# ╟─bb4fa686-908b-11eb-0aa0-27ab50bb0a4e
# ╠═1c31637a-8666-11eb-2ea6-17b0d4507e10
# ╟─e0bbf162-fea4-427b-9133-569b51632746
# ╠═28c0a218-8666-11eb-291a-71fd93c60da5
# ╠═943550b1-7a39-4e99-9054-4b8e1c584b6f
# ╠═0bddff38-99ad-475e-9515-579de07e8ab3
# ╠═90fea6e7-9803-4833-8f66-f339ed9dd9b2
# ╠═77a68d3b-8fc5-4ef6-bb20-3140779e6535
# ╠═bb2e22d4-a54e-4519-9a65-89dd1f82a77b
# ╠═2cd0e291-da66-489a-90af-38ce273c950d
# ╠═7f081083-3b93-4100-ac1f-4866708fad77
# ╠═586958c7-c69a-483a-b7c8-eb065db96f3b
# ╠═4589e4cc-dd61-4239-b027-b135c3525f42
# ╠═30371eba-26ac-429a-bb2e-30de65091288
# ╠═b78e625e-4c71-4f1c-8711-cfe90c27519b
# ╠═f0b9c27a-9623-48ee-b096-6a3fb56bc35c
# ╠═013420a8-9129-11eb-0d7f-ab8fb7c70786
# ╠═fe18d00d-3fe9-4bf4-b926-2dcc12d513ab
# ╠═55046ba7-f5a6-4d9a-a767-3dcdb85bb1dc
# ╠═53f112f6-ba9f-4089-ae96-83baad162fcc
# ╟─990240b5-b30d-4c50-911c-42f01d7c98c3
# ╠═d0c7c8d5-e94b-434b-8920-c81688dfa5de
# ╠═877170f1-4f5c-438d-8949-854d3e892a2e
# ╠═a80e6642-908c-11eb-3ac3-6737f11f4eff
# ╠═9cb2c012-608a-4cff-b312-88bd2071da1f
# ╠═cd707e24-df49-4f74-9890-0cdfcc1b9296
# ╠═482b76cf-b3fb-4a30-8678-56ebe539fc5e
# ╠═a35dd220-5ab1-4f70-901c-6a55f8e50123
# ╠═aaf97203-1412-493f-94cb-75591ff981f7
# ╠═4ed47097-0d51-4639-a19f-38bbaa6c508d
# ╠═45fd4819-12bd-449e-9f37-5885cf3604ca
# ╠═7f595bfa-71d6-42f3-913c-b7574c0d4fc2
# ╠═ccb649f8-0516-4a64-8ccc-5af13a7fe6e1
# ╟─a6db60d8-371d-4b93-b164-38d44f00d17a
# ╟─96639a5e-81c2-4eac-837c-fff5bc7ea5aa
# ╟─639922e2-c0e5-48ee-a2b2-9e3afd543e1e
# ╟─205870db-bf67-4ece-8de4-239aaa5edc47
# ╠═21c6d2f8-519e-41ff-8895-ae3c77eb7a02
# ╠═8b918c37-dcc2-4f23-b0e5-bced29d3c654
# ╠═5347f077-66b4-45e6-a2f1-b0e6f52b0768
# ╠═087da75f-f6a3-42f7-a5ee-c36d25326dcd
# ╠═ea543fe2-fc1d-42fe-b563-8fca369b1974
# ╠═eda71706-2771-47a6-8435-0274ef4b77cd
# ╠═a7f980e2-4910-464f-8b83-155fcfba5501
# ╠═afe6ab7a-9133-11eb-1f9f-dd7c0b27e4d0
# ╠═e2525056-d023-49c8-b086-f9a6debecb9e
# ╠═38979d58-a88f-4e3c-ba0e-e0ed7675054d
# ╠═67b2b186-9139-11eb-2f29-f9f36be38b6f
# ╠═a7bd7368-6bf3-4bac-9483-3b61528f20c7
# ╠═47e4f919-3b72-47b9-8389-4a62587bb76b
# ╠═f96e9b12-9139-11eb-08cb-5319bca7a01d
# ╠═4624886f-f61b-4eaa-8a08-f275d74e9ac5
# ╠═79af9066-8680-11eb-0d91-15759e995f6c
# ╟─2f259d3a-8688-11eb-1733-13978ae756ce
# ╟─56d72e3c-86fa-11eb-2928-798e913a5976
# ╠═920bf9b9-7422-46cf-a772-1ee05e814474
# ╠═cd53b04b-749a-47b3-93fe-0ad194c4f053
# ╠═e2d92deb-8ac2-4991-a842-d9906316d6a8
# ╠═ae0c0e93-94eb-4249-a1a5-0e06071d2620
# ╠═fc629ff0-471e-46eb-8c5e-888f02801d8e
# ╠═48753733-a844-4bbf-9ad2-8ffe7a2086f9
# ╠═ce01c6b0-ad7f-49da-86a5-92ce67a7ee9a
# ╠═b0b10062-93c9-11eb-1b87-0f76a1735fa1
# ╠═5eaa4bc1-25d4-4a0e-837d-4c1f24f868c1
# ╠═c126f2d9-0ccb-4b5c-9600-e14ad96e9591
# ╠═7ccc9401-2fee-4d5e-a6ca-9a863c60d5a6
# ╟─a40e5ab4-ff6c-4cbc-a377-db60b76d1963
# ╠═e4abed36-f47b-4d77-9d10-e7455edd9325
# ╠═202f831f-8e54-4d17-9f27-aea1a07ab919
# ╠═33b05bd8-12a6-47bf-b229-f4b97ea279d0
# ╟─92196df1-1448-41aa-afd1-0c4de5dc0dc8
# ╟─78f98759-f580-4d42-85d1-fde4662b666c
# ╟─b8fe9096-adfd-46f9-84ce-9be60ec49eff
# ╠═f273fd6a-e1eb-4d80-86ab-ecdb43be3e53
# ╠═df379f4e-c27b-422d-a1f3-1373cd1f7e76
# ╠═3d982cfb-7c5e-422a-8120-b6048b3cd658
# ╟─da582c1b-0f85-48cf-afd1-61e8a83e6c48
# ╠═241a2416-6bc7-4804-b6ae-874d9336a02a
# ╠═f45d93f2-cc4d-4395-b9f3-f9a49983b8e2
# ╠═12d79258-842a-45dd-bc09-fcad9f9ad7be
# ╠═8f52fa68-7f07-4d06-8a5e-9c89708f29cd
# ╠═cbbc1b49-967f-4ea1-8a0f-7afed8e8858c
# ╠═e66c1dab-c4eb-4d49-9785-5f9a72a3b16f
# ╠═24800b8e-34f4-4169-848d-1e2c52e0600b
# ╠═23a4ab0e-36ce-4c89-81b7-81aa7da74bff
# ╠═d7c27f64-fccb-4ca1-8e25-a8d4a91ff21c
# ╠═dd73c8a9-a9fd-491d-9b6d-2ee926588d7e
# ╠═28a7d41b-159e-4753-ae99-f4d1a72a7e54
# ╠═85dddc12-0a69-4745-a895-2508ea3ffc5c
# ╠═e72d34ff-1a0a-4dd4-a088-49e671d01408
# ╠═60badadd-4df8-426a-b367-86e7bdc0c25f
# ╠═6a5f6e6f-06cb-46cd-bab4-f3a0778af2b6
# ╠═29ad555a-932b-4ee2-96c0-1e5550f7a1fd
# ╠═67709bf1-298a-4371-8f91-1b6fab98ef6b
# ╠═f613b185-0f04-4926-8e7d-f1d9d19046b5
# ╠═3ab61c40-38a1-4654-99c5-6d5f6f36ea85
# ╠═74ad1f8d-1a38-48de-a158-7f4c09bc7031
# ╠═6b95832a-7c78-4ec9-891b-aaae5e382411
# ╠═53dd2eaf-a0ac-41aa-a1c8-eb8c6daa9b17
# ╠═6f9c9bfd-b533-4e0a-9a55-8cebaf773660
# ╠═d5b81d54-ba97-43d8-852b-e511f65949ac
# ╠═8925c81a-6c09-4ff4-bd05-ea4137a5752a
# ╠═31c984a7-2a4f-4941-8824-3042d8bdfcb0
# ╠═a5e0febb-1bb2-4f01-acfa-927c2dcb433b
# ╠═30156129-845e-4261-b549-e6c576c749ef
# ╠═bc05ba67-4fdc-48bb-8985-6a9f62718f2e
# ╠═dc9773bd-b1f0-4262-9a6a-d8636fb5fe2b
# ╠═ea5bd4b0-9a3a-461c-9389-463c375b7c98
# ╠═31b1020e-a8e2-44fd-a08b-e49dca478d22
# ╠═b71c2f66-52ae-466c-8241-36fdfdeacb88

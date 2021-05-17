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
md"# TOPATU FBI 

- Analysis of the fluorimeter response of Fbi G1 & G2 in solution at different concentrations. 
- Calculations of the expected response in the TOPATU setup at 405 and 325 nm
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
    ns, μs, ms, ps, s, minute, hr, d, yr, Hz,
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
md"### Lasers available in TOPATU

- Blue laser, with wavelength 405 nm and power up to 100 mW
- UV   laser, with wavelength 325 nm and power up to 10 mW
"

# ╔═╡ c0fd3834-893c-45d5-8530-8758c8661e1e
md" - Set laser power (in mW)"

# ╔═╡ 990240b5-b30d-4c50-911c-42f01d7c98c3
@bind lp NumberField(10^-1:10^3; default=0.1)

# ╔═╡ 0519d9f0-a35a-4f5c-a70b-692809bee03e
md"- power = $lp mW"

# ╔═╡ e9076d72-8648-11eb-0c02-9b00b0605e72
l400 = lfi.LabFbi.Laser(405.0nm, lp*mW)

# ╔═╡ dad410f8-cca9-4f2e-8735-b2e1fb135ca8
l325 = lfi.LabFbi.Laser(325.0nm, lp*mW)

# ╔═╡ 052941b7-6a8e-4930-bce7-7deec4c13ff0
function plot_photon_energy_vs_wavelength(wl=250:800)
	Lx = collect(wl) * nm
	Ex = lfi.LabFbi.photon_energy.(Lx)
	pewl = lpi.LabPlots.plot_xy(Lx/nm, Ex/eV, 
	                            "λ (nm)", "E (eV)", "Photon energy")
	return pewl
end


# ╔═╡ 9f9b7846-2120-4c8e-bb45-61e00115d71b
function plot_number_photons_vs_power(laser, title)
	Px = collect(lti.LabTools.logrange(10^-2, 10^3, 100)) *mW
	Np = lfi.LabFbi.n_photons.((laser.λ,), Px)
	pnp = lpi.LabPlots.loglog_xy(Px/mW, Np/Hz, 
	        "P (mW)", "N", title)
	return pnp
end

# ╔═╡ c5f5887d-157c-42cb-a0ed-e3f82b553eec
begin
	pewl = plot_photon_energy_vs_wavelength()
	pnp405 = plot_number_photons_vs_power(l400, "Number of photons, laser 405 nm")
	pnp325 = plot_number_photons_vs_power(l325, "Number of photons, laser 325 nm")
	plot(pewl,pnp405,pnp325, layout = (3, 1), legend=false, fmt = :png)
end

# ╔═╡ bb4fa686-908b-11eb-0aa0-27ab50bb0a4e
md"#### delivered energy: Set the time below"

# ╔═╡ 1c31637a-8666-11eb-2ea6-17b0d4507e10
@bind tt NumberField(1:10; default=1)

# ╔═╡ e0bbf162-fea4-427b-9133-569b51632746
md"- time = $tt second"

# ╔═╡ 28c0a218-8666-11eb-291a-71fd93c60da5
de = uconvert(mJ, lfi.LabFbi.delivered_energy(l400, tt * 1.0s))

# ╔═╡ 943550b1-7a39-4e99-9054-4b8e1c584b6f
md"### Objective"

# ╔═╡ 0bddff38-99ad-475e-9515-579de07e8ab3
obj = lfi.LabFbi.Objective("Topatu", 0.6, 75.0);

# ╔═╡ 375afbf8-47bd-470d-93f2-715f13760501
md"## Laser setup"

# ╔═╡ 5bd97b62-5d4d-4483-8c05-2eec1d1e6b57
lT = lfi.LabFbi.LaserSetup(Dict("l405" => l400, "l325" => l325), lp*mW, obj)

# ╔═╡ f4fd8031-1b85-4f98-898f-eb09835ddc2d
function laser_info_md(lT, obj)
	md = " ### TOPATU setup \n
- Objective NA = $(obj.NA)\n"
	for lname in keys(lT.lasers)
		md = string(md, "- Laser name = $lname \n")
		lp = lti.LabTools.to_fstr(lT.lasers[lname].P/mW,"%5.2g")
		wst= lti.LabTools.to_fstr(lT.glasers[lname].w0/nm,"%5.2g")
		zr = lti.LabTools.to_fstr(lT.glasers[lname].zr/nm,"%5.2g")
		dl = lti.LabTools.to_fstr(lT.dls[lname]/nm,"%5.2g")
		fa = lti.LabTools.to_fstr(lT.fovs[lname].a/μm^2,"%5.2g")
		I  = lti.LabTools.to_fstr(lT.Is[lname]/(Hz* cm^-2),"%5.2g")
		md = string(md, "  - laser power = $lp mW \n")
		md = string(md, "  - laser waist = $wst nm \n")
		md = string(md, "  - laser zr    = $zr   nm \n")
		md = string(md, "  - diff limit  = $dl   nm \n")
		md = string(md, "  - Fov  area   = $fa   μm2 \n")
		md = string(md, "  - photon ρ    = $I   Hz/cm2 \n")
	end
	return md
end

# ╔═╡ 62da60bf-df76-41f9-9027-78953fbbcdb8
begin
	linfo = laser_info_md(lT, obj)
	Markdown.parse("$linfo")
end

# ╔═╡ 30371eba-26ac-429a-bb2e-30de65091288
md"### Spatial distribution 425 nm beam"

# ╔═╡ b78e625e-4c71-4f1c-8711-cfe90c27519b
begin
	R = collect(lti.LabTools.logrange(1, 500, 100))
	gir0(z) = lfi.LabFbi.gI(lT.glasers["l405"], z*nm, 0*nm)
	giz0(r) = lfi.LabFbi.gI(lT.glasers["l405"], 0*nm, r*nm)
	piz0 = lpi.LabPlots.plot_xy(R, giz0.(R), 
		   "R(nm)", "IG(z0)", "Beam profile in R at Z= 0");
	pir0 = lpi.LabPlots.plot_xy(R, gir0.(R), 
		   "Z(nm)", "IG(r0)", "Beam profile in Z at R= 0");
	plot(piz0,pir0, layout = (1, 2), legend=false, fmt = :png)
end

# ╔═╡ dd35eedb-9ed0-4bc2-9636-3192c1b034ab
md"### Spatial distribution 325 nm beam"

# ╔═╡ 3d182094-a3bb-413a-a1fc-6026839557c0
begin
	gir0325(z) = lfi.LabFbi.gI(lT.glasers["l325"], z*nm, 0*nm)
	giz0325(r) = lfi.LabFbi.gI(lT.glasers["l325"], 0*nm, r*nm)
	piz0325 = lpi.LabPlots.plot_xy(R, giz0325.(R), 
		   "R(nm)", "IG(z0)", "Beam profile in R at Z= 0");
	pir0325 = lpi.LabPlots.plot_xy(R, gir0325.(R), 
		   "Z(nm)", "IG(r0)", "Beam profile in Z at R= 0");
	plot(piz0325,pir0325, layout = (1, 2), legend=false, fmt = :png)
end

# ╔═╡ f0b9c27a-9623-48ee-b096-6a3fb56bc35c
md" - Notice that essentially 100 % of the beam is contained in the FoV"

# ╔═╡ 013420a8-9129-11eb-0d7f-ab8fb7c70786
md"#### Photon density"

# ╔═╡ fe18d00d-3fe9-4bf4-b926-2dcc12d513ab
begin
	Px = collect(lti.LabTools.logrange(10^-2, 10^3, 100)) *mW
	Ip = lfi.LabFbi.photon_density.((405nm,), Px, (lT.fovs["l405"].a,))
	Ip325 = lfi.LabFbi.photon_density.((325nm,), Px, (lT.fovs["l325"].a,))
	println("hi")
end

# ╔═╡ 55046ba7-f5a6-4d9a-a767-3dcdb85bb1dc
begin
	pdfov =lpi.LabPlots.loglog_xy(Px/mW, Ip * cm^2/Hz, 
	        "P (mW)", "ρ (Hz/cm2)", "Photon ρ in FoV (400 nm)");
	pdfov325 =lpi.LabPlots.loglog_xy(Px/mW, Ip325 * cm^2/Hz, 
	        "P (mW)", "ρ (Hz/cm2)", "Photon ρ in FoV (325 nm) ");
println("hi")
end


# ╔═╡ 41e148a6-5dec-49a2-bffa-c4b9d5816d77
plot(pdfov,pdfov325, layout = (1, 2), legend=false, fmt = :png)

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

# ╔═╡ fe77ba58-c130-4870-a7c0-180bd9af22bb
fbiFluo = lfi.LabFbi.fbi_fluorophores(adf, ["g1", "g2"], [325,405], [0.67,0.67])

# ╔═╡ 2861c9e2-ae5a-4a9b-91ad-88a1ea51b450
efluo = lfi.LabFbi.emitted_fluorescence(fbiFluo, lT)

# ╔═╡ 1e9740d1-1fad-4a0c-bfd4-a189857bd9f0
function efluo_info_md(efluo, title)
	FBI=Dict()
	FBIBA=Dict()
	gs = ["g1", "g2"]
	ls = [325,405]
	
	for gn in gs
		for l in ls
			lfbi = string("l", l, gn)
			FBI[lfbi] =  lti.LabTools.to_fstr(efluo.fbi[lfbi]/Hz,"%5.2g")
			FBIBA[lfbi] = lti.LabTools.to_fstr(efluo.fbiba[lfbi]/Hz,"%5.2g")
	   	end
	end

	md = """ ### $title \n
	
| λ (nm) | Fbi G1 (Hz) | FbiBa G1 (Hz)|Fbi G2 (Hz) | FbiBa G2 (Hz)
|:------:|:----------: |:------------:|:------------:|:------------:|
| 325    | $(FBI["l325g1"])| $(FBIBA["l325g1"]) | $(FBI["l325g2"]) | $(FBIBA["l325g2"]) |
| 405    | $(FBI["l405g1"])| $(FBIBA["l405g1"]) | $(FBI["l405g2"]) |$(FBIBA["l405g2"]) |

"""

	return md
end

# ╔═╡ cfe3217e-57b3-4cdb-ba09-5df3f194ceac
begin
	mdef = efluo_info_md(efluo, "Emitted fluorescence per molecule in FoV")
	Markdown.parse("$mdef")
end

# ╔═╡ ccb649f8-0516-4a64-8ccc-5af13a7fe6e1
md""" ### Concentration in solution"""

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

# ╔═╡ 61de7676-2d4f-40da-9250-92d00f8d780b
function molfov_info_md(cs, mfov)
	
	mdcs = lti.LabTools.to_fstr(cs, "%5.2g")
	mdfov = Dict()
	for lname in keys(mfov)
		mdfov[lname] = lti.LabTools.to_fstr(mfov[lname], "%5.2g")
	end
	
	md = """ ### Number of molecules \n
- The concentration of FBI in solution is $mdcs M
- This corresponds to:
  - laser 325 nm: $(mdfov["l325"]) molecules in the FoV volume
  - laser 405 nm: $(mdfov["l405"]) molecules in the FoV volume
"""
	return md
end

# ╔═╡ 5347f077-66b4-45e6-a2f1-b0e6f52b0768
solfbi = lfi.LabFbi.Solution("FBI solution", cs*M)

# ╔═╡ 8add778e-bffd-4934-a867-6f9658b4be89
mfov = lfi.LabFbi.molecules_in_fov(lT, solfbi)

# ╔═╡ 2acdf979-7a0a-4538-92bd-f31df54c503e
begin
	mdmfov = molfov_info_md(cs, mfov)
	Markdown.parse("$mdmfov")
end

# ╔═╡ 7d2e1f10-6a8b-46af-b2f9-9438ff3f2f89
efluo.fbiba

# ╔═╡ 32503979-9ef1-4b7a-aeb9-1a63b33f60cd
ffov = lfi.LabFbi.emitted_fluorescence_fov(efluo, mfov)

# ╔═╡ 1aefe05a-e218-473f-a9df-56eb3fd844cb
efluo

# ╔═╡ 65833e10-898f-44a9-8731-3e3df800a2f6
begin
	mdffov = efluo_info_md(ffov, "Emitted fluorescence ALL molecules in FoV")
	Markdown.parse("$mdffov")
end

# ╔═╡ 8e16733d-9b6f-4170-a05f-e23cc63b4a05
ofilter = lfi.LabFbi.transmission(obj);

# ╔═╡ 38979d58-a88f-4e3c-ba0e-e0ed7675054d
md""" ### TOPATU objective:
- NA           = $(obj.NA) 
- Transmision  = $(lti.LabTools.to_fstr(ofilter, "%5.2g"))
"""

# ╔═╡ b938a16d-3a15-478b-a2db-c304dea0d84a
fobj =  lfi.LabFbi.filtered_fluorescence(ffov, ofilter);

# ╔═╡ 34826b8a-6ff0-479a-80eb-0cd9ac5b9080
begin
	mdfobj = efluo_info_md(fobj, "Fluorescence filtered by objective")
	Markdown.parse("$mdfobj")
end

# ╔═╡ 19e03f4c-8620-4987-959c-f4977f19dd90
md"## Fluorescence as a function of concentration"

# ╔═╡ f96e9b12-9139-11eb-08cb-5319bca7a01d
md"## Topatu setup
- Laser
- Objective
- Filters
- CCD
"

# ╔═╡ 79af9066-8680-11eb-0d91-15759e995f6c
orca = lfi.LabFbi.ccd()

# ╔═╡ 2f259d3a-8688-11eb-1733-13978ae756ce
begin
	xL = collect(350.:5.:1000.)
	xEff =orca.(xL)
	porca = plot(xL, xEff, leg=false, lw=2);

	xlabel!("λ (nm)")
	ylabel!("Efficiency")
	title!("Efficiency of CCD ")
	
	
	naL = collect(0.1:0.05:1)
	obL = [lfi.LabFbi.Objective("Nikon-SLVD", na, 100) for na in naL]
	trL = [lfi.LabFbi.transmission(ob) for ob in obL]

	ptr = plot(naL, trL, leg=false, lw=2);

	xlabel!("numerical aperture")
	ylabel!("Transmission")
	title!("Tranmission of the objective ")
true
end

# ╔═╡ ee9fa653-3662-43d0-95af-331e61a12b92
plot(ptr,porca, layout = (1, 2), legend=false, fmt = :png)

# ╔═╡ 56d72e3c-86fa-11eb-2928-798e913a5976
md"### Filters"

# ╔═╡ 920bf9b9-7422-46cf-a772-1ee05e814474
begin
	λmin = 350.0nm
	λmax = 800.0nm
end

# ╔═╡ d9b5e50a-7125-4223-983f-dc2a2236abd8
filterd = lfi.LabFbi.get_filters(lfi.LabFbi.fnames, datadir("filters"))

# ╔═╡ cd53b04b-749a-47b3-93fe-0ad194c4f053
#Filters = lfi.LabFbi.load_filter.(lfi.LabFbi.fnames, (datadir("filters"),));

# ╔═╡ e2d92deb-8ac2-4991-a842-d9906316d6a8
#fdfs, fints = collect(zip(Filters...));

# ╔═╡ ae0c0e93-94eb-4249-a1a5-0e06071d2620
#lfi.LabFbi.fnames

# ╔═╡ fc629ff0-471e-46eb-8c5e-888f02801d8e
#fbp405, flp425, fbp430, flp450, fdn405_522, fn405 = fints

# ╔═╡ dd86edf9-9766-41d1-8e73-53cba6fe853d
fw = [390:420,420:500,420:520,440:540, 390:560, 390:450]

# ╔═╡ 981a0fb1-4a08-46c2-ba48-a95dd8b0f523
Pft = lpi.LabPlots.plot_filters(lfi.LabFbi.fnames, filterd, fw);

# ╔═╡ 81adc270-e6d8-46c5-a802-aaee7d78fc3b
plot(Pft..., layout = (3, 2), legend=false, fmt = :png)

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

# ╔═╡ 4af7c8cc-70fd-4b20-afde-acd6047ad75c
filterd

# ╔═╡ 5eaa4bc1-25d4-4a0e-837d-4c1f24f868c1
ftbp430 = lfi.LabFbi.filterBP430(filterd)

# ╔═╡ 2dfa4921-1560-4855-b936-d9b1fe8d93b2
ftlp450 = lfi.LabFbi.filterLP450(filterd)

# ╔═╡ 7e97d32a-7fb3-4722-a098-1f850cc0ed2d
cflt =lfi.LabFbi.common_filters(filterd)

# ╔═╡ c126f2d9-0ccb-4b5c-9600-e14ad96e9591
#filterSetBand430(λ) = map(x->flp425(x)^2 * fdn405(x) * fbp430(x) * orca(x), λ)

# ╔═╡ 8d4306e6-fe3e-43d6-b5cf-d5e74a3c509c
filterSetBand430(λ) = map(x->ftbp430(x) * orca(x), λ)

# ╔═╡ 7ccc9401-2fee-4d5e-a6ca-9a863c60d5a6
filterSetLP450(λ) = map(x->ftlp450(x) * orca(x), λ)

# ╔═╡ 789f7956-a8b9-4c1d-b07b-c57f48e51c8d
filterSetAll(λ) = map(x->cflt(x) * orca(x), λ)

# ╔═╡ e4abed36-f47b-4d77-9d10-e7455edd9325
pfilterSetBand430 =lpi.LabPlots.plot_filterset(filterSetBand430, "filterSetBand430", 400.0:500.0);

# ╔═╡ 202f831f-8e54-4d17-9f27-aea1a07ab919
pfilterSetLP450 =lpi.LabPlots.plot_filterset(filterSetLP450, "filterSetLP450", 400.0:700.0);

# ╔═╡ d05cd5b3-e5b5-4af5-b00a-3d2fe9eb7de3
pfilterSetAll =lpi.LabPlots.plot_filterset(filterSetAll, "filterSetALL", 400.0:700.0);

# ╔═╡ 33b05bd8-12a6-47bf-b229-f4b97ea279d0
plot(pfilterSetBand430, pfilterSetLP450, pfilterSetAll,layout = (1, 3), legend=false, fmt = :png)

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

# ╔═╡ 44911bef-7064-4a74-bb8d-7cf6fb66bac2
function dfhead(df, fbi=true)
	nx = split.(names(df),"_")
	nz = [n for n in nx[2:end]] 
	if fbi
		return nz[1][1], nz[1][3], [n[2] for n in nz]
	else
		return nz[1][1]*nz[1][2], nz[1][4], [n[3] for n in nz]
	end
end

# ╔═╡ 36cdec3a-b57e-428d-8cee-c232dcc55703
function dflm(df)
	l = df[!, "λ"]
	return l[1], l[end]
end

# ╔═╡ 67940e67-673f-4e64-901b-2a6b91d5a96b
function cnames(fbidf, nx)
	return names(fbidf)[nx]
end

# ╔═╡ 7a7f700b-f918-48ba-9085-f77b6751707b
struct FbiDfInfo
	name::String 
	nλ::String
	C::Array{String}
	CN::String
	λi::Float64
	λf::Float64
end

# ╔═╡ accf2ea3-2c2d-4166-8ca2-e618dcd0fdf3
function fbidfinfo(df, fbi=true,nx=2:6)
	fbili, fbilf, = dflm(df)
	fbiN, fbiL, fbiC = dfhead(df, fbi)
	fbiCs = lti.LabTools.vect_to_fstr(fbiC, "%s")
	Cnames = cnames(df, nx)
	return FbiDfInfo(fbiN,fbiL,Cnames,fbiCs,fbili,fbilf)
end
	

# ╔═╡ a4ce5e2e-6eb9-40f3-8cd3-7834e8020ede
function fbidfinfomd(dfi)
	md ="\n
- df name = $(dfi.name)
- df λ    = $(dfi.nλ)
- Concentration (in M) = $(dfi.CN)
- λmin = $(dfi.λi) nm, λmax = $(dfi.λf) nm
	"
	return md
end

# ╔═╡ f45d93f2-cc4d-4395-b9f3-f9a49983b8e2
begin
    wf = 350.
    we = 800.
    ws = 2.
    wl = wf:ws:we
end

# ╔═╡ 59d8cd90-6aba-4385-8409-c8a4c483a8a6
md"#### FbiG1Em405"

# ╔═╡ 158ab3e3-1650-4ab4-ae91-5ff8a73e42a3
fbi405df = lfi.LabFbi.load_df_from_csv(datadir("fluorimeter/fluorescence"), "FbiG1Em405.csv", lfi.LabFbi.spG);

# ╔═╡ bb3105cb-1af2-417f-8931-0a2426ddb622
begin
	fbig1405dfi = fbidfinfo(fbi405df)
	fbig1405dfimd = fbidfinfomd(fbig1405dfi)
	Markdown.parse("$fbig1405dfimd")
end

# ╔═╡ 517b220a-2dfe-4444-acfa-7c5b225de565
fbig1405dfi.CN

# ╔═╡ ce03c302-ec83-4a53-9c62-55584b08e7fb
function ncs(CN, ws=1:5)
	sC = lstrip.(split(CN,","))
	Cs =parse.((Float64,), sC)
	return (Cs / Cs[1])[ws]
end

# ╔═╡ fadf9e50-791d-4ac5-95bf-9ce72bc798d1
nCs = ncs(fbig1405dfi.CN, 1:5)

# ╔═╡ 22a8226c-aae4-489b-be73-77e1aba7eeda
Cs = fbig1405dfi.C

# ╔═╡ b1043565-10f9-4df0-b1b0-c031661c94ba
begin
	pfbig1405 = lpi.LabPlots.plotdf_xys(fbi405df, "λ", 
		                                fbig1405dfi.C, false, fbig1405dfi.C,
	                                    markercolors, 
					                    "λ(nm)", "I (au)", 
					                    "FBI G1", legend=:topright) ;
	pfbig1405N = lpi.LabPlots.plotdf_xys(fbi405df, "λ",   
		                                 fbig1405dfi.C, true, fbig1405dfi.C,
	                                     markercolors, 
					                     "λ(nm)", "I (au)", 
					                     "FBI G1", legend=:topright) ;
	true
end

# ╔═╡ 1906c991-cf1c-4c08-b0f4-8deea433cc08
plot(pfbig1405, pfbig1405N, layout = (1, 2), legend=:topright, fmt = :png)

# ╔═╡ c73995fc-accb-4936-ad9b-27be581c5009
md"#### FbiBaG1Em405"

# ╔═╡ de2195d2-d5c8-4979-8d8d-4233c1eb0195
fbiba405df = lfi.LabFbi.load_df_from_csv(datadir("fluorimeter/fluorescence"), "FbiBaG1Em405.csv", lfi.LabFbi.spG);

# ╔═╡ 42e3b3bc-c670-4f52-b0e6-826689183986
begin
	fbibag1405dfi = fbidfinfo(fbiba405df, false)
	fbibag1405dfimd = fbidfinfomd(fbibag1405dfi)
	Markdown.parse("$fbibag1405dfimd")
end

# ╔═╡ d4837bf0-a0ca-455f-9d27-f67634db5a11
begin
	pfbibag1405 = lpi.LabPlots.plotdf_xys(fbiba405df, "λ", 
		                                  fbibag1405dfi.C, false, fbibag1405dfi.C,
	                					  markercolors, 
										  "λ(nm)", "I (au)", 
										  "FBI G1", legend=:topright) ;
	pfbibag1405N = lpi.LabPlots.plotdf_xys(fbiba405df, "λ", 
										   fbibag1405dfi.C, true, fbibag1405dfi.C,
										   markercolors, 
										   "λ(nm)", "I (au)", 
										   "FBI G1", legend=:topright) ;
	true
end

# ╔═╡ 80b1d5de-18c8-4d38-8c75-deb015ea5f66
plot(pfbibag1405, pfbibag1405N, layout = (1, 2), legend=:topright, fmt = :png)

# ╔═╡ bd916cbd-10c1-4e34-97f4-147e805f9024
md"#### FbiG1Em325"

# ╔═╡ 630e87d8-7d15-4699-8340-33c5485e5145
fbi325df = lfi.LabFbi.load_df_from_csv(datadir("fluorimeter/fluorescence"), "FbiG1Em325.csv", lfi.LabFbi.spG);

# ╔═╡ baa16010-727c-4aba-b965-f5cf310cb15e
begin
	fbig1325dfi = fbidfinfo(fbi325df)
	fbig1325dfimd = fbidfinfomd(fbig1325dfi)
	Markdown.parse("$fbig1325dfimd")
end

# ╔═╡ e1c1853e-b2b8-4df3-bbf3-6d33e6d0ed25
begin
	pfbig1325 = lpi.LabPlots.plotdf_xys(fbi325df, "λ", fbig1325dfi.C, 
										false, fbig1325dfi.C,
										markercolors, 
										"λ(nm)", "I (au)", 
										"FBI G1 (325 nm)", legend=:topright) ;
	pfbig1325N = lpi.LabPlots.plotdf_xys(fbi325df, "λ", fbig1325dfi.C, 
										 true, fbig1325dfi.C,
										 markercolors, 
										 "λ(nm)", "I (au)", 
										 "FBI G1 (325 nm)", legend=:topright) ;
	true
end

# ╔═╡ 98433f26-dab9-4a37-938f-c37b7dfc99ec
plot(pfbig1325, pfbig1405N, layout = (1, 2), legend=:topright, fmt = :png)

# ╔═╡ 5e5a4301-f746-4cb7-9802-9891feee2562
md"#### FbiBaG1Em325"

# ╔═╡ 793a7a8a-bd2f-4bf1-b172-1f2913124581
fbiBa325df = lfi.LabFbi.load_df_from_csv(datadir("fluorimeter/fluorescence"), "FbiBaG1Em325.csv", lfi.LabFbi.spG);

# ╔═╡ df42a91b-05f7-4222-9fc0-d9cda9a85f77
begin
	fbibag1325dfi = fbidfinfo(fbiBa325df, false)
	fbibag1325dfimd = fbidfinfomd(fbibag1325dfi)
	Markdown.parse("$fbibag1325dfimd")
end

# ╔═╡ 0efa00f4-cde3-421c-877b-a08edcb391ad
begin
	pfbibag1325 = lpi.LabPlots.plotdf_xys(fbiBa325df, "λ", 
										  fbibag1325dfi.C, false, fbibag1325dfi.C,
										  markercolors, 
										  "λ(nm)", "I (au)", 
										  "FBI G1", legend=:topright) ;
	pfbibag1325N = lpi.LabPlots.plotdf_xys(fbiBa325df, "λ", 
										   fbibag1325dfi.C, true, fbibag1325dfi.C,
										   markercolors, 
										   "λ(nm)", "I (au)", 
										   "FBI G1", legend=:topright) ;
	true
end

# ╔═╡ e7b7deaa-ae96-444e-ab46-eaf8b23ba72c
plot(pfbibag1325, pfbibag1325N, layout = (1, 2), legend=:topright, fmt = :png)

# ╔═╡ 5021b4e5-8590-4f1c-b647-a0f84f4c2fba
md"### Effect of concentration in FBI-G1

- Intensity decreases with concentration, but the response is not linear at high concentration.

- Molecule only chelated for high concentration (association constant is not large enough to chelate for smaller concentrations. 
"

# ╔═╡ 154a0653-3e59-462e-a9c7-6aba75c6e2c1
md"#### FbiG2Em325"

# ╔═╡ 708120f0-8962-4bc4-aac4-f025af8f9374
fbiG2325df = lfi.LabFbi.load_df_from_csv(datadir("fluorimeter/fluorescence"), "FbiG2Em325.csv", lfi.LabFbi.spG);

# ╔═╡ ca9eee65-84e0-4a0d-a4ab-1bafbc8d5f77
begin
	fbig2325dfi = fbidfinfo(fbiG2325df)
	fbig2325dfimd = fbidfinfomd(fbig2325dfi)
	Markdown.parse("$fbig2325dfimd")
end

# ╔═╡ 54ab8f4f-e40c-4479-9fe4-e79490a6b60d
begin
	pfbig2325 = lpi.LabPlots.plotdf_xys(fbiG2325df, "λ", 
										fbig2325dfi.C, false, fbig2325dfi.C,
										markercolors, 
										"λ(nm)", "I (au)", 
										"FBI G2 (325 nm)", legend=:topright) ;
	pfbig2325N = lpi.LabPlots.plotdf_xys(fbiG2325df, "λ", 
										 fbig2325dfi.C, true, fbig2325dfi.C,
										markercolors, 
										"λ(nm)", "I (au)", 
										"FBI G2 (325 nm)", legend=:topright) ;
	true
end

# ╔═╡ bd0f5ac3-2b57-4ea7-b2b5-d77b515b93b3
plot(pfbig2325, pfbig2325N, layout = (1, 2), legend=false, fmt = :png)

# ╔═╡ acc09228-daca-47bb-8d3e-bdf7d52adbcf
md"#### FbiBaG2Em325"

# ╔═╡ 15bc7ce6-2066-4ac4-96ad-6297dfbaa126
fbiBaG2325df = lfi.LabFbi.load_df_from_csv(datadir("fluorimeter/fluorescence"), "FbiBaG2Em325.csv", lfi.LabFbi.spG);

# ╔═╡ acf23156-8da7-40f0-96b5-e5bfc32759b8
begin
	fbibag2325dfi = fbidfinfo(fbiBaG2325df, false)
	fbibag2325dfimd = fbidfinfomd(fbibag2325dfi)
	Markdown.parse("$fbibag2325dfimd")
end

# ╔═╡ 593a1ba2-ec68-4679-adef-b3dee12e8a08
pfbibag2325 = lpi.LabPlots.plotdf_xys(fbiBaG2325df, "λ", 
									 fbibag2325dfi.C, false, fbibag2325dfi.C,
									 markercolors, 
									 "λ(nm)", "I (au)", 
									 "FBI Ba G2", legend=:topright);

# ╔═╡ bc0b04f1-32c4-4e0b-8bd9-f9f1b27dfd7d
pfbibag2325N = lpi.LabPlots.plotdf_xys(fbiBaG2325df, "λ", fbibag2325dfi.C, 
									  true, fbibag2325dfi.C,
									  markercolors, 
									  "λ(nm)", "I (au)", 
									  "FBI Ba G2", legend=:topright);

# ╔═╡ 34f3a8a3-ff2b-4c94-a4d2-50b6bd0b1376
plot(pfbibag2325, pfbibag2325N, layout = (1, 2), legend=false, fmt = :png)

# ╔═╡ eb41faf1-9ca1-4c90-9e74-23c5c6d5494a
plot(pfbig2325N, pfbibag2325N, layout = (1, 2), legend=false, fmt = :png)

# ╔═╡ 141a14d6-08dd-4263-82ff-1895f4f85713
pfbig2n = lpi.LabPlots.merge_plots!(pfbig2325N, pfbibag2325N);

# ╔═╡ d130fc32-6283-4cba-9958-05f7641e4753
pfbig2 = lpi.LabPlots.merge_plots!(pfbig2325, pfbibag2325);

# ╔═╡ eb3056e0-42b7-4981-bf66-c636e8e72d37
plot(pfbig2, pfbig2n, layout = (1, 2), legend=false, fmt = :png)

# ╔═╡ 9ce59fce-7b74-419e-a5f5-06684689792a
md"### PDF functions"

# ╔═╡ 68fb50e3-e690-4a59-88c5-c1973a109c00
begin
	gfbi405 = lfi.LabFbi.dftogf.((fbig1405dfi.λi:2.0:fbig1405dfi.λf,), (fbi405df,), fbig1405dfi.C)
	gfbiba405 = lfi.LabFbi.dftogf.((fbibag1405dfi.λi:2.0:fbibag1405dfi.λf,), (fbiba405df,), fbibag1405dfi.C)
end

# ╔═╡ b573b01c-e0b4-47cf-beb3-e6bf8c5dc453
begin
	pgfbi405 =lpi.LabPlots.plotdf_gfs(gfbi405, wl,fbig1405dfi.C, markercolors, 
	                    "λ (nm)", "I (au)", "FBIG1",
		                legend=:topright);
	pgfbiba405 =lpi.LabPlots.plotdf_gfs(gfbiba405, wl,fbibag1405dfi.C, markercolors, 
	                    "λ (nm)", "I (au)", "FBIBaG1",
		                legend=:topright);
true
end

# ╔═╡ 44db087f-95ae-4744-bead-580f06c667f4
plot(pgfbi405, pgfbiba405, layout = (1, 2), legend=false, fmt = :png)

# ╔═╡ e8a88d6b-bd09-4426-b24d-a56719a82b3c
begin
	gfbi325 = lfi.LabFbi.dftogf.((fbig1325dfi.λi:2.0:fbig1325dfi.λf,), (fbi325df,), fbig1325dfi.C)
	gfbiba325 = lfi.LabFbi.dftogf.((fbibag1325dfi.λi:2.0:fbibag1325dfi.λf,), (fbiBa325df,), fbibag1325dfi.C)
end

# ╔═╡ cd430d52-c37e-471a-a4d3-1acc328fb59b
begin
	pgfbi325 =lpi.LabPlots.plotdf_gfs(gfbi325, wl,fbig1325dfi.C, markercolors, 
	                    "λ (nm)", "I (au)", "FBIG1",
		                legend=:topright);
	pgfbiba325 =lpi.LabPlots.plotdf_gfs(gfbiba325, wl,fbibag1325dfi.C, markercolors, 
	                    "λ (nm)", "I (au)", "FBIBaG1",
		                legend=:topright);
true
end

# ╔═╡ 9b3376e4-b675-483e-afd1-b61f0c00d03c
plot(pgfbi325, pgfbiba325, layout = (1, 2), legend=false, fmt = :png)

# ╔═╡ 98c6064f-8389-498e-a034-728ba4d8a683
begin
	gfbi325g2 = lfi.LabFbi.dftogf.((fbig2325dfi.λi:2.0:fbig2325dfi.λf,), (fbiG2325df,), fbig2325dfi.C)
	gfbiba325g2 = lfi.LabFbi.dftogf.((fbibag2325dfi.λi:2.0:fbibag2325dfi.λf,), (fbiBaG2325df,), fbibag2325dfi.C)
end

# ╔═╡ 783aeb63-5783-4cf4-a2ec-6f5c9f1f9d7f
begin
	pgfbi325g2 =lpi.LabPlots.plotdf_gfs(gfbi325g2, wl,fbig2325dfi.C, markercolors, 
	                    "λ (nm)", "I (au)", "FBIG2",
		                legend=:topright);
	pgfbiba325g2 =lpi.LabPlots.plotdf_gfs(gfbiba325g2, wl,fbibag2325dfi.C, markercolors, 
	                    "λ (nm)", "I (au)", "FBIBaG2",
		                legend=:topright);
true
end

# ╔═╡ edee404b-39fb-4ff8-8bee-0dfc25d1f6e1
plot(pgfbi325g2, pgfbiba325g2, layout = (1, 2), legend=false, fmt = :png)

# ╔═╡ 1c61d0e8-e197-4fd0-9aa1-7b0728179209
lpi.LabPlots.merge_plots!(pgfbi325g2, pgfbiba325g2)

# ╔═╡ 66523f06-fe31-4642-9e46-0ecb1ed0fbf4
md"""
### Filtering
"""

# ╔═╡ b161ba71-3625-4d72-bd5d-a58446327d8e
begin
	gFbi   = Dict{String, Any}()
	gFbiBa = Dict{String, Any}()
	gFbi["l405g1"]= [g.pdf for g in gfbi405]
	gFbiBa["l405g1"]= [g.pdf for g in gfbiba405]
	gFbi["l325g1"]= [g.pdf for g in gfbi325]
	gFbiBa["l325g1"]= [g.pdf for g in gfbiba325]
	gFbi["l325g2"]= [g.pdf for g in gfbi325g2]
	gFbiBa["l325g2"]= [g.pdf for g in gfbiba325g2]
end

# ╔═╡ 8412be2d-bd07-4ae3-82f1-b35ab5949f3e
struct GBTransport
	f430BP::Function 
	f450LP::Function 
end

# ╔═╡ 4ead6795-7944-4445-b83f-1a4e7d8ab70f
function gb_transport(gf=gfbi["l405g1"], fluo=fobj.fbi["l405g1"]/Hz, nC = nCs)
	f430BP(λ) = map(x->filterSetBand430(x) * gf(x) * fluo * nC, λ)
	f450LP(λ)  = map(x->filterSetLP450(x) * gf(x) * fluo * nC,  λ)
	return GBTransport(f430BP,f450LP)
end

# ╔═╡ 7b6046b7-473b-4809-a8e5-7cfc9dd5bf11
fbiG1gbt = gb_transport.(gFbi["l405g1"], fobj.fbi["l405g1"]/Hz, nCs)

# ╔═╡ 5fa23bef-9e69-45b5-8ed2-e3206d97c92e
fbibaG1gbt = gb_transport.(gFbiBa["l405g1"], fobj.fbi["l405g1"]/Hz, nCs)

# ╔═╡ a264a212-f6a3-41ac-a80b-7567ded2af0d


# ╔═╡ 368c4e24-f367-4383-a64d-c104ef346776
function filtered_fluorescence(fbiGbt=fbiG1gbt, 
		                       fbiBaGbt=fbibaG1gbt, 
		                       wr=400.0:700.0,
		                       wbp=(420.0,440.0), 
		                       wlp=(450.0,700.0))
	
	pFbi430BP   = []
	pFbi450LP   = []
	pFbiBa430BP = []
	pFbiBa450LP = []
	nFbi430BP   = []
	nFbi450LP   = []
	nFbiBa430BP = []
	nFbiBa450LP = []

	for i in [1,3,4,5]
		push!(pFbi430BP, lpi.LabPlots.plot_filterset(fbiGbt[i].f430BP, 
												"filteredFbi430BP", wr))
			
		push!(nFbi430BP, lfi.LabFbi.qpdf(fbiGbt[i].f430BP, wbp[1],wbp[2]))
			
		push!(pFbi450LP, lpi.LabPlots.plot_filterset(fbiGbt[i].f450LP, 
												"filteredFbi450LP", wr))
			
		push!(nFbi450LP, lfi.LabFbi.qpdf(fbiGbt[i].f450LP,  wlp[1],wlp[2]))
			
		push!(pFbiBa430BP, lpi.LabPlots.plot_filterset(fbiBaGbt[i].f430BP, 
												"filteredFbiBa430BP", wr))	
			
		push!(nFbiBa430BP,lfi.LabFbi.qpdf(fbiBaGbt[i].f430BP, wbp[1],wbp[2]))
		
		push!(pFbiBa450LP, lpi.LabPlots.plot_filterset(fbiBaGbt[i].f450LP, 
												"filteredFbiba450LP", wr))
			
		push!(nFbiBa450LP,lfi.LabFbi.qpdf(fbiBaGbt[i].f450LP, wlp[1],wlp[2]))
	end
	pffFbi = Dict("pFbi430BP" => pFbi430BP, "pFbi450LP" =>pFbi450LP,
				  "pFbiBa430BP" => pFbiBa430BP, "pFbiBa450LP" =>pFbiBa450LP)
	nffFbi = Dict("nFbi430BP" => nFbi430BP, "nFbi450LP" =>nFbi450LP,
				  "nFbiBa430BP" => nFbiBa430BP, "nFbiBa450LP" =>nFbiBa450LP)
	return 	nffFbi, pffFbi			
end

# ╔═╡ 4d8b56d4-1ac2-4e37-8f2b-e59173514fbd
nffFbi, pffFbi = filtered_fluorescence()

# ╔═╡ b05bee10-c034-4e92-8233-20d388aaf9ca
plot(pffFbi["pFbi430BP"]..., layout = (2, 2), legend=false, fmt = :png)

# ╔═╡ ac2da991-20fa-4922-99c0-e28d4f395300
plot(pffFbi["pFbi450LP"]..., layout = (2, 2), legend=false, fmt = :png)

# ╔═╡ e2de29a2-61ac-412e-beb3-bf019d437214
plot(pffFbi["pFbiBa430BP"]..., layout = (2, 2), legend=false, fmt = :png)

# ╔═╡ 8933349f-91c0-4329-b3e1-7aac01ffb557
plot(pffFbi["pFbiBa450LP"]..., layout = (2, 2), legend=false, fmt = :png)

# ╔═╡ edb845bc-7e73-45fa-bd88-de9f38559407
md"""
### On the shape of FBI

- Except for very large concentrations, the normalised pdfs for the free molecule is the same for all the concentrations.
- In the case of G1, we only observe chelation for very high concentration, and thus the response depends strongly from the concentration.
"""

# ╔═╡ 562ebf36-19a7-48af-8825-cbab5ce7e3b7
md"## Predicted rates in CCD"

# ╔═╡ 29fe8937-216b-4c30-9808-dc95d6208b1f
function rfluo_info_md(nffFbi, Cs, i=1)
	
	sffFbi = Dict("Fbi430BP" => lti.LabTools.to_fstr(nffFbi["nFbi430BP"][i], 
			                    "%5.2g"),
		          "Fbi450LP" =>lti.LabTools.to_fstr(nffFbi["nFbi450LP"][i],
			                    "%5.2g"),
				  "FbiBa430BP" => lti.LabTools.to_fstr(nffFbi["nFbiBa430BP"][i],
			                   "%5.2g"),
		          "FbiBa450LP" =>lti.LabTools.to_fstr(nffFbi["nFbiBa450LP"][i],
			                    "%5.2g"))


	md = """ ### Fluorescence recorded by the CCD: $(Cs[i]) \n
	
| Filter | Fbi G1 (Hz) | FbiBa G1 (Hz)
|:---------- | ---------- |:------------:|
| BP430    | $(sffFbi["Fbi430BP"])| $(sffFbi["FbiBa430BP"]) | 
| LP450    | $(sffFbi["Fbi450LP"])| $(sffFbi["FbiBa450LP"])   |

"""

	return md
end

# ╔═╡ 63a54c1b-3b55-4b79-aa61-8ffa75c0429c
begin
	rfmdc1 = rfluo_info_md(nffFbi, Cs, 1)
	Markdown.parse("$rfmdc1")
end

# ╔═╡ 639ca6d6-e75f-41b8-870c-bc736c0f758b
begin
	rfmdc2 = rfluo_info_md(nffFbi, Cs, 2)
	Markdown.parse("$rfmdc2")
end

# ╔═╡ c4268ea2-fdfd-47ff-bf99-206f67cdaec7
begin
	rfmdc3 = rfluo_info_md(nffFbi, Cs, 3)
	Markdown.parse("$rfmdc3")
end

# ╔═╡ f71fac30-37a2-49d8-a12b-94fd71534d57
begin
	rfmdc4 = rfluo_info_md(nffFbi, Cs, 4)
	Markdown.parse("$rfmdc4")
end

# ╔═╡ db3c7664-a532-466f-b65a-15800d5c6363
md"""## The case of Monolayers

- We do not know what is the molar extinction coefficient (cross section) in dry medium. 
- We don't know either if the shape in a mono-layer is the same than in solution.
- Instead, knowing the efficiency of the TOPATU setup, we can measure those cross sections.

The number of photons produced per fluorophore is: 
$ \gamma = I \cdot \sigma$, 

where $I$ is the beam density in the spot and $\sigma$ the effective cross section, already multiplied by the quantum efficiency of the fluorophore.

The number of molecules contained in the spot (e.g, at the diffractive limit) depends on one hand of the layer density (roughly: $10^6$ molecules per layer) and how many of these molecules are excited by the laser. In the case of the unchelated species, one expects that all molecules are excited by the laser. Thus if $m$ is the number of molecules in the spot, the total fluorescence emitted is:

$f = \gamma \cdot m = I \cdot \sigma \cdot m$

The whole setup (including filters, objective and CCD) has an optical efficiency of $\epsilon$. Thus, the recorded fluorescence, $F$, is:

$F = f \cdot \epsilon = I \cdot \sigma \cdot m \cdot \epsilon$

And from here we can find the cross section:

$\sigma = \frac{F}{I\cdot m \cdot \epsilon}$

"""

# ╔═╡ da36394d-a8b1-4914-99e4-826d10b1a168
ml = lfi.LabFbi.Mlayer("FBI",1.0/nm^2)

# ╔═╡ 6b73fcb5-535d-498e-b500-abbc27b98173
nml = lfi.LabFbi.nofa(ml, lT.fovs["l405"].a)

# ╔═╡ 07932901-f0c4-4edc-9fc3-0c8d7af25282
md"- emitted fluorescence (assuming same cross section than in solution)"

# ╔═╡ e1787274-9153-47e9-8ba0-7d68283a4dfc
efluo.fbi["l405g1"]

# ╔═╡ 85737ce4-50d2-40fa-955c-ebf3488ca8bf
f = nml * efluo.fbi["l405g1"]

# ╔═╡ 4095cee2-f2f0-4a59-aa24-1f520f8c1360
md"- Fluorescence transmitted by the objective"

# ╔═╡ cf7d64e7-4af8-47e9-af48-d58595afec61
ffm = ofilter * f

# ╔═╡ d9847db0-c937-45df-b3c1-0e976898b091
ofilter

# ╔═╡ 4c9befae-555d-4a51-a372-225085b90dc8
f

# ╔═╡ beefb08e-e8f2-4e65-b6f6-1562ba63847d
fbix = gFbi["l405g1"][2]

# ╔═╡ 3fb7a6b2-e55a-47b4-8578-f2b32071ca2f
FBIX(λ) = map(x->filterSetAll(x) * fbix(x), λ)

# ╔═╡ 3370d88d-b2e2-446f-b8ee-16e4b9c1dbec
FBIX(700.)

# ╔═╡ 54dd1897-faae-4920-aa12-712996254385
lpi.LabPlots.plot_filterset(FBIX, "filteredFbi", 400.0:700.0)

# ╔═╡ c91a0c19-1c98-424f-a37b-049083db3acb
tflt = lfi.LabFbi.qpdf(FBIX, 400.0, 700.0)

# ╔═╡ 5d467edc-74de-4b8e-a13b-d104c052b8ac
fluoccd = ofilter * f * tflt

# ╔═╡ 0ff2061a-7f07-4f58-ba34-4a9536bafed5
dr = 13μm

# ╔═╡ 0edd6a6d-41b4-4523-84c4-5667fd9bf126
ar = π*dr^2

# ╔═╡ 7106c47a-d6c2-47b5-a899-347268faf94c
ars = lti.LabTools.to_fstr(ar/μm^2,"%5.2g")

# ╔═╡ 938eb752-fbed-4545-bb7b-4e3b6821b946
md"- This number needs to be normalised to the ROI area taken in the CCD to integrate the charge. The radius of the ROI is $dr and thus, the area is $ars μm2"

# ╔═╡ ee5a9a7c-38f7-4b3c-9591-8a34cf772fec
io = fluoccd / ar

# ╔═╡ 06930724-252f-42f9-8c73-7fcfcbed9b76
ios = lti.LabTools.to_fstr(io/(Hz * μm^-2),"%5.2g");

# ╔═╡ c2434331-c6ee-45bf-8d9e-049c1a827bd2
md"- Therefore, the observed irradiance would be: $ios Hz/μm2"

# ╔═╡ Cell order:
# ╟─79cfd2fc-9046-11eb-2b13-1b877d57645d
# ╟─4f35dcfc-8644-11eb-1c86-75608b38bf5f
# ╟─3207f446-8643-11eb-37ba-c9aec47fcb8f
# ╟─5115917a-8644-11eb-19fc-0528741ca75d
# ╟─fc8b79a2-8728-11eb-2da7-e3ffa3ceef08
# ╟─68e738e7-88bd-41c2-89e4-594f07d64ddc
# ╟─0f2f4c78-8729-11eb-2bab-27812ce8c47e
# ╟─5ee27d52-86fd-11eb-365e-9f2b1a095575
# ╟─621ec96c-86fd-11eb-1c41-379cc17180dc
# ╟─9b853f27-4288-42a6-8f12-ca004e1773b7
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
# ╟─c0fd3834-893c-45d5-8530-8758c8661e1e
# ╟─990240b5-b30d-4c50-911c-42f01d7c98c3
# ╟─0519d9f0-a35a-4f5c-a70b-692809bee03e
# ╠═e9076d72-8648-11eb-0c02-9b00b0605e72
# ╠═dad410f8-cca9-4f2e-8735-b2e1fb135ca8
# ╟─052941b7-6a8e-4930-bce7-7deec4c13ff0
# ╟─9f9b7846-2120-4c8e-bb45-61e00115d71b
# ╟─c5f5887d-157c-42cb-a0ed-e3f82b553eec
# ╟─bb4fa686-908b-11eb-0aa0-27ab50bb0a4e
# ╠═1c31637a-8666-11eb-2ea6-17b0d4507e10
# ╟─e0bbf162-fea4-427b-9133-569b51632746
# ╠═28c0a218-8666-11eb-291a-71fd93c60da5
# ╟─943550b1-7a39-4e99-9054-4b8e1c584b6f
# ╠═0bddff38-99ad-475e-9515-579de07e8ab3
# ╟─375afbf8-47bd-470d-93f2-715f13760501
# ╟─5bd97b62-5d4d-4483-8c05-2eec1d1e6b57
# ╟─f4fd8031-1b85-4f98-898f-eb09835ddc2d
# ╟─62da60bf-df76-41f9-9027-78953fbbcdb8
# ╟─30371eba-26ac-429a-bb2e-30de65091288
# ╟─b78e625e-4c71-4f1c-8711-cfe90c27519b
# ╟─dd35eedb-9ed0-4bc2-9636-3192c1b034ab
# ╟─3d182094-a3bb-413a-a1fc-6026839557c0
# ╟─f0b9c27a-9623-48ee-b096-6a3fb56bc35c
# ╟─013420a8-9129-11eb-0d7f-ab8fb7c70786
# ╠═fe18d00d-3fe9-4bf4-b926-2dcc12d513ab
# ╠═55046ba7-f5a6-4d9a-a767-3dcdb85bb1dc
# ╠═41e148a6-5dec-49a2-bffa-c4b9d5816d77
# ╟─a80e6642-908c-11eb-3ac3-6737f11f4eff
# ╟─9cb2c012-608a-4cff-b312-88bd2071da1f
# ╠═cd707e24-df49-4f74-9890-0cdfcc1b9296
# ╟─482b76cf-b3fb-4a30-8678-56ebe539fc5e
# ╠═fe77ba58-c130-4870-a7c0-180bd9af22bb
# ╟─2861c9e2-ae5a-4a9b-91ad-88a1ea51b450
# ╟─1e9740d1-1fad-4a0c-bfd4-a189857bd9f0
# ╟─cfe3217e-57b3-4cdb-ba09-5df3f194ceac
# ╠═ccb649f8-0516-4a64-8ccc-5af13a7fe6e1
# ╟─a6db60d8-371d-4b93-b164-38d44f00d17a
# ╟─96639a5e-81c2-4eac-837c-fff5bc7ea5aa
# ╟─639922e2-c0e5-48ee-a2b2-9e3afd543e1e
# ╟─205870db-bf67-4ece-8de4-239aaa5edc47
# ╠═21c6d2f8-519e-41ff-8895-ae3c77eb7a02
# ╟─61de7676-2d4f-40da-9250-92d00f8d780b
# ╟─2acdf979-7a0a-4538-92bd-f31df54c503e
# ╠═5347f077-66b4-45e6-a2f1-b0e6f52b0768
# ╟─8add778e-bffd-4934-a867-6f9658b4be89
# ╠═7d2e1f10-6a8b-46af-b2f9-9438ff3f2f89
# ╠═32503979-9ef1-4b7a-aeb9-1a63b33f60cd
# ╠═1aefe05a-e218-473f-a9df-56eb3fd844cb
# ╟─65833e10-898f-44a9-8731-3e3df800a2f6
# ╟─38979d58-a88f-4e3c-ba0e-e0ed7675054d
# ╠═8e16733d-9b6f-4170-a05f-e23cc63b4a05
# ╠═b938a16d-3a15-478b-a2db-c304dea0d84a
# ╟─34826b8a-6ff0-479a-80eb-0cd9ac5b9080
# ╠═19e03f4c-8620-4987-959c-f4977f19dd90
# ╟─f96e9b12-9139-11eb-08cb-5319bca7a01d
# ╠═79af9066-8680-11eb-0d91-15759e995f6c
# ╟─2f259d3a-8688-11eb-1733-13978ae756ce
# ╠═ee9fa653-3662-43d0-95af-331e61a12b92
# ╟─56d72e3c-86fa-11eb-2928-798e913a5976
# ╠═920bf9b9-7422-46cf-a772-1ee05e814474
# ╠═d9b5e50a-7125-4223-983f-dc2a2236abd8
# ╠═cd53b04b-749a-47b3-93fe-0ad194c4f053
# ╠═e2d92deb-8ac2-4991-a842-d9906316d6a8
# ╠═ae0c0e93-94eb-4249-a1a5-0e06071d2620
# ╠═fc629ff0-471e-46eb-8c5e-888f02801d8e
# ╠═dd86edf9-9766-41d1-8e73-53cba6fe853d
# ╠═981a0fb1-4a08-46c2-ba48-a95dd8b0f523
# ╠═81adc270-e6d8-46c5-a802-aaee7d78fc3b
# ╟─b0b10062-93c9-11eb-1b87-0f76a1735fa1
# ╠═4af7c8cc-70fd-4b20-afde-acd6047ad75c
# ╠═5eaa4bc1-25d4-4a0e-837d-4c1f24f868c1
# ╠═2dfa4921-1560-4855-b936-d9b1fe8d93b2
# ╠═7e97d32a-7fb3-4722-a098-1f850cc0ed2d
# ╠═c126f2d9-0ccb-4b5c-9600-e14ad96e9591
# ╠═8d4306e6-fe3e-43d6-b5cf-d5e74a3c509c
# ╠═7ccc9401-2fee-4d5e-a6ca-9a863c60d5a6
# ╠═789f7956-a8b9-4c1d-b07b-c57f48e51c8d
# ╠═e4abed36-f47b-4d77-9d10-e7455edd9325
# ╠═202f831f-8e54-4d17-9f27-aea1a07ab919
# ╠═d05cd5b3-e5b5-4af5-b00a-3d2fe9eb7de3
# ╠═33b05bd8-12a6-47bf-b229-f4b97ea279d0
# ╟─92196df1-1448-41aa-afd1-0c4de5dc0dc8
# ╟─78f98759-f580-4d42-85d1-fde4662b666c
# ╟─b8fe9096-adfd-46f9-84ce-9be60ec49eff
# ╠═f273fd6a-e1eb-4d80-86ab-ecdb43be3e53
# ╠═df379f4e-c27b-422d-a1f3-1373cd1f7e76
# ╠═3d982cfb-7c5e-422a-8120-b6048b3cd658
# ╟─da582c1b-0f85-48cf-afd1-61e8a83e6c48
# ╠═241a2416-6bc7-4804-b6ae-874d9336a02a
# ╠═44911bef-7064-4a74-bb8d-7cf6fb66bac2
# ╠═36cdec3a-b57e-428d-8cee-c232dcc55703
# ╠═67940e67-673f-4e64-901b-2a6b91d5a96b
# ╠═7a7f700b-f918-48ba-9085-f77b6751707b
# ╠═accf2ea3-2c2d-4166-8ca2-e618dcd0fdf3
# ╠═a4ce5e2e-6eb9-40f3-8cd3-7834e8020ede
# ╠═f45d93f2-cc4d-4395-b9f3-f9a49983b8e2
# ╠═59d8cd90-6aba-4385-8409-c8a4c483a8a6
# ╠═158ab3e3-1650-4ab4-ae91-5ff8a73e42a3
# ╠═bb3105cb-1af2-417f-8931-0a2426ddb622
# ╠═517b220a-2dfe-4444-acfa-7c5b225de565
# ╠═ce03c302-ec83-4a53-9c62-55584b08e7fb
# ╠═fadf9e50-791d-4ac5-95bf-9ce72bc798d1
# ╠═22a8226c-aae4-489b-be73-77e1aba7eeda
# ╠═b1043565-10f9-4df0-b1b0-c031661c94ba
# ╠═1906c991-cf1c-4c08-b0f4-8deea433cc08
# ╠═c73995fc-accb-4936-ad9b-27be581c5009
# ╠═de2195d2-d5c8-4979-8d8d-4233c1eb0195
# ╠═42e3b3bc-c670-4f52-b0e6-826689183986
# ╠═d4837bf0-a0ca-455f-9d27-f67634db5a11
# ╠═80b1d5de-18c8-4d38-8c75-deb015ea5f66
# ╠═bd916cbd-10c1-4e34-97f4-147e805f9024
# ╠═630e87d8-7d15-4699-8340-33c5485e5145
# ╠═baa16010-727c-4aba-b965-f5cf310cb15e
# ╠═e1c1853e-b2b8-4df3-bbf3-6d33e6d0ed25
# ╠═98433f26-dab9-4a37-938f-c37b7dfc99ec
# ╠═5e5a4301-f746-4cb7-9802-9891feee2562
# ╠═793a7a8a-bd2f-4bf1-b172-1f2913124581
# ╠═df42a91b-05f7-4222-9fc0-d9cda9a85f77
# ╠═0efa00f4-cde3-421c-877b-a08edcb391ad
# ╠═e7b7deaa-ae96-444e-ab46-eaf8b23ba72c
# ╠═5021b4e5-8590-4f1c-b647-a0f84f4c2fba
# ╠═154a0653-3e59-462e-a9c7-6aba75c6e2c1
# ╠═708120f0-8962-4bc4-aac4-f025af8f9374
# ╠═ca9eee65-84e0-4a0d-a4ab-1bafbc8d5f77
# ╟─54ab8f4f-e40c-4479-9fe4-e79490a6b60d
# ╠═bd0f5ac3-2b57-4ea7-b2b5-d77b515b93b3
# ╠═acc09228-daca-47bb-8d3e-bdf7d52adbcf
# ╠═15bc7ce6-2066-4ac4-96ad-6297dfbaa126
# ╠═acf23156-8da7-40f0-96b5-e5bfc32759b8
# ╟─593a1ba2-ec68-4679-adef-b3dee12e8a08
# ╟─bc0b04f1-32c4-4e0b-8bd9-f9f1b27dfd7d
# ╠═34f3a8a3-ff2b-4c94-a4d2-50b6bd0b1376
# ╠═eb41faf1-9ca1-4c90-9e74-23c5c6d5494a
# ╠═141a14d6-08dd-4263-82ff-1895f4f85713
# ╠═d130fc32-6283-4cba-9958-05f7641e4753
# ╠═eb3056e0-42b7-4981-bf66-c636e8e72d37
# ╠═9ce59fce-7b74-419e-a5f5-06684689792a
# ╠═68fb50e3-e690-4a59-88c5-c1973a109c00
# ╠═b573b01c-e0b4-47cf-beb3-e6bf8c5dc453
# ╠═44db087f-95ae-4744-bead-580f06c667f4
# ╟─e8a88d6b-bd09-4426-b24d-a56719a82b3c
# ╟─cd430d52-c37e-471a-a4d3-1acc328fb59b
# ╠═9b3376e4-b675-483e-afd1-b61f0c00d03c
# ╟─98c6064f-8389-498e-a034-728ba4d8a683
# ╟─783aeb63-5783-4cf4-a2ec-6f5c9f1f9d7f
# ╠═edee404b-39fb-4ff8-8bee-0dfc25d1f6e1
# ╠═1c61d0e8-e197-4fd0-9aa1-7b0728179209
# ╠═66523f06-fe31-4642-9e46-0ecb1ed0fbf4
# ╠═b161ba71-3625-4d72-bd5d-a58446327d8e
# ╠═8412be2d-bd07-4ae3-82f1-b35ab5949f3e
# ╠═4ead6795-7944-4445-b83f-1a4e7d8ab70f
# ╠═7b6046b7-473b-4809-a8e5-7cfc9dd5bf11
# ╠═5fa23bef-9e69-45b5-8ed2-e3206d97c92e
# ╠═a264a212-f6a3-41ac-a80b-7567ded2af0d
# ╠═368c4e24-f367-4383-a64d-c104ef346776
# ╠═4d8b56d4-1ac2-4e37-8f2b-e59173514fbd
# ╠═b05bee10-c034-4e92-8233-20d388aaf9ca
# ╠═ac2da991-20fa-4922-99c0-e28d4f395300
# ╠═e2de29a2-61ac-412e-beb3-bf019d437214
# ╠═8933349f-91c0-4329-b3e1-7aac01ffb557
# ╟─edb845bc-7e73-45fa-bd88-de9f38559407
# ╟─562ebf36-19a7-48af-8825-cbab5ce7e3b7
# ╟─29fe8937-216b-4c30-9808-dc95d6208b1f
# ╠═63a54c1b-3b55-4b79-aa61-8ffa75c0429c
# ╠═639ca6d6-e75f-41b8-870c-bc736c0f758b
# ╠═c4268ea2-fdfd-47ff-bf99-206f67cdaec7
# ╠═f71fac30-37a2-49d8-a12b-94fd71534d57
# ╠═db3c7664-a532-466f-b65a-15800d5c6363
# ╠═da36394d-a8b1-4914-99e4-826d10b1a168
# ╠═6b73fcb5-535d-498e-b500-abbc27b98173
# ╠═07932901-f0c4-4edc-9fc3-0c8d7af25282
# ╠═e1787274-9153-47e9-8ba0-7d68283a4dfc
# ╠═85737ce4-50d2-40fa-955c-ebf3488ca8bf
# ╠═4095cee2-f2f0-4a59-aa24-1f520f8c1360
# ╠═cf7d64e7-4af8-47e9-af48-d58595afec61
# ╠═d9847db0-c937-45df-b3c1-0e976898b091
# ╠═4c9befae-555d-4a51-a372-225085b90dc8
# ╠═beefb08e-e8f2-4e65-b6f6-1562ba63847d
# ╠═3fb7a6b2-e55a-47b4-8578-f2b32071ca2f
# ╠═3370d88d-b2e2-446f-b8ee-16e4b9c1dbec
# ╠═54dd1897-faae-4920-aa12-712996254385
# ╠═c91a0c19-1c98-424f-a37b-049083db3acb
# ╠═5d467edc-74de-4b8e-a13b-d104c052b8ac
# ╠═0ff2061a-7f07-4f58-ba34-4a9536bafed5
# ╠═0edd6a6d-41b4-4523-84c4-5667fd9bf126
# ╠═7106c47a-d6c2-47b5-a899-347268faf94c
# ╠═938eb752-fbed-4545-bb7b-4e3b6821b946
# ╠═ee5a9a7c-38f7-4b3c-9591-8a34cf772fec
# ╠═06930724-252f-42f9-8c73-7fcfcbed9b76
# ╠═c2434331-c6ee-45bf-8d9e-049c1a827bd2

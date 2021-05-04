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
LabPlots.plot_xy(Lx/nm, Ex/eV, 
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

# ╔═╡ 586958c7-c69a-483a-b7c8-eb065db96f3b
md"## Field of View (FOV)
- Corresponds to the region (area, volume) iluminated by the laser
- It is defined by a diameter and a length (or thickness). 
- For example, for a gaussian laser beam, the FOV can be defined in terms of the σ of the beam (in xy and in z)
"

# ╔═╡ 013420a8-9129-11eb-0d7f-ab8fb7c70786
md"#### FOV and photon density"

# ╔═╡ bfa4e6ef-469d-4552-b1b3-8bc1240dc1eb
md"- Define the diameter of the FoV"

# ╔═╡ 274494c9-7ec3-4635-8e7d-6b5bbf27cef0
@bind fovd NumberField(1:10; default=1)

# ╔═╡ c89247fa-e827-4cea-b080-50df89bcfc20
md"- Define the z of the FoV"

# ╔═╡ 8fe41c56-b0c1-4635-8977-713f5e6b496e
@bind fovz NumberField(1:10; default=1)

# ╔═╡ 162d5452-8741-49aa-a219-c9f5aec83b82
md"- Define the units"

# ╔═╡ 35357f1b-f1fa-4a88-b6be-176dfb695585
39182/6778

# ╔═╡ 0f84598d-32ad-43bc-8a9f-eb78a08d81b7


# ╔═╡ 910d2777-321a-4548-bafe-3314c8f29436


# ╔═╡ ef55c4be-1cfd-48f3-901e-301ae3c3b607
@bind fovu TextField(default="μm")

# ╔═╡ 7d25636f-a9e6-446a-879b-a053db2c310b
fstr(2, "%5.3g")

# ╔═╡ c1288245-be28-4216-9ca3-d6723dd6ff34
function unt(val, unit)
	ex = quote
	$val * $unit
	end
	eval(ex)
end

# ╔═╡ 513be477-3279-4169-b254-4058c1a2baa2
typeof(unt(2, mm))

# ╔═╡ 53f112f6-ba9f-4089-ae96-83baad162fcc
md"- Define the power of the laser (in mW)"

# ╔═╡ 990240b5-b30d-4c50-911c-42f01d7c98c3
@bind lp NumberField(10^-1:10^3; default=0.1)

# ╔═╡ 7ff6fe7c-9129-11eb-1805-0dd052ec01a8
begin
	fovnf = lfi.LabFbi.Fov(fovd*mm, 1μm)
	I = lfi.LabFbi.photon_density(405nm, lp*mW, fovnf.a)
end

# ╔═╡ d0c7c8d5-e94b-434b-8920-c81688dfa5de
begin
	fovas = lti.LabTools.to_fstr(fovnf.a/mm^2,"%5.2g")
	Is    = lti.LabTools.to_fstr(I/(Hz*mm^-2), "%5.2g")
	lps   = lti.LabTools.to_fstr(lp, "%5.2g")
	end

# ╔═╡ 877170f1-4f5c-438d-8949-854d3e892a2e
md" 
- FoV diameter (mm)       = $fovd 
- Fov area (``mm^{2}``)         = $fovas
- Laser power (mW)        = $lps
- Laser density (Hz/``mm^{2}``) = $Is
"

# ╔═╡ a80e6642-908c-11eb-3ac3-6737f11f4eff
md"#### Number of photons produced in solution:
- Concentration, Fov and laser power are tuned by parameters
- Use the measured cross section for FBI"

# ╔═╡ 9cb2c012-608a-4cff-b312-88bd2071da1f
md"#### Molar extinction coefficient
- λ in nm
- ϵ in M``^{-1}``cm``^{-1}``
"

# ╔═╡ cd707e24-df49-4f74-9890-0cdfcc1b9296
adf = LabFbi.load_df_from_csv(datadir("fbi"),
	                          "molar_extinction_coefficient_G1G2.csv", 
	                          LabFbi.enG) 


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
md"##### Emitted fluorescence"

# ╔═╡ 4ed47097-0d51-4639-a19f-38bbaa6c508d
begin
	γfbi325G1 = lfi.LabFbi.fluorescence(fbi325G1, I)
	γfbiba325G1 = lfi.LabFbi.fluorescence(fbiba325G1, I)
	γfbi405G1 = lfi.LabFbi.fluorescence(fbi405G1, I)
	γfbiba405G1 = lfi.LabFbi.fluorescence(fbiba405G1, I)
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

# ╔═╡ dc32d330-8cb5-498a-a9f5-a9cd2359baf4
md"## Solution"

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

# ╔═╡ 8b918c37-dcc2-4f23-b0e5-bced29d3c654
md"- The concentration of FBI in solution is $(lti.LabTools.to_fstr(cs, \"%7.1g\")) M"

# ╔═╡ 5347f077-66b4-45e6-a2f1-b0e6f52b0768
solfbi = lfi.LabFbi.Solution("FBI solution", cs*M)

# ╔═╡ 844b4b75-4285-4449-9635-64d00e2443b1
ρfbi = nofv(solfbi, vol::Unitful.Volume)

# ╔═╡ e9185f30-23d5-493c-9ac3-7b08201d8354
md"- Number of molecules per unit volume"

# ╔═╡ e39b6c3b-1b35-4e3d-a77f-ba1f8a574a86
md"## Silica"

# ╔═╡ dd9c1fd6-742f-4548-ad05-22977ccb2608
md"- The following two fields are used to define the concentration of the FBI solution in silica:
- First field: set the value in front of exponent of concentration in mmol/mg
- Second field: set the power"

# ╔═╡ 933a08d7-d9c9-4921-a9f8-cb6a90eb6940
@bind c1 NumberField(2:0.1:3; default=2.27)

# ╔═╡ ae44aa4d-54bd-42d7-be33-58f7980832ad
@bind c1p NumberField(1e-8:1e-5; default=1e-5)

# ╔═╡ 8f910b3e-bfa8-4774-a279-02db079e275c
cc = c1 * c1p ;

# ╔═╡ 1ee7d4cf-b25f-4b04-a194-7d729fa86beb
md"- The concentration of FBI in the silica powder is $(lti.LabTools.to_fstr(cc, \"%7.1g\")) mmol/mg"

# ╔═╡ 991eb30d-c5eb-4ed6-978e-1556961df89a
fbiSiPo = lfi.LabFbi.Powder("FBI-SiPo", cc*mmol/mg, 30mg/cm^2)

# ╔═╡ 0973ae7d-cb35-4e0d-95cb-fe0258114e27
nmFovFbiSiPo = lfi.LabFbi.nofa(fbiSiPo, fovnf.a)

# ╔═╡ 17e37f70-41da-4400-9c23-f2c85ee8927c
md"- The numberof FBI molecules in a spot of  the silica powder is $(lti.LabTools.to_fstr(nmFovFbiSiPo, \"%7.2g\"))"

# ╔═╡ c620fc65-fa0e-47fe-a8dc-b04e2d7064b7
begin
γfovfbi325G1    =  γfbi325G1 * nmFovFbiSiPo
γfovfbiba325G1  =  γfbiba325G1 * nmFovFbiSiPo
γfovfbi405G1    =  γfbi405G1 * nmFovFbiSiPo
γfovfbiba405G1  =  γfbiba405G1 * nmFovFbiSiPo
end

# ╔═╡ 6f6e61db-e9b1-43bc-a68a-38aba9fccfb9
begin
	γfovfbi325G1s = lti.LabTools.to_fstr(γfovfbi325G1/Hz,"%7.3g")
	γfovfbiba325G1s = lti.LabTools.to_fstr(γfovfbiba325G1/Hz,"%7.3g")
	γfovfbi405G1s = lti.LabTools.to_fstr(γfovfbi405G1/Hz,"%7.3g")
	γfovfbiba405G1s = lti.LabTools.to_fstr(γfovfbiba405G1/Hz,"%7.3g")
end

# ╔═╡ b844a604-a9b4-4b83-8d0c-c4711c3bc170
md""" #### Fluorescence emitted by Silica in FoV

| λ (nm) | Fbi G1 (Hz) | FbiBa G1 (Hz)
|:---------- | ---------- |:------------:|
| 325    | $γfovfbi325G1s| $γfovfbiba325G1s | 
| 405    | $γfovfbi405G1s| $γfovfbiba405G1s   | 
"""

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
ga = lfi.LabFbi.geometrical_acceptance(d_target_ccd/mm, D_ccd/mm)

# ╔═╡ 3d52a298-7ee2-4237-9f1f-39448d624867
md"- The geometrical acceptance for the unfocused setup (Laser light is collected by a lens of $(D_ccd/mm) mm located at $(d_target_ccd/mm) mm from the target) is ga = $(lti.LabTools.to_fstr(ga, \"%7.2g\")) "

# ╔═╡ a7bd7368-6bf3-4bac-9483-3b61528f20c7
begin
nfovfbi325G1    =  γfbi325G1 * nmFovFbiSiPo   * ga
nfovfbiba325G1  =  γfbiba325G1 * nmFovFbiSiPo * ga
nfovfbi405G1    =  γfbi405G1 * nmFovFbiSiPo   * ga
nfovfbiba405G1  =  γfbiba405G1 * nmFovFbiSiPo * ga
end

# ╔═╡ e98c14f8-433a-4f58-a2b2-376b27570b6f
begin
	nfovfbi325G1s = lti.LabTools.to_fstr(nfovfbi325G1/Hz,"%7.3g")
	nfovfbiba325G1s = lti.LabTools.to_fstr(nfovfbiba325G1/Hz,"%7.3g")
	nfovfbi405G1s = lti.LabTools.to_fstr(nfovfbi405G1/Hz,"%7.3g")
	nfovfbiba405G1s = lti.LabTools.to_fstr(nfovfbiba405G1/Hz,"%7.3g")
end

# ╔═╡ 1f0e38e5-3583-4d41-b95d-a7934ccb6a4c
md""" #### Fluorescence Reaching CCD, no objective only geometrical acceptance

| λ (nm) | Fbi G1 (Hz) | FbiBa G1 (Hz)
|:---------- | ---------- |:------------:|
| 325    | $nfovfbi325G1s| $nfovfbiba325G1s | 
| 405    | $nfovfbi405G1s| $nfovfbiba405G1s   | 
"""

# ╔═╡ f96e9b12-9139-11eb-08cb-5319bca7a01d
md"## Setup
- Setup includes an objective (if present), a CCD to collect light and a collection of filters"

# ╔═╡ b4badd12-8671-11eb-1f1b-69b36be8f7e4
begin
	naL = collect(0.1:0.05:1)
	obL = [lfi.LabFbi.Objective("Nikon-SLVD", na, 100) for na in naL]
	trL = [lfi.LabFbi.transmission(ob) for ob in obL]

	plot(naL, trL, leg=false, lw=2)

	xlabel!("numerical aperture")
	ylabel!("Transmission")
	title!("Tranmission as a function of NA")
end

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

	-The 405 nm BP filter (not needed for this calculation, since the beam is filtered before reaching the optical path, the filter contributes, together with other factors to a normalisation factor regarding power that needs to be computed)
	-The 425 nm LP filter, which enters twice in the optical path and eliminates high frequencies (below 425 nm).
	-The double noth which selects a region between 405 and 522 nm 
	-The 430 nm BP filter, which selects the signal (FbIBa2+ peak) region

- LP-450: Select light in the long pass Band 450 nm. This is done using the same combination than before except the last filter which is a 450 nm LP filter. This selects primarily the region of background (FBI).

- In addition one needs to add the efficiency of the CCD to have the complete transport function for both cases
"

# ╔═╡ 5eaa4bc1-25d4-4a0e-837d-4c1f24f868c1
fdn405 = fdn405_522

# ╔═╡ c126f2d9-0ccb-4b5c-9600-e14ad96e9591
filterSetBand430(λ) = map(x->flp425(x)^2 * fdn405(x) * fbp430(x) * orca(x), λ)

# ╔═╡ 7ccc9401-2fee-4d5e-a6ca-9a863c60d5a6
filterSetLP450(λ) = map(x->flp425(x)^2 * fdn405(x) * flp450(x) * orca(x), λ)

# ╔═╡ 9bd87f70-24e6-4476-acf7-7647face1b3e
flp425(430.0)^2 * fdn405(430.0) * fbp430(430.0) * orca(430.) ≈ filterSetBand430(430.0)

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

# ╔═╡ 23fabeb9-a947-4d44-9da2-a929cf9db6a4
md"
- Get the spectrum from FBI and FBIBa powder
- Notice that the discrimination between both species depends both of their relative brightness and their shape (e.g, chromatic separation between species). In Silica powder, with standard concentrations, the response of both species is described by the data below. 
- However, the response expected in a monolayer cannot be extrapolated from the silica measurements, since concentrations of free and chelated molecules can have an impact in their relative brightness. 
"

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

# ╔═╡ 95ce45d1-d000-46c4-9996-2aa6db2876cf
md"### Expected response"

# ╔═╡ cf7dc98f-da76-4017-a967-8c33455406bc
md"#### Filtered spectra"

# ╔═╡ e6149377-2aff-4c87-8ab0-226c1c65c378
lfi.LabFbi.fom(fbi::Gf, fbiba::Gf, λmin::Number, λmax::Number)

# ╔═╡ 30156129-845e-4261-b549-e6c576c749ef
begin
	fi430bp = quadgk(filteredFbi430BP, 375.0,  700.0)[1]
	fi450lp = quadgk(filteredFbi450LP, 375.0,  700.0)[1]
	rfbi = fi430bp / fi450lp
end

# ╔═╡ 2dc64611-1215-4ed0-a0e3-80536ff87b21
@test rfbi ≈ fbi.fbiratio(filteredFbi430BP, filteredFbi450LP, 375.0,  700.0,375.0,  700.0)

# ╔═╡ 3325a9f2-004f-4b5b-84bc-bb2618a51a43
begin
	rfbiba = fbi.fbiratio(filteredFbiBa430BP, filteredFbiBa450LP, 375.0, 700.0,
		                                                      375.0,  700.0)
	r2 = rfbiba / rfbi
end

# ╔═╡ fa2b678f-9721-4820-9823-c0bddc03ba30
md"- Filtered spectrum ratios and double ratio:  
- rfbi =$(utils.float_to_fstr(rfbi, \"%7.1g\")) 
- rfbiba =$(utils.float_to_fstr(rfbiba, \"%7.1g\"))
- r2 =$(utils.float_to_fstr(r2, \"%7.2g\"))"

# ╔═╡ 1a1b4a30-3966-40ab-8469-1509539defa2
md"#### Unfiltered specta"

# ╔═╡ e1c62a09-b62f-423f-a8b3-ba67bf308b92
begin
	rufbi = fbi.fbiratio(fbipdf, fbipdf, 425.0,  435.0, 450.0,  700.0)
	rufbiba = fbi.fbiratio(fbibapdf, fbibapdf, 425.0,  435.0, 450.0,  700.0)
	ur2 = rufbiba / rufbi
end

# ╔═╡ 8e5ea441-fe39-4e32-ae80-986331f786d9
md"- Unfiltered spectrum ratios and double ratio:  
- urfbi =$(utils.float_to_fstr(rufbi, \"%7.1g\")) 
- urfbiba =$(utils.float_to_fstr(rufbiba, \"%7.1g\"))
- ur2 =$(utils.float_to_fstr(ur2, \"%7.2g\"))"

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

# ╔═╡ 7dc8b84e-461a-46fb-8da0-e7ae12a74d78
nγga

# ╔═╡ 02e26731-8762-4680-8d15-3da475b3f3ce
nγFbiccdLP450 = nccd(nγga, filteredFbi450LP, 450.0, 700.0)

# ╔═╡ 1cfebc3a-4bbf-4f88-9084-c78e90861deb
nγFbiccdBP430 = nccd(nγga, filteredFbi430BP, 425.0, 435.0)

# ╔═╡ 7e189333-e1d4-4f84-9d6c-46c68359c2d7
nγFbiBaccdLP450 = nccd(nγga, filteredFbiBa450LP, 450.0, 700.0)

# ╔═╡ 6f4aebde-ffe8-4a58-b586-956a9323502e
nγFbiBaccdBP430 = nccd(nγga, filteredFbiBa430BP, 425.0, 435.0)

# ╔═╡ 98af40b7-836f-4fa8-8774-ca95ebafc4a6
r2BP430 = nγFbiBaccdBP430 / nγFbiccdBP430

# ╔═╡ de377d32-9c1a-4caa-a370-000b2fb34c83
md" #### Summary of results:

- Assume that the number of emitted photons is the same for both species (adjust otherwise).

- Rate of FBI at the CCD LP450   = $(utils.float_to_fstr(nγFbiccdLP450/Hz, \"%7.2g\"))
- Rate of FBI at the CCD BP430   = $(utils.float_to_fstr(nγFbiccdBP430/Hz, \"%7.2g\"))
- Rate of FBIBa at the CCD LP450 = $(utils.float_to_fstr(nγFbiBaccdLP450/Hz, \"%7.2g\"))
- Rate of FBIBa at the CCD BP430 = $(utils.float_to_fstr(nγFbiBaccdBP430/Hz, \"%7.2g\"))
- R2 CCD  = $(utils.float_to_fstr(r2BP430, \"%7.2g\"))
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
# ╠═586958c7-c69a-483a-b7c8-eb065db96f3b
# ╟─013420a8-9129-11eb-0d7f-ab8fb7c70786
# ╠═bfa4e6ef-469d-4552-b1b3-8bc1240dc1eb
# ╠═274494c9-7ec3-4635-8e7d-6b5bbf27cef0
# ╠═c89247fa-e827-4cea-b080-50df89bcfc20
# ╠═8fe41c56-b0c1-4635-8977-713f5e6b496e
# ╠═162d5452-8741-49aa-a219-c9f5aec83b82
# ╠═35357f1b-f1fa-4a88-b6be-176dfb695585
# ╠═0f84598d-32ad-43bc-8a9f-eb78a08d81b7
# ╠═910d2777-321a-4548-bafe-3314c8f29436
# ╠═ef55c4be-1cfd-48f3-901e-301ae3c3b607
# ╠═7d25636f-a9e6-446a-879b-a053db2c310b
# ╠═c1288245-be28-4216-9ca3-d6723dd6ff34
# ╠═513be477-3279-4169-b254-4058c1a2baa2
# ╟─53f112f6-ba9f-4089-ae96-83baad162fcc
# ╟─990240b5-b30d-4c50-911c-42f01d7c98c3
# ╠═7ff6fe7c-9129-11eb-1805-0dd052ec01a8
# ╠═d0c7c8d5-e94b-434b-8920-c81688dfa5de
# ╠═877170f1-4f5c-438d-8949-854d3e892a2e
# ╟─a80e6642-908c-11eb-3ac3-6737f11f4eff
# ╠═9cb2c012-608a-4cff-b312-88bd2071da1f
# ╠═cd707e24-df49-4f74-9890-0cdfcc1b9296
# ╠═482b76cf-b3fb-4a30-8678-56ebe539fc5e
# ╠═a35dd220-5ab1-4f70-901c-6a55f8e50123
# ╠═aaf97203-1412-493f-94cb-75591ff981f7
# ╠═4ed47097-0d51-4639-a19f-38bbaa6c508d
# ╟─45fd4819-12bd-449e-9f37-5885cf3604ca
# ╠═7f595bfa-71d6-42f3-913c-b7574c0d4fc2
# ╠═dc32d330-8cb5-498a-a9f5-a9cd2359baf4
# ╟─ccb649f8-0516-4a64-8ccc-5af13a7fe6e1
# ╟─a6db60d8-371d-4b93-b164-38d44f00d17a
# ╟─96639a5e-81c2-4eac-837c-fff5bc7ea5aa
# ╟─639922e2-c0e5-48ee-a2b2-9e3afd543e1e
# ╟─205870db-bf67-4ece-8de4-239aaa5edc47
# ╠═21c6d2f8-519e-41ff-8895-ae3c77eb7a02
# ╟─8b918c37-dcc2-4f23-b0e5-bced29d3c654
# ╠═5347f077-66b4-45e6-a2f1-b0e6f52b0768
# ╠═844b4b75-4285-4449-9635-64d00e2443b1
# ╠═e9185f30-23d5-493c-9ac3-7b08201d8354
# ╠═e39b6c3b-1b35-4e3d-a77f-ba1f8a574a86
# ╠═dd9c1fd6-742f-4548-ad05-22977ccb2608
# ╠═933a08d7-d9c9-4921-a9f8-cb6a90eb6940
# ╠═ae44aa4d-54bd-42d7-be33-58f7980832ad
# ╠═8f910b3e-bfa8-4774-a279-02db079e275c
# ╠═1ee7d4cf-b25f-4b04-a194-7d729fa86beb
# ╠═991eb30d-c5eb-4ed6-978e-1556961df89a
# ╠═0973ae7d-cb35-4e0d-95cb-fe0258114e27
# ╠═17e37f70-41da-4400-9c23-f2c85ee8927c
# ╠═c620fc65-fa0e-47fe-a8dc-b04e2d7064b7
# ╠═6f6e61db-e9b1-43bc-a68a-38aba9fccfb9
# ╠═b844a604-a9b4-4b83-8d0c-c4711c3bc170
# ╟─afe6ab7a-9133-11eb-1f9f-dd7c0b27e4d0
# ╠═67b2b186-9139-11eb-2f29-f9f36be38b6f
# ╠═55dc84a0-9139-11eb-1694-bba96f7bad9e
# ╟─3d52a298-7ee2-4237-9f1f-39448d624867
# ╠═a7bd7368-6bf3-4bac-9483-3b61528f20c7
# ╠═e98c14f8-433a-4f58-a2b2-376b27570b6f
# ╠═1f0e38e5-3583-4d41-b95d-a7934ccb6a4c
# ╟─f96e9b12-9139-11eb-08cb-5319bca7a01d
# ╟─b4badd12-8671-11eb-1f1b-69b36be8f7e4
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
# ╟─9bd87f70-24e6-4476-acf7-7647face1b3e
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
# ╠═23fabeb9-a947-4d44-9da2-a929cf9db6a4
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
# ╠═67709bf1-298a-4371-8f91-1b6fab98ef6b
# ╠═f613b185-0f04-4926-8e7d-f1d9d19046b5
# ╟─95ce45d1-d000-46c4-9996-2aa6db2876cf
# ╟─cf7dc98f-da76-4017-a967-8c33455406bc
# ╠═e6149377-2aff-4c87-8ab0-226c1c65c378
# ╠═30156129-845e-4261-b549-e6c576c749ef
# ╠═2dc64611-1215-4ed0-a0e3-80536ff87b21
# ╠═3325a9f2-004f-4b5b-84bc-bb2618a51a43
# ╟─fa2b678f-9721-4820-9823-c0bddc03ba30
# ╟─1a1b4a30-3966-40ab-8469-1509539defa2
# ╠═e1c62a09-b62f-423f-a8b3-ba67bf308b92
# ╟─8e5ea441-fe39-4e32-ae80-986331f786d9
# ╟─1e325c2e-893a-4136-9302-45a61aa469df
# ╟─5a48f61e-0de6-4b36-9f3e-6e539269ec03
# ╠═2197fe18-811a-4c2e-a1fc-daac14ff1fa0
# ╠═14d074ec-1880-449f-abfa-eab2f431fd76
# ╟─0bb2ec9b-1741-451b-8dd2-31fde86098ad
# ╠═7dc8b84e-461a-46fb-8da0-e7ae12a74d78
# ╠═02e26731-8762-4680-8d15-3da475b3f3ce
# ╠═1cfebc3a-4bbf-4f88-9084-c78e90861deb
# ╠═7e189333-e1d4-4f84-9d6c-46c68359c2d7
# ╠═6f4aebde-ffe8-4a58-b586-956a9323502e
# ╠═98af40b7-836f-4fa8-8774-ca95ebafc4a6
# ╟─de377d32-9c1a-4caa-a370-000b2fb34c83

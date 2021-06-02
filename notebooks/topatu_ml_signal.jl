### A Pluto.jl notebook ###
# v0.14.4

using Markdown
using InteractiveUtils

# ╔═╡ 26badf86-379f-4e66-b3c9-9dc122294a91
import Pkg;

# ╔═╡ c5119b0c-66e3-4fcc-b8f4-253b6e2d0c89
Pkg.add("LaTeXStrings")

# ╔═╡ f3ec8d6e-a8bf-11eb-1f98-636c338d90e7
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
	using Images
	using ImageIO
	using ImageView
end

# ╔═╡ 955d45a1-9d25-4a23-86c2-cc279d14e710
using DrWatson

# ╔═╡ c60716c4-108c-4c1e-85c1-024520065952


# ╔═╡ 0447a463-0ceb-49e9-a565-a47d9e881455
md"# Expected signal in the TOPATU setup using Monolayers"


# ╔═╡ 173989ad-8f8e-4f7e-a88a-9647bf63c844
@quickactivate "LabFBI"

# ╔═╡ 312e7f8b-9df5-4752-8de8-afff2484dad7
projectdir()

# ╔═╡ fe595465-577c-4352-ae1a-66ce5b5edbf8
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

# ╔═╡ 21712b78-ccf6-47fa-b450-13186d35cd5f
begin
	setup = ingredients(srcdir("setup.jl"))
	plt   = ingredients(srcdir("plotDF.jl"))
	dff   = ingredients(srcdir("dffunctions.jl"))
	utl   = ingredients(srcdir("utils.jl"))
	lfi = ingredients(srcdir("LabFbi.jl"))
	lti = ingredients(srcdir("LabTools.jl"))
	lpi = ingredients(srcdir("LabPlots.jl"))
end

# ╔═╡ 14f34300-f974-432d-8a53-91da2a0dc298
import Unitful:
    nm, μm, mm, cm, m, km, inch, ft, mi,
    ac,
    mg, g, kg,
    Ra, °F, °C, K,
    rad, °,
    ns, μs, ms, s, minute, hr, d, yr, Hz,
    eV,
    μJ, mJ, J,
	mW, W, μW,
    A, N, mol, mmol, V, L, M

# ╔═╡ 68c4adce-58c0-4807-b87b-10c4e8d4b9b4
markercolors = [:green :orange :black :purple :red  :yellow :brown :white]

# ╔═╡ 4c9c1952-9a1b-4bfa-b665-a05ad399e112
md"""## Calculation: 

Let us consider:

- A Multilayer of fluorophores (MLFn). By convention, we define $n=1$ (and speak of a monolayerm MLF1 or simply ML) when the fluorophore densitiy, $alpha$, is $\alpha = 1$ fluorophore/$nm^2$.  

- A laser of wavelength $\lambda$ and power P. 

- A setup, which includes an Objective, a set of filters and a CCD. We assume that the laser is focused, through the objective, in a spot limited by diffraction. 

Our goal is to compute the rate of photons recorded by the CCD due to the irradiance of the spot (which we will call FOV or field-of-view). For this we need to know:

- The FOV area, which yields the number of fluorophores in the spot, for a given $\alpha$. Since we assume that the FOV is limited by diffraction, the area is simply:
$\pi d^2$, where $d$ is the diameter of the diffractive spot:

$$d = 1.22 * \frac{\lambda}{2A}$$

and A is the numerical apeture (NA) of the obective. 

- The absorption cross section ($\sigma$) which gives the probability, per molecule, that a photon is absorbed by the fluorophore.

- The quantum efficiency Q of the fluorophore, which gives the fraction of the time that a fluorophore which has absorbed an excitation photon, emits a fluorescence photon. 

The number of fluorescence photons per molecule is:

$$N = \sigma Q I$$

Where $I$ is the photon density (photons/unit area) in the diffractive spot.

The total fluorescence emitted by the spot is:

$$F = N m$$

where $m = \alpha \pi d^2$ is the total number of molecules in the spot.

The fluorescence recorded by the CCD is:

$f = F \epsilon_o \epsilon_f \epsilon_c$

where $\epsilon_o$ is the efficiency due to the objective (computed to zero order as the geometrical transmitance of the objective), $\epsilon_f$ is the transmitance of the filter set and  $\epsilon_c$ is the efficiency of the CCD. Notice that, in the most general case, $f = f(\lambda)$

All the results shown below are computed in the notebook
"""

# ╔═╡ bad2ab3d-ac3d-412a-8cfa-0802baa31cb4
md"#### Number of photons
- Left Plot: Photon energy as a function of λ
- Right Plot: Number of photons as a function of power, for a laser of 405nm
"

# ╔═╡ bcd854b0-9918-478e-8db6-656b87652eed
md"#### Molar extinction coefficient
Table below gives the molar extinction coefficient (ϵFbiG1,...) with errors (σFbiG1,...) for two different wavelenghts
- λ in nm
- ϵ in M``^{-1}``cm``^{-1}``
"

# ╔═╡ bc2b045a-e211-4294-b79e-0cb37bc5cd41
adf = dff.load_df_from_csv(datadir("fbi"),
	                          "molar_extinction_coefficient_G1G2.csv", 
	                          dff.enG)

# ╔═╡ af3fbdc8-0d4e-47f9-9137-2d0d11590300
md"## Fluorescence recorded by the CCD"

# ╔═╡ d1c1bab7-20ba-48ec-ac7c-50fa937a7491
md"### Topatu setup"

# ╔═╡ 6d130ae9-4d2b-43cb-a0a7-fab6a411376c
load("laserFBI.png")

# ╔═╡ 3ad8f3f1-4d22-41c4-8719-425c443cdc4f
md"### Fluorescence Transmitted through the objective"

# ╔═╡ 56ada632-74bf-420f-ad7c-412d4cd29e9f
md"### Filters"

# ╔═╡ c0abf811-c32d-4d03-b44f-f0baa4783bcc
md"## Code for the calculation"

# ╔═╡ b81f7da7-0489-4468-a0d1-fa825472ff94
begin
    wf = 350.
    we = 800.
    ws = 2.
    wl = wf:ws:we
	LL = collect(wl)
	println("")
end

# ╔═╡ 71488556-1030-4cf4-8a4e-c3617c31b65e
md"#### Laser setup"

# ╔═╡ 06a98663-2546-45fb-b2d4-d6eae5335260
md"#### Fluorescence per molecule for FBI"

# ╔═╡ ac8a4e49-e1ec-4fa9-99c6-3d5829d4fbcf
md"#### Emitted fluorescence in the diffractive spot"

# ╔═╡ 7e2fac7a-f435-423b-98e2-6facbdcfd5d7
md"#### Objective tranmission and CCD transmission"

# ╔═╡ d8823f39-9023-4f78-9b80-8a0352dd5e2f
begin
	T = setup.transmission(0.6)
	orca = setup.ccd()
	tccd = orca(405.0)
	tx = utl.to_fstr(T,"%5.3g")
	tcx = utl.to_fstr(tccd,"%5.3g")
	println(" ")
end

# ╔═╡ d2219407-6114-4be3-b369-15d7bb12a1be
md"The NA of the objective is 0.6. As it can be seen in the graphic this corresponds to a transmission of $tx "

# ╔═╡ 99b59d73-2303-4f53-8b12-c40e6726a302
plt.plot_xy(LL, orca.(LL), 
	        "λ(nm)", "ϵ", "CCD efficiency as a function of λ")

# ╔═╡ f1592f8b-7436-42d8-b110-6fdba86e808d
md"#### Filters"

# ╔═╡ 13e869f8-132d-4ac9-9e8d-73744884da47
begin
	filterd = lfi.LabFbi.get_filters(lfi.LabFbi.fnames, datadir("filters"))
	fw = [390:420,420:500,420:520,440:540, 390:560, 390:450]
	Pft = lpi.LabPlots.plot_filters(lfi.LabFbi.fnames, filterd, fw);
	f425 = filterd["LongPass425"]
	f425(500.)
	tf425 = utl.to_fstr(f425(500.),"%5.3g")
	println("")
end

# ╔═╡ 014ace2d-3634-4377-9d97-eb109be95338
plot(Pft..., layout = (3, 2), legend=false, fmt = :png)

# ╔═╡ 09209559-64a4-4227-8997-38d46e64391f
md"#### Common filters multiplied by CCD efficiency"

# ╔═╡ 7db1d260-304f-4326-8ef4-9ae9e9311aaa
begin
	cflt =lfi.LabFbi.common_filters(filterd)
	filterSetAll(λ) = map(x->cflt(x) * orca(x), λ)
end

# ╔═╡ b4dd686d-9246-48c5-9504-9bd3b9e9f223
pfilterSetAll =lpi.LabPlots.plot_filterset(filterSetAll, 
	           "filterSetALL", 400.0:700.0);


# ╔═╡ e71cea3f-7dd4-4d20-94b4-aa37d781a14d
md"#### FBI spectrum"

# ╔═╡ 5a7bd84d-d1cc-4d1a-850a-de42e6d7c429
fbi405df = lfi.LabFbi.load_df_from_csv(datadir("fluorimeter/fluorescence"), "FbiG1Em405.csv", lfi.LabFbi.spG);


# ╔═╡ e49bf918-1a33-45bd-a020-f8c42e0146f4
md"#### Final calculation"

# ╔═╡ f4d2c123-4213-467e-ad8e-fd0ca90d2ee1
md"## Functions"

# ╔═╡ 93205bc6-1a2c-42d5-aba6-4c0845617f35
struct LaserSetup
	laser   ::setup.Laser
	obj     ::setup.Objective
	nphot   ::typeof(1.0Hz)
	phote   ::typeof(1.0eV)
	dlimit  ::typeof(1.0nm)
	fov     ::setup.Fov
	Iphot   ::typeof(1.0Hz*μm^-2)

	function LaserSetup(laser, obj)
		nphot   = setup.n_photons(laser)
		phote   = setup.photon_energy(laser.λ)
	    dlimit  = setup.diffraction_limit(laser, obj)
		fov     = setup.Fov(dlimit, dlimit)
		Iphot   = setup.photon_density(laser, fov)
		new(laser, obj, nphot, phote, dlimit, fov, Iphot)
	end
	
end

# ╔═╡ 00d79f5c-8e2f-4760-b6ea-562975bc363f
function laser_info_md(lsp)
	md = " ### Laser Setup \n
- Objective NA = $(lsp.obj.NA)\n"
	w =  utl.to_fstr(lsp.laser.λ/nm,"%5.3g")
	lp = utl.to_fstr(lsp.laser.P/mW,"%5.3g")
	pe = utl.to_fstr(lsp.phote/eV,"%5.3g")
	np = utl.to_fstr(lsp.nphot/Hz,"%5.3g")
	dl = utl.to_fstr(lsp.dlimit/nm,"%5.3g")
	fa = utl.to_fstr(lsp.fov.a/nm^2,"%5.3g")
	I =  utl.to_fstr(lsp.Iphot/(Hz* cm^-2),"%5.3g")
	
	md = string(md, "- laser λ        = $w nm \n")
	md = string(md, "- laser power    = $lp mW \n")
	md = string(md, "- photon energy  = $pe eV \n")
	md = string(md, "- photon rate    = $np Hz \n")
	md = string(md, "- diff limit (d) = $dl nm \n")
	md = string(md, "- Fov  area      = $fa nm2 \n")
	md = string(md, "- photon density = $I  Hz/cm2 \n")
	
	return md
end


# ╔═╡ bc4458c9-99cd-45f1-8c57-15bb21863319
begin
	P    = 0.1mW
	l405 = setup.Laser(405.0nm, P)
	tobj = setup.Objective("Topatu objective", 0.6, 100.0)

	lst   = LaserSetup(l405, tobj)
	linfo = laser_info_md(lst)
	println("")
end

# ╔═╡ 162a84d5-0830-4702-a974-095f84382eb5
Markdown.parse("$linfo")

# ╔═╡ 6dff150f-7a6a-4273-93ee-8c98751b3ed3
begin
	fbiFluo = lfi.LabFbi.fbi_fluorophores(adf, ["g1", "g2"], 
		                                  [325,405], [0.67,0.67])
	fbi405nm = fbiFluo.fbi["l405g1"]
	fbiba405nm = fbiFluo.fbiba["l405g1"]
	Nfbi = lfi.LabFbi.fluorescence(fbi405nm, lst.Iphot)
	Nfbiba = lfi.LabFbi.fluorescence(fbiba405nm, lst.Iphot)
	println(" ")
end

# ╔═╡ a8beb6de-7aab-443b-bdfe-f91cc4db9b92
begin
	nfluo = lst.fov.a/nm^2
	Ffbi = nfluo * Nfbi
	Ffbiba = nfluo * Nfbiba
	nfluox = utl.to_fstr(nfluo,"%5.3g")
	Ffbix = utl.to_fstr(Ffbi/Hz,"%5.3g")
	Ffbibax = utl.to_fstr(Ffbiba/Hz,"%5.3g")
	println(" ")
end

# ╔═╡ 88e59d90-8a3a-4023-bd2a-e33630ab097c
md"## Fluorescence emitted by the spot

If we assume a monolayer (1 fluorophore per nm2) then the number of fluorophores in the diffractive spot equals the fov area. Thus:

nfluo = $nfluox

The total fluorescence emitted by the spot is the product of the fluorescence per molecule and the number of molecules.

- For FBI: F = $Ffbix Hz
- For FBIBa: F = $Ffbibax Hz
"

# ╔═╡ 36ff3df2-fdc2-4f76-9f3a-72caeb5bb0f6
function fluo_md(fbiFluo, select="l405g1")
	fbi405nm = fbiFluo.fbi[select]
	fbiba405nm = fbiFluo.fbiba[select]
	Nfbi = lfi.LabFbi.fluorescence(fbi405nm, lst.Iphot)
	Nfbiba = lfi.LabFbi.fluorescence(fbiba405nm, lst.Iphot)
	
	md = " ### Photons per fluorophore \n"
	fbie =  utl.to_fstr(fbi405nm.ϵ/(cm^-1*M^-1),"%5.3g")
	fbiq =  utl.to_fstr(fbi405nm.Q,"%5.3g")
	fbis =  utl.to_fstr(fbi405nm.σ/(cm^2),"%5.3g")
	fbae =  utl.to_fstr(fbiba405nm.ϵ/(cm^-1*M^-1),"%5.3g")
	fbaq =  utl.to_fstr(fbiba405nm.Q,"%5.3g")
	fbas =  utl.to_fstr(fbiba405nm.σ/(cm^2),"%5.3g")
	nfbi =  utl.to_fstr(Nfbi/(Hz),"%5.3g")
	nfbiba =  utl.to_fstr(Nfbiba/(Hz),"%5.3g")
	
	
	md = string(md, "- FBI ϵ   = $fbie ``cm^{-1} M^{-1}`` \n")
	md = string(md, "- FBI Q   = $fbiq  \n")
	md = string(md, "- FBI σ   = $fbis ``cm^2``  \n")

	md = string(md, "- FBIBa ϵ = $fbae ``cm^{-1} M^{-1}`` \n")
	md = string(md, "- FBIBa Q = $fbaq  \n")
	md = string(md, "- FBIBa σ = $fbas ``cm^2``  \n\n")
	
	md = string(md, "- FBI: photon rate per fluorophore = $nfbi Hz  \n")
	md = string(md, "- FBIBa: photon rate per fluorophore = $nfbiba Hz  \n")

	return md
end

# ╔═╡ 149a48a0-d5bd-432b-936d-7aa87b7b847f
begin
	finfo = fluo_md(fbiFluo)
	Markdown.parse("$finfo")
end

# ╔═╡ 59a4d332-39a4-42f0-bdeb-cb907939fde2
function plot_photon_energy_vs_wavelength(wl=250:800)
	Lx = collect(wl) * nm
	Ex = setup.photon_energy.(Lx)
	pewl = plt.plot_xy(Lx/nm, Ex/eV, 
	                            "λ (nm)", "E (eV)", "Photon energy")
	return pewl
end

# ╔═╡ 710ad954-750d-43f2-bdbf-ffc322d2f1d8
function plot_number_photons_vs_power(laser, title)
	Px = collect(utl.logrange(10^-2, 10^3, 100)) *mW
	Np = setup.n_photons.((laser.λ,), Px)
	pnp = plt.loglog_xy(Px/mW, Np/Hz, 
	        "P (mW)", "N", title)
	return pnp
end

# ╔═╡ 05f575a6-b242-4ac6-b3de-a7653d1f3e77
begin
	pewl = plot_photon_energy_vs_wavelength()
	pnp405 = plot_number_photons_vs_power(l405, "Number of γ, laser 405 nm")
	println(" ")
end

# ╔═╡ a7924698-401d-458b-a93f-8db833cbef6e
plot(pewl,pnp405, layout = (1, 2), legend=false, fmt = :png)

# ╔═╡ eca4fcfd-0997-4bef-bdea-bf4b4844bdd2
begin
	Nx = collect(0.1:0.005:1.0)
	Tx = setup.transmission.(Nx)
end

# ╔═╡ 6d39c50b-e465-4ac3-8903-10fc1207c8a8
plt.plot_xy(Nx, Tx, 
	        "NA", "T", "Transmission as a function of Numerical Aperture")

# ╔═╡ 396020ee-9857-4cd2-ab15-dce3f2c688a9
function dfhead(df, fbi=true)
	nx = split.(names(df),"_")
	nz = [n for n in nx[2:end]] 
	if fbi
		return nz[1][1], nz[1][3], [n[2] for n in nz]
	else
		return nz[1][1]*nz[1][2], nz[1][4], [n[3] for n in nz]
	end
end

# ╔═╡ 416cd5a0-444f-40a7-b919-31e1b020b3c4
function dflm(df)
	l = df[!, "λ"]
	return l[1], l[end]
end

# ╔═╡ 83a658b5-452b-422a-93fd-34953a73db5e
function cnames(fbidf, nx)
	return names(fbidf)[nx]
end

# ╔═╡ b94df023-2852-4693-bcea-8100d85d9492
struct FbiDfInfo
	name::String 
	nλ::String
	C::Array{String}
	CN::String
	λi::Float64
	λf::Float64
end

# ╔═╡ 64c77fa9-e9bb-4240-9ba0-e03424a826be
function fbidfinfo(df, fbi=true,nx=2:6)
	fbili, fbilf, = dflm(df)
	fbiN, fbiL, fbiC = dfhead(df, fbi)
	fbiCs = lti.LabTools.vect_to_fstr(fbiC, "%s")
	Cnames = cnames(df, nx)
	return FbiDfInfo(fbiN,fbiL,Cnames,fbiCs,fbili,fbilf)
end

# ╔═╡ b7cbbb87-c0a0-4064-9029-b47f2756b6bd
function fbidfinfomd(dfi)
	md ="\n
- df name = $(dfi.name)
- df λ    = $(dfi.nλ)
- Concentration (in M) = $(dfi.CN)
- λmin = $(dfi.λi) nm, λmax = $(dfi.λf) nm
	"
	return md
end

# ╔═╡ 0a474bc8-21fb-4434-8c77-6650872c6b6f
function ncs(CN, ws=1:5)
	sC = lstrip.(split(CN,","))
	Cs =parse.((Float64,), sC)
	return (Cs / Cs[1])[ws]
end

# ╔═╡ 01c93d1a-5fa9-44f9-bccb-32a94c5c66f6
begin
	fbig1405dfi = fbidfinfo(fbi405df)
	fbig1405dfimd = fbidfinfomd(fbig1405dfi)
	nCs = ncs(fbig1405dfi.CN, 1:5)
	Cs = fbig1405dfi.C
	gfbi405 = lfi.LabFbi.dftogf(fbig1405dfi.λi:2.0:fbig1405dfi.λf, fbi405df,        
		                        fbig1405dfi.C[1])
	println("")
end

# ╔═╡ 2afe901e-bce0-4da4-ac61-43022b30fe9a
pfbi405 =lpi.LabPlots.plotdf_gfs([gfbi405], 350.:2.0:800.0,
	                             [fbig1405dfi.C[1]], markercolors, 
	                    "λ (nm)", "I (au)", "FBIG1",
		                legend=:topright);

# ╔═╡ c034c35a-89b5-44fd-9a0e-b817d1f8fc89
begin
	FBI(λ)  = map(x->cflt(x) * orca(x) * gfbi405.pdf(x),  λ)
	efft = lfi.LabFbi.qpdf(FBI, 400.0, 700.0)
	efftx = utl.to_fstr(efft,"%5.3g")
	println("")
end

# ╔═╡ dbca0013-8668-472f-8d3d-df75c220e7ec
md"### FBI spectrum

- The final ingredient for the calculation is the spectrum of FBI (and FBIBa).
- Left plot shows the efficiency of the filter as a function of λ (including the common filters and the CCD).

- Center  plot shows the normalized FBI spectrum (pdf FBI) 
- Right plot shows the normalized FBI spectrum convoluted with the filter

To compute the efficiency of transmission (the effect of the filters convoluted with the efficiency of the CCD), we integrate the topright plot. The result is:

eff_filters = $efftx

"

# ╔═╡ cea72215-3c24-46fb-9649-5254138ff08f
pfbi = lpi.LabPlots.plot_xy(LL, FBI.(LL), "λ(nm)", "I (norm)", "Filtered FBI");

# ╔═╡ 0fa51c71-bb2a-4bfa-bdae-c5e470ca2745
plot(pfilterSetAll, pfbi405, pfbi, layout = (1, 3), legend=false, fmt = :png)

# ╔═╡ 501c7cef-f1a5-4da2-b013-69ba05d9a1a8
begin
	f = Ffbi * T * efft
	fx = utl.to_fstr(f/Hz ,"%5.3g")
	println("")
end

# ╔═╡ 3b773ae3-c100-407b-8fc6-41128d28b451
md"### Fluorescence recorded in the CCD

$$f = F T_o T_f$$

Where F is the fluorescence, $T_o$ the objective tranmission and $T_f$ the transmission of filter + CCD. For FBI, we have computed the following values:

- F = $Ffbix Hz
- To = $tx
- Tf = $efftx

Thus the calculation yields that the efficiency recorded in the CCD is

- f = $fx Hz
"

# ╔═╡ Cell order:
# ╠═26badf86-379f-4e66-b3c9-9dc122294a91
# ╠═c5119b0c-66e3-4fcc-b8f4-253b6e2d0c89
# ╠═c60716c4-108c-4c1e-85c1-024520065952
# ╠═f3ec8d6e-a8bf-11eb-1f98-636c338d90e7
# ╠═955d45a1-9d25-4a23-86c2-cc279d14e710
# ╟─0447a463-0ceb-49e9-a565-a47d9e881455
# ╟─173989ad-8f8e-4f7e-a88a-9647bf63c844
# ╠═312e7f8b-9df5-4752-8de8-afff2484dad7
# ╟─fe595465-577c-4352-ae1a-66ce5b5edbf8
# ╟─21712b78-ccf6-47fa-b450-13186d35cd5f
# ╠═14f34300-f974-432d-8a53-91da2a0dc298
# ╠═68c4adce-58c0-4807-b87b-10c4e8d4b9b4
# ╟─4c9c1952-9a1b-4bfa-b665-a05ad399e112
# ╟─bad2ab3d-ac3d-412a-8cfa-0802baa31cb4
# ╠═a7924698-401d-458b-a93f-8db833cbef6e
# ╠═162a84d5-0830-4702-a974-095f84382eb5
# ╟─bcd854b0-9918-478e-8db6-656b87652eed
# ╠═bc2b045a-e211-4294-b79e-0cb37bc5cd41
# ╠═149a48a0-d5bd-432b-936d-7aa87b7b847f
# ╟─88e59d90-8a3a-4023-bd2a-e33630ab097c
# ╟─af3fbdc8-0d4e-47f9-9137-2d0d11590300
# ╟─d1c1bab7-20ba-48ec-ac7c-50fa937a7491
# ╠═6d130ae9-4d2b-43cb-a0a7-fab6a411376c
# ╟─3ad8f3f1-4d22-41c4-8719-425c443cdc4f
# ╠═6d39c50b-e465-4ac3-8903-10fc1207c8a8
# ╟─d2219407-6114-4be3-b369-15d7bb12a1be
# ╠═99b59d73-2303-4f53-8b12-c40e6726a302
# ╠═56ada632-74bf-420f-ad7c-412d4cd29e9f
# ╠═014ace2d-3634-4377-9d97-eb109be95338
# ╟─dbca0013-8668-472f-8d3d-df75c220e7ec
# ╠═0fa51c71-bb2a-4bfa-bdae-c5e470ca2745
# ╟─3b773ae3-c100-407b-8fc6-41128d28b451
# ╠═c0abf811-c32d-4d03-b44f-f0baa4783bcc
# ╠═b81f7da7-0489-4468-a0d1-fa825472ff94
# ╟─71488556-1030-4cf4-8a4e-c3617c31b65e
# ╠═bc4458c9-99cd-45f1-8c57-15bb21863319
# ╟─06a98663-2546-45fb-b2d4-d6eae5335260
# ╠═6dff150f-7a6a-4273-93ee-8c98751b3ed3
# ╟─ac8a4e49-e1ec-4fa9-99c6-3d5829d4fbcf
# ╠═a8beb6de-7aab-443b-bdfe-f91cc4db9b92
# ╟─7e2fac7a-f435-423b-98e2-6facbdcfd5d7
# ╠═d8823f39-9023-4f78-9b80-8a0352dd5e2f
# ╠═f1592f8b-7436-42d8-b110-6fdba86e808d
# ╠═13e869f8-132d-4ac9-9e8d-73744884da47
# ╟─09209559-64a4-4227-8997-38d46e64391f
# ╠═7db1d260-304f-4326-8ef4-9ae9e9311aaa
# ╠═b4dd686d-9246-48c5-9504-9bd3b9e9f223
# ╟─e71cea3f-7dd4-4d20-94b4-aa37d781a14d
# ╠═5a7bd84d-d1cc-4d1a-850a-de42e6d7c429
# ╠═01c93d1a-5fa9-44f9-bccb-32a94c5c66f6
# ╠═2afe901e-bce0-4da4-ac61-43022b30fe9a
# ╠═cea72215-3c24-46fb-9649-5254138ff08f
# ╟─e49bf918-1a33-45bd-a020-f8c42e0146f4
# ╠═c034c35a-89b5-44fd-9a0e-b817d1f8fc89
# ╠═501c7cef-f1a5-4da2-b013-69ba05d9a1a8
# ╠═f4d2c123-4213-467e-ad8e-fd0ca90d2ee1
# ╠═36ff3df2-fdc2-4f76-9f3a-72caeb5bb0f6
# ╠═93205bc6-1a2c-42d5-aba6-4c0845617f35
# ╠═00d79f5c-8e2f-4760-b6ea-562975bc363f
# ╠═05f575a6-b242-4ac6-b3de-a7653d1f3e77
# ╠═59a4d332-39a4-42f0-bdeb-cb907939fde2
# ╠═710ad954-750d-43f2-bdbf-ffc322d2f1d8
# ╠═eca4fcfd-0997-4bef-bdea-bf4b4844bdd2
# ╠═396020ee-9857-4cd2-ab15-dce3f2c688a9
# ╠═416cd5a0-444f-40a7-b919-31e1b020b3c4
# ╠═83a658b5-452b-422a-93fd-34953a73db5e
# ╠═b94df023-2852-4693-bcea-8100d85d9492
# ╠═64c77fa9-e9bb-4240-9ba0-e03424a826be
# ╠═b7cbbb87-c0a0-4064-9029-b47f2756b6bd
# ╠═0a474bc8-21fb-4434-8c77-6650872c6b6f

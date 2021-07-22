### A Pluto.jl notebook ###
# v0.14.7

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

# ╔═╡ bd59a5b8-8dac-4982-92e3-6d410b69e31c
begin
	import Pkg
	Pkg.add("StatsPlots")
end

# ╔═╡ eb573495-dfeb-471c-a0b2-11b2c76f6b25
Pkg.add("SpecialFunctions")

# ╔═╡ 25b471c6-e7b4-41a9-b7e6-bd316d78b0c7
Pkg.add("HCubature")

# ╔═╡ cc9d45d5-81d4-4525-add1-68f5adbbcc42
Pkg.add("GR")

# ╔═╡ f5a9c61a-c216-11eb-1ddc-6f24685e81bc
begin
	using Markdown
	using InteractiveUtils
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

# ╔═╡ 7c296957-7de3-411e-be97-bb3272c29102
begin
	using PhysicalConstants.CODATA2018
	using StatsPlots
	using Formatting
	using Printf
	using StrLiterals
	using StrFormat
	using Format
	using SpecialFunctions
	using HCubature
end

# ╔═╡ 65c05bfa-a565-41c1-aff7-71f96c7b9cbd
using DrWatson

# ╔═╡ 50df4350-0159-4a32-a4a3-2b196734853c
md"# Calibration of TOPATU with solutions: unfocused setup"

# ╔═╡ 10715061-650c-48a9-87af-d726b3444768
Pkg.update("GR")

# ╔═╡ 189dd02c-5243-49c3-88e4-fdb8f79380c9
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

# ╔═╡ ac29bcd2-37be-40f8-bbcf-095e3af08a7b
function void()
	println("")
end

# ╔═╡ a0c6b3d3-a43d-45c3-b511-9a21ce8e2390
@quickactivate "LabFBI"

# ╔═╡ 8c9a2ea3-6a7e-45ae-99a1-6db448064e35
projectdir()

# ╔═╡ 223d92a7-6db0-4383-848e-72f592c8627e
begin
	lfi = ingredients(srcdir("LabFbi.jl"))
	lti = ingredients(srcdir("LabTools.jl"))
	lpi = ingredients(srcdir("LabPlots.jl"))
end

# ╔═╡ 160eb718-f450-494a-a35b-1d3771d5b35e
markercolors = [:green :orange :black :purple :red  :yellow :brown :white]

# ╔═╡ 73f50e92-1c10-4d52-8fe3-f6fcd8e5b42a
gr(size=(500,500), xtickfontsize=12, ytickfontsize=12, xguidefontsize=14, yguidefontsize=14, legendfontsize=10, dpi=100, grid=(:y, :gray, :solid, 2, 0.4));


# ╔═╡ eb50ea90-be5e-4301-98b6-3871add430ff


# ╔═╡ 27c53ff7-e67e-4ce8-bcc1-adaeef8c2129
md"## Units"

# ╔═╡ 5869c147-6e13-4269-942a-fd84340011ba
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

# ╔═╡ 8077956e-84c6-4565-8647-c8bb197ae65e
md"## Gaussian Laser"

# ╔═╡ c449a041-0106-4cb3-a3b0-4b5689c5391d
load(projectdir("notebooks/gaussianBeam1.png"))

# ╔═╡ 189b1cae-a002-4b1e-9ca2-fd8e4ae00826
md"### Gaussian beam 

$g(r,z) = e^{-2(r/w(z))^2}$

$w(z) = w_0 \sqrt{1 + (\frac{z}{z_r})^2}$

$I(r,z) = I_0 [\frac{w_0}{w(z)}]^2 g(r,z)$

$z_r = \pi  \frac{w_0^2}{\lambda}$
The formulas verify the following limits:

$I(r, 0) = I_0 e^{-2(r/w_0)^2}$

$I(0, z) = I_0 [\frac{w_0}{w(z)}]^2 = \frac{I_0}{1 + (z/z_r)^2}$

Which has maximum value ($I_0$) at $z=0$ and then decreases as z increases. For $z=z_r$ the intensity along the axis is reduced by half. 

for $z << z_r, w(z) = w_0$, and:

$I(0, z) = I_0$

$I(r, z) = I_0 e^{-2(r/w_0)^2}$


Thus, the beam intensity falls exponentially with the radius r and is a constant ($I_0$) along the longitudinal axis, for $r=0$.

The beam power is obtained integrating the irradiance at any transverse plane. For simplicity we can choose $z=0$:

$P = \int_0^\infty I(r,0) 2 \pi r dr = 2 \pi I_0 \int e^{-2(r/w_0)^2} r dr$ 

Define $y = \frac{2 r^2}{w_0^2}$. Then $dy = \frac{4 r}{w_0^2} dr$, $dr = \frac{w_0^2}{4 r} dy$ and:

$P = 2 \pi I_0 \frac{w_0^2}{4} \int_0^\infty e^{-y} dy$ 

Since $\int_0^\infty e^{-y} dy = 1$, we find:

$P = I_0 \frac{\pi w_0^2}{2}$

or:

$I_0 = \frac{2 P}{\pi w_0^2}$ 

"

# ╔═╡ 780fb489-af23-4d82-bc6b-e00ed51d0b26
md"## Parameters"

# ╔═╡ 05584d03-2987-45bf-b24d-a188a56fe8d5
md" ##### Set laser wavelength (in nm)"

# ╔═╡ e53540b8-93b6-489f-b533-c34bf680e4a1
@bind ll NumberField(10.0^2:10.0^3; default=405.0)

# ╔═╡ e43cefef-23da-408a-b205-a7ee4296fb58
md" #####  Set laser power (in mW)"

# ╔═╡ de34e9c1-48e9-4d30-9c69-f3e5358881ae
@bind lp NumberField(10.0^-1:10.0^3; default=0.1)

# ╔═╡ 324cab76-4809-4ccb-b44b-361ae009eb4e
md" ##### Set $w_0$ at the objective lens or at the object plane (in mm)"

# ╔═╡ cb4e1c63-65c6-41f6-af91-319728046f85
@bind lw0 NumberField(10.0^-1:10.0^1; default=0.5)

# ╔═╡ 450184dd-8c86-47ca-b7b2-0ccdfa6e84ef
md" ##### Set thickness of the vase (in mm)"

# ╔═╡ 5b5ef870-3a41-48b1-99b8-6a95e90ad125
@bind lz NumberField(1.0:100.0; default=12.43)

# ╔═╡ 7c343c65-ac30-4662-b133-9055cb1ed8d5
md" ##### Set location of the CCD (in mm)"

# ╔═╡ f5f1f76b-9dc6-48c2-8e67-87a565cdf41b
@bind dccd NumberField(10.0:10.0^3; default=526.0)

# ╔═╡ 82c72876-4a85-403c-9642-96d2505aba32
md" ##### Set diameter of the CCD (in mm)"

# ╔═╡ ea8a0d64-4329-4080-80e0-0fa12ec1554b
@bind Dccd NumberField(10.0:20.0; default=13.3)

# ╔═╡ a8be2b0d-26f0-408d-a069-f92cc09d5514
md" ##### Set conversion ADC to PES"

# ╔═╡ d89bc3d6-8642-496b-a3f4-17f5c801512c
@bind adcpes NumberField(0.0:1.0; default=0.47)

# ╔═╡ dac46998-e6a2-41fb-adf7-c1d20d7493ff
md""" ##### Set concentration of test solution"""

# ╔═╡ 9602abd4-50cf-47dd-8bcd-450fc3ccf587
md" ##### Set Q (oxygenated) for fluorophore"

# ╔═╡ f9c273a2-b808-4e0e-b949-fe29d05c16d4
@bind Q NumberField(0.0:1.0; default=0.1)

# ╔═╡ c3a6a828-b18d-461c-bb17-2b13feaa593e
md"- Set units (in M)"

# ╔═╡ 0ac0eeb8-1bde-4a47-874d-f68649c19320
@bind csu NumberField(1.0:10.0; default=1.0)

# ╔═╡ 6d57e4ea-857a-414e-a848-dd5ef4541a91
md"- Set concentration power (in M)"

# ╔═╡ f3c687c7-03a2-4502-bbab-d6366121d747
@bind csp NumberField(1e-5:1e-9; default=1e-5)

# ╔═╡ b0e199e1-c1b9-48ce-bd25-9b48cd5ef9d1
cs = csu * csp 

# ╔═╡ 89a3cd82-cb35-4f68-89ac-06c50dbf71a9
md"- choose set of solutions (FBI or IrRu)"

# ╔═╡ e2b13a80-17bf-4647-b083-f3f9252fa1e4
@bind sset TextField(default="FBI")

# ╔═╡ 904130b2-f0e8-4f82-8afc-06977fc04801
md"#### Define solution"

# ╔═╡ b91be939-f092-4bf3-914c-fbe1b02dc9fd
solfbi = lfi.LabFbi.Solution("C = $cs*M", cs*M)

# ╔═╡ a969225a-1254-48ef-9775-8189ccfb9110
Cmol_mm3 = lfi.LabFbi.nofv(solfbi, 1*mm^3) / mm^3

# ╔═╡ 371dfb8c-4b37-427a-8be5-9448262ef65d
begin
	λ = ll * nm
	P = lp * mW
	w0 = lw0 * mm
	lzr = lz * mm
	void()
end

# ╔═╡ 0092e011-1abc-4cdb-a529-d926ae51c81e
gl = lfi.LabFbi.GLaser(λ, P, w0);

# ╔═╡ ba8e38da-c4ff-4638-a362-10ddad2e296e
md"""### Calibration setup

Consider now the problem of calibrating the TOPATU setup using a fluorescent solution of known fluorophores cross section (σ) and quantum yield (Q), with a given concentration (C). The solution is contained in a vial of quartz. The laser beam is collimated to a waist of $(gl.w0) mm and crosses the vial, which has a thickness $(lz) mm. A CCD camera of diameter $Dccd mm, located at a position $dccd mm, with respect to the vial, along the z axis, records the photons emitted by the solution. 

For a waist of $(gl.w0), zr = $(lti.LabTools.to_fstr(gl.zr/mm,"%5.2g")) and thus any point in the vial has a small z compared with zr. Therefore:

$w(z) = w_0$

$I(r, z) = I_0 e^{-2(r/w_0)^2}$

When the beam crosses the vial it iluminates all the fluorophores. The emitted rate of photons is found by integrating the irradiance and multiplying by the cross section, the quantum yield and the concentration:

$N = \sigma\, Q\, C \int_z \int_r I(r)\, 2 \pi r \, dr \,  dz$

The integral in $z$ is trivial and yields the length of the vial, $L$. The integral in $r$ can be precisely approximated by the integral to all the space, since the dimensions in r of the vial are much larger than the beam spread. Thus:

$N = I_0\, \sigma\, Q\, C\, L\, 2 \pi \int e^{-2(r/w_0)^2} r dr$

$N =  \frac{ I_0\, \sigma\, Q\, C\, L\, \pi w_0^2}{2} = P \sigma\, Q\, C\, L$

For a gaussian beam > 99.9 % of the power is contained in a radius $\rho = 2 w(z)$. In our case $w(z) \sim w_0$, and thus all the beam power is contained in $2 w_0$. Thus:

- Beam power contained in $(2*gl.w0)

This means that the photons are emitted along a "wire" of $(2*gl.w0) mm diameter and  $lz length. The dimensions of the fluorescent volume are small compared with the solid angle, since the CCD is at a distance d = $dccd mm much larger than t = $lz mm. Therefore, we can approximate that the emission occurs from a single point located in the center of the vase. 

Notice also that for given fluorophore and setup, the product $\sigma\, Q\, C\, L$ is a constant, and thus the slope of the power scan is simply P.

The geometrical acceptance is: 

$A = \frac{1}{2}(1 - \frac{d}{\sqrt{(d^2 + (D/2)^2}})$

where $D$ is the aperture of the CCD and $d$ is the distance to the interaction point.
"""

# ╔═╡ a22bb851-cc1d-4c8b-9a5c-505bffe9d39e
md"### CCD and Fluorophore"

# ╔═╡ 571c0b4d-c554-440c-927f-2996710e2246
md"### Laser Measurement, Ru++ solution C = 2.5E-5 M"

# ╔═╡ 6d946d69-33c0-40d2-8ddf-ebe9b7957ab6
md"## Notebook"

# ╔═╡ 3c0476c6-da9b-44e6-964b-37ea70854ba6
md"- ipes is the intensity in photoelectrons, obtained from the measured data multiplying by the CCD conversion constant ($adcpes pes per adc count)"

# ╔═╡ 9ca51bea-10ac-43e9-8da6-f3d960c634b8
md""" ### Predicted power scan:
""" 

# ╔═╡ 3b1d3646-eef1-4a2b-a6fa-e631307eca40
PP = collect(0.0:10.0:5000.0) * μW

# ╔═╡ 2f74a681-310c-4914-9b25-ee780b513417
md"#### Read Power Scan data"

# ╔═╡ f887dc67-8095-442c-aa3b-9f57e5ec7362
rul = lfi.LabFbi.load_df_from_csv(datadir("solutionLaser"),
	                          "20210603_Ru_Solution_2p5em5_M_Power_ramp.csv", 
	                          lfi.LabFbi.enG)

# ╔═╡ a0f45e02-644c-4b96-ba6e-f2e1c32e60a0
pm = rul.P

# ╔═╡ 817be26b-897e-4996-a452-5c55dd7b45bf
ipes = rul.I[1] *adcpes

# ╔═╡ 6557e731-e074-4d42-b4df-7137be2ba358
md""" ### Measured rate:

- The measured rate at a power of $(rul.P[1]) μW is $(lti.LabTools.to_fstr(ipes,"%5.2g")) Hz
""" 

# ╔═╡ 4a194d15-566d-404a-a022-5347447cf1b6
arul = lfi.LabFbi.load_df_from_csv(datadir("solutionLaser"),
	                          "20210608_Ru_Solution_ALL_M_Power_ramp.csv", 
	                          lfi.LabFbi.enG)

# ╔═╡ 3483e954-f3d0-419e-87c8-1aca0788f36c
md" - Fit to a straight line, y = x1 + x2*x. Then assign the value of x1 to the background ambient light and subtract" 

# ╔═╡ b98f48b8-9445-49b9-94f7-65cf4c72b2dc
ruf = lti.LabTools.lfit(rul.P, rul.I * adcpes)

# ╔═╡ db46f5fd-20d6-4177-9c4c-1db73fab1b94
x1, x2 = coef(ruf)

# ╔═╡ db72332d-9c76-4611-a853-663cbb9a79af
xgI = lfi.LabFbi.Iρ(gl);

# ╔═╡ 984e422b-cc31-44c0-8aec-e3b7ae40b0dc
xr = xgI(1.0mm, 0.0mm) / xgI(0.0mm, 0.0mm)

# ╔═╡ 1081c803-2fef-4141-8540-5d2cdd70fd39
xr2 = xgI(0.5mm, 0.0mm) / xgI(0.0mm, 0.0mm)

# ╔═╡ 6976b4d9-3d80-429c-8b19-4323acf6ae63
begin
	xRmm =collect(-2.0:0.01:2.0)*mm
	xprz0 = lpi.LabPlots.plot_xy(xRmm/mm, xgI.(xRmm,(0.0*mm),)/(Hz*cm^-2), "r (mm)", 
	                "I(r,0) (Hz*cm^-2) ",
	                "I(r,0) for beam");
	xprz1 = lpi.LabPlots.plot_xy(xRmm/mm, xgI.(xRmm,(1.0*mm),)/(Hz*cm^-2), "r (mm)", 
	                "I(r,1mm) (Hz*cm^-2) ",
	                "I(r,1mm) for beam");
	xprz2 = lpi.LabPlots.plot_xy(xRmm/mm, xgI.(xRmm,(5.0*mm),)/(Hz*cm^-2), "r (mm)", 
	                "I(r,5mm) (Hz*cm^-2) ",
	                "I(r,5mm) for beam");
	xprz3 = lpi.LabPlots.plot_xy(xRmm/mm, xgI.(xRmm,(10.0*mm),)/(Hz*cm^-2), "r (mm)", 
	                "I(r,10mm) (Hz*cm^-2) ",
	                "I(r,10mm) for beam");
	void()
end

# ╔═╡ 5c77e6ad-878c-4e90-a3bf-f5fc87b8e4fe
plot(xprz0,xprz1,xprz2,xprz3, layout = (2, 2), legend=false, fmt = :png)

# ╔═╡ 42fa897b-36d6-4568-a9ae-98e5ac359f75
md"## Prediction, Ru++ solution C = 2.5E-5 M"

# ╔═╡ 3ca78d4a-52f2-4b1f-9c2c-d1f3cfbcd14d
md"#### Define CCD"

# ╔═╡ e0bd1a1e-6274-4062-aab3-886ae3146eb3
orca = lfi.LabFbi.ccd()

# ╔═╡ 8e112380-2e85-4c16-a7ff-417f4d2c9d0b
begin
	xL = collect(350.:5.:1000.)
	xEff =orca.(xL)
	porca = plot(xL, xEff, leg=false, lw=2);

	xlabel!("λ (nm)")
	ylabel!("Efficiency")
	title!("Efficiency of CCD ")
	void()
end

# ╔═╡ 0777a0d1-673a-4db3-8eb6-46d9b35e7701
md"#### Read the Ru DF and interpolate function to data"

# ╔═╡ 209747a6-5a77-431c-97ad-f1ad3b185a7f
ru405df = lfi.LabFbi.load_df_from_csv(datadir("fluorimeter/RuIrPy"), "RuppSol.csv", lfi.LabFbi.spG);


# ╔═╡ 527bf9ce-62d6-45f9-a833-9692c6fa36f8
md"#### plot spectrum"

# ╔═╡ b0c0fa14-c303-4709-9c29-f1fad2a38758
md"#### Compute efficiency of CCD integrated over all spectrum"

# ╔═╡ 755e3bb2-9fdc-4ab1-8206-a122fa0e050a
md"### Filters"

# ╔═╡ e48da31e-b52a-4790-bb88-26142c3c32e1
filterd = lfi.LabFbi.get_filters(lfi.LabFbi.fnames, datadir("filters"))

# ╔═╡ 8259abcd-5eab-43de-a345-86624e79e696
fw = [390:420,420:500,420:520,440:540, 390:560, 390:450]

# ╔═╡ eedd1f33-a337-4239-8237-11d83a95fbde
Pft = lpi.LabPlots.plot_filters(lfi.LabFbi.fnames, filterd, fw);

# ╔═╡ 6464b1bb-6c3f-4a6a-af39-69b6718bda96
plot(Pft..., layout = (3, 2), legend=false, fmt = :png)

# ╔═╡ b56f168a-992b-4043-a06a-8620819416a7
cflt =lfi.LabFbi.common_filters(filterd)

# ╔═╡ c14ed4e6-d688-4451-9c80-de16d9e9e204
md" ## Functions"

# ╔═╡ f8fb5bd8-408f-4203-8eec-d94146a96dc3
struct Setup
	C::typeof(1.0M)         # concentration (in M)
	c::typeof(1.0mm^-3)    # concentration (in molecules/mm3)
	t::typeof(1.0mm)       # thickness of vase
	d::typeof(1.0mm)       # diameter of CCD
	D::typeof(1.0mm)       # Distance between CCD and vase
end

# ╔═╡ c13d2440-2944-45e6-b88c-3c74c418cc3e
setup = Setup(cs*M, Cmol_mm3, lzr, dccd*mm, Dccd*mm)

# ╔═╡ bc31b8c8-d5c6-4724-97b8-485075463549
function get_fit(fglm)
	function sline(x)
		return c[1] + c[2] * x
	end
	
	c = coef(fglm)
	return sline
end
	

# ╔═╡ 228a018c-dcae-4804-95d0-a1a6db9dcfd2
prf = get_fit(ruf)

# ╔═╡ 74ee59b9-4f4d-4b9a-b0fa-cfaf1d61cd69
IF = prf.(PP/μW)

# ╔═╡ cc999b93-2a71-4f8e-bae9-4d5b0b6f616d
function fluorescence_rate(gl, fl, t, C)
	ir = lfi.LabFbi.Iρ(gl)
	i0 = ir(0.0mm, 0.0mm)
	return 0.5 * i0 * fl.σ * fl.Q * π * gl.w0^2 * t * C
end

# ╔═╡ 40714bf8-b620-4229-a326-30b8c8a77551
function fluorescence_rate2(gl, fl, t, C)
	p = lfi.LabFbi.n_photons(gl.λ, gl.P)
	return uconvert(Hz, p * fl.σ * fl.Q * t * C)
end

# ╔═╡ 48c61758-bbe1-4832-9159-abe589a59894
function rate_prediction(gl,setup, sset, eff_ccdFilters)
	if sset == "FBI"
		cdf = lfi.LabFbi.load_df_from_csv(datadir("fbi"),
	                          "molar_extinction_coefficient_G1G2.csv", 
	                          lfi.LabFbi.enG) 
		fl = lfi.LabFbi.fbi_fluorophores(cdf, ["g1", "g2"], [325,405], [0.67,0.67])
		fr = fluorescence_rate2(gl, fl.fbi["l405g1"], setup.t, setup.c)

	else
		cdf = lfi.LabFbi.load_df_from_csv(datadir("fbi"),
	                          "molar_extinction_coefficient_RuIr.csv", 
	                          lfi.LabFbi.enG) 
		fl = lfi.LabFbi.iru_fluorophores(cdf, 
	                       ["IrF", "Ir", "Ru", "IrF+","Ir++","Ru++"],
				           [325,405],
				           [Q,Q,Q,Q,Q, Q])
		fr = fluorescence_rate2(gl, fl["l405Ru++"], setup.t, setup.c)
	end
	ga = lfi.LabFbi.geometrical_acceptance(setup.d/mm, setup.D/mm)
	return fr * ga * eff_ccdFilters
	
			
end

# ╔═╡ bdf0eb6f-3c18-408b-9151-89c725237a4f
function rate_prediction_power(pp, setup, sset, eff_ccdFilters)
	GL = [lfi.LabFbi.GLaser(λ, p, w0) for p in pp]
	IFR = [rate_prediction(gl,setup, sset, eff_ccdFilters) for gl in GL]
	return IFR
end

# ╔═╡ 55e57114-347a-434b-be67-2ddb14d43067
function r_θ(θr, z)
	return θr * z * 2.0
end

# ╔═╡ 61f995bc-d2cb-4d93-a05c-b98041c3c564
rxθ = r_θ(gl.θr, lzr) ;

# ╔═╡ 11c323c3-760e-49ea-bf5e-b291a4667075
md"""#### Beam divergence for unfocused beam

- The beam divergence after traversing the thickness of the solution container
(lz = $(lti.LabTools.to_fstr(lzr/mm,"%5.2g")) mm) is very small,
rθ = $(lti.LabTools.to_fstr(rxθ/μm,"%5.2g")) μm. The power at 1 mm from the axis is $(lti.LabTools.to_fstr(xr,"%5.2g")) of the beam power. As shown in the plots below, the power does not change either along z.
"""

# ╔═╡ 45a00153-71a1-4072-a548-507fcc241575
function glaser_md(gl)
	np = lfi.LabFbi.n_photons(gl.λ, gl.P)
	md = " ### Laser \n"
	ll = lti.LabTools.to_fstr(gl.λ/nm,"%5.2g")
	lp = lti.LabTools.to_fstr(gl.P/mW,"%5.2g")
	nnp = lti.LabTools.to_fstr(np/Hz,"%5.2g")
	lw0 = lti.LabTools.to_fstr(gl.w0/mm,"%5.2g")
	lw0 = lti.LabTools.to_fstr(gl.w0/mm,"%5.2g")
	lzr = lti.LabTools.to_fstr(gl.zr/mm,"%5.2g")
	lθ  = lti.LabTools.to_fstr(gl.θr*1e+3,"%5.2g")
	rho = lti.LabTools.to_fstr(gl.ρ0/(Hz* cm^-2),"%5.2g")
	md = string(md, "  - laser λ     = $ll  nm \n")
	md = string(md, "  - laser P     = $lp  mW \n")
	md = string(md, "  - photon rate = $nnp  np  Hz \n")
	md = string(md, "  - laser w0    = $lw0 mm \n")
	md = string(md, "  - laser zr    = $lzr mm \n")
	md = string(md, "  - laser θr    = $lθ  mrad \n")
	md = string(md, "  - photon ρ    = $rho Hz/cm2 \n")
	return md
end

# ╔═╡ b3370280-5321-4a05-816f-ce41acb09952
begin
ginfo = glaser_md(gl)	
Markdown.parse("$ginfo")
end

# ╔═╡ 1e6f38e5-f521-4467-9f8e-e4a62751f8bd
function setup_md(setup, sset)
	md = " ### setup \n"
	ll = lti.LabTools.to_fstr(setup.t/mm,"%5.2g")
	cc = lti.LabTools.to_fstr(setup.C/M,"%5.2g")
	cmol = lti.LabTools.to_fstr(setup.c/mm^-3,"%5.2g")
	dc = lti.LabTools.to_fstr(setup.D/mm,"%5.2g")
	Dc = lti.LabTools.to_fstr(setup.d/mm,"%5.2g")
	md = string(md, "  - molecule            = $sset \n")
	md = string(md, "  - Thickness vase  = $ll  mm \n")
	md = string(md, "  - Concentration   = $cc  M \n")
	md = string(md, "  - Concentration   = $cmol molec/mm^3 \n")
	md = string(md, "  - diameter CCD    = $dc mm  \n")
	md = string(md, "  - distance to CCD = $Dc mm  \n")
	return md
end

# ╔═╡ 73642c7c-e095-4eee-bb07-c86dbcce5c7d
begin
	sinfo = setup_md(setup, sset)
	Markdown.parse("$sinfo")
end

# ╔═╡ d7f6f631-d03f-4c8e-a33e-31ad218fa43c
struct DfInfo
	name::String
	λex ::String
	λi::Float64
	λf::Float64
	CN::Array{String}
	C::Array{typeof(1.0M)}
end

# ╔═╡ 472ba287-5f97-432d-bb3a-077f42365214
function dfinfo(df,nx=2:7)
	function cnames(df, nx)
		return names(df)[nx]
	end
	
	function dflm(df)
		l = df[!, "λ"]
		return l[1], l[end]
	end
	
	function dfc(df, nx=2:7)
		nz = split.(cnames(df,nx),"_")
		cs = [n[2] for n in nz]
		return parse.(Float64, cs) * M
	end
	
	λi, λf = dflm(df)
	cs = dfc(df, nx)
	cn = cnames(df, nx)
	cn1s = split(cn[1],"_")
	lexc = cn1s[3]
	dfname = cn1s[1]
	
	return DfInfo(dfname, lexc, λi, λf, cn, cs)
end

# ╔═╡ a6177074-45db-4a41-8946-30a4a434df22
dfi = dfinfo(ru405df,2:7)

# ╔═╡ 71e7af35-3070-4979-b416-aebf41c996ec
gfru405 = lfi.LabFbi.dftogf(dfi.λi:2.0:dfi.λf, ru405df, dfi.CN[3])

# ╔═╡ 080f125c-1172-427e-8ce6-5c130f28343f
ϵccd(λ) = map(x->orca(x) * gfru405.pdf(x), λ)

# ╔═╡ ad7fbfdd-b1d5-4f7a-813e-0916b868cada
ϵccdFilters(λ) = map(x->cflt(x) * orca(x)* gfru405.pdf(x), λ)

# ╔═╡ e662c24a-2e0b-4e64-b1b0-db0d2084b0bc
begin
	pru405 = lpi.LabPlots.plotdf_xys(ru405df, "λ", 
		                             dfi.CN[3:3], false,  dfi.CN[3:3],
	                                    markercolors, 
					                    "λ(nm)", "I (au)", 
					                    "Ru++", legend=:topright) 
pgfru405 =lpi.LabPlots.plotdf_gfs([gfru405], dfi.λi:2.0:dfi.λf,[dfi.CN[3]], 
	                              pdf=false,
	                              [:red], "λ (nm)", "I (au)", "FBIG1",
		                          legend=:topright)
void()
end

# ╔═╡ da6100ac-e6b6-40a7-b34c-a6d120a7e417
prux = lpi.LabPlots.merge_plots!(pru405, pgfru405)

# ╔═╡ e54c6c15-a386-47b1-a368-4466bb274b8a
plot(porca,prux, layout = (1, 2), legend=false, fmt = :png)

# ╔═╡ f2320190-8caf-4539-8360-e2bff21adb06
lfi.LabFbi.qpdf(gfru405.pdf, dfi.λi, dfi.λf)

# ╔═╡ b967a7d2-83b8-4904-9eb8-2a524de44351
eff_ccd = lfi.LabFbi.qpdf(ϵccd, dfi.λi, dfi.λf)

# ╔═╡ 9462c377-93c8-43bd-9b83-6aed6bdbf379
eff_ccdFilters = lfi.LabFbi.qpdf(ϵccdFilters, dfi.λi, dfi.λf)

# ╔═╡ a96245fd-4ced-4314-b48e-ae2eec93abc8
IP = rate_prediction_power(uconvert.(mW,PP), setup, sset, eff_ccdFilters) 

# ╔═╡ c7ba4361-5a9e-4961-87d8-ad934d0bd31e
begin
pru = scatter(rul.P, rul.I * adcpes, yerror= rul.σI, label="data",
		       shape = :circle, color = :black, markersize = 5, legend=false)
pru2 = plot!(pru, PP/μW, IF, color = :red,
				 lw = 2, label="fit to data",
				 linestyle = :solid)
	
pru3 = plot!(pru2, PP/μW, IP/Hz, color = :blue,
				 lw = 2, label="model",
				 linestyle = :solid, legend=:topleft)
xlabel!("P (μW)")
ylabel!("I (photoelectrons)")
title!("I vs P")

end

# ╔═╡ de85876d-edc2-4a32-baf8-35804dbdc8dc
pip = lpi.LabPlots.plot_xy(PP/μW, IP/Hz, "P (μW)", 
	                "I(Hz)",
	                "Predicted Power scan")

# ╔═╡ 534d40ce-70cc-4f87-ac48-bfe7330c5acb
begin
	pre = lti.LabTools.lfit(PP/μW, IP/Hz)
	xp1, xp2 = coef(pre)
	r2 = x2 / xp2	
	void()
end

# ╔═╡ 8999a1fd-d27c-4b60-acfb-927e0ff6d3f2
md""" - The ratio between the slopes of data and model is $(lti.LabTools.to_fstr(r2,"%5.2g"))"""

# ╔═╡ a99e1894-f7c9-4132-9bf4-ddc91dde26f5
rateccd = rate_prediction(gl,setup, sset, eff_ccdFilters)

# ╔═╡ cdd8ebd3-5491-448e-9346-85a8f3f75d9c
md""" ### Predicted rate:

- The integrated rate expected in the CCD for the above parameters  with a concentration of $cs M is $(lti.LabTools.to_fstr(rateccd/Hz,"%5.2g")) Hz
""" 

# ╔═╡ 5ee1d01d-54e3-40c7-a4dc-e3c40c625b07
function dfinfomd(dfi)
	md =""" ### Fluorophore \n
- df name  = $(dfi.name)
- df λ exc = $(dfi.λex)
- Concentration (in M) = $(lti.LabTools.vect_to_fstr(dfi.C./M, "%s"))
- λmin = $(dfi.λi) nm, λmax = $(dfi.λf) nm
	"""
	return md
end

# ╔═╡ 08c21279-c41c-4d45-b02c-e849cf818cc4
begin
	runfo = dfinfomd(dfi)
	Markdown.parse("$runfo")
end

# ╔═╡ Cell order:
# ╟─50df4350-0159-4a32-a4a3-2b196734853c
# ╠═eb573495-dfeb-471c-a0b2-11b2c76f6b25
# ╠═25b471c6-e7b4-41a9-b7e6-bd316d78b0c7
# ╠═cc9d45d5-81d4-4525-add1-68f5adbbcc42
# ╠═10715061-650c-48a9-87af-d726b3444768
# ╠═f5a9c61a-c216-11eb-1ddc-6f24685e81bc
# ╠═bd59a5b8-8dac-4982-92e3-6d410b69e31c
# ╠═7c296957-7de3-411e-be97-bb3272c29102
# ╠═65c05bfa-a565-41c1-aff7-71f96c7b9cbd
# ╠═189dd02c-5243-49c3-88e4-fdb8f79380c9
# ╠═ac29bcd2-37be-40f8-bbcf-095e3af08a7b
# ╠═a0c6b3d3-a43d-45c3-b511-9a21ce8e2390
# ╟─8c9a2ea3-6a7e-45ae-99a1-6db448064e35
# ╠═223d92a7-6db0-4383-848e-72f592c8627e
# ╠═160eb718-f450-494a-a35b-1d3771d5b35e
# ╠═73f50e92-1c10-4d52-8fe3-f6fcd8e5b42a
# ╠═eb50ea90-be5e-4301-98b6-3871add430ff
# ╟─27c53ff7-e67e-4ce8-bcc1-adaeef8c2129
# ╠═5869c147-6e13-4269-942a-fd84340011ba
# ╠═8077956e-84c6-4565-8647-c8bb197ae65e
# ╟─c449a041-0106-4cb3-a3b0-4b5689c5391d
# ╟─189b1cae-a002-4b1e-9ca2-fd8e4ae00826
# ╟─ba8e38da-c4ff-4638-a362-10ddad2e296e
# ╟─780fb489-af23-4d82-bc6b-e00ed51d0b26
# ╟─05584d03-2987-45bf-b24d-a188a56fe8d5
# ╟─e53540b8-93b6-489f-b533-c34bf680e4a1
# ╟─e43cefef-23da-408a-b205-a7ee4296fb58
# ╟─de34e9c1-48e9-4d30-9c69-f3e5358881ae
# ╟─324cab76-4809-4ccb-b44b-361ae009eb4e
# ╟─cb4e1c63-65c6-41f6-af91-319728046f85
# ╟─450184dd-8c86-47ca-b7b2-0ccdfa6e84ef
# ╠═5b5ef870-3a41-48b1-99b8-6a95e90ad125
# ╟─7c343c65-ac30-4662-b133-9055cb1ed8d5
# ╠═f5f1f76b-9dc6-48c2-8e67-87a565cdf41b
# ╠═82c72876-4a85-403c-9642-96d2505aba32
# ╠═ea8a0d64-4329-4080-80e0-0fa12ec1554b
# ╠═a8be2b0d-26f0-408d-a069-f92cc09d5514
# ╠═d89bc3d6-8642-496b-a3f4-17f5c801512c
# ╟─dac46998-e6a2-41fb-adf7-c1d20d7493ff
# ╟─9602abd4-50cf-47dd-8bcd-450fc3ccf587
# ╟─f9c273a2-b808-4e0e-b949-fe29d05c16d4
# ╟─c3a6a828-b18d-461c-bb17-2b13feaa593e
# ╟─0ac0eeb8-1bde-4a47-874d-f68649c19320
# ╟─6d57e4ea-857a-414e-a848-dd5ef4541a91
# ╟─f3c687c7-03a2-4502-bbab-d6366121d747
# ╟─b0e199e1-c1b9-48ce-bd25-9b48cd5ef9d1
# ╠═89a3cd82-cb35-4f68-89ac-06c50dbf71a9
# ╠═e2b13a80-17bf-4647-b083-f3f9252fa1e4
# ╠═904130b2-f0e8-4f82-8afc-06977fc04801
# ╠═b91be939-f092-4bf3-914c-fbe1b02dc9fd
# ╠═a969225a-1254-48ef-9775-8189ccfb9110
# ╠═371dfb8c-4b37-427a-8be5-9448262ef65d
# ╠═0092e011-1abc-4cdb-a529-d926ae51c81e
# ╠═b3370280-5321-4a05-816f-ce41acb09952
# ╠═c13d2440-2944-45e6-b88c-3c74c418cc3e
# ╠═73642c7c-e095-4eee-bb07-c86dbcce5c7d
# ╠═a22bb851-cc1d-4c8b-9a5c-505bffe9d39e
# ╠═e54c6c15-a386-47b1-a368-4466bb274b8a
# ╠═08c21279-c41c-4d45-b02c-e849cf818cc4
# ╠═cdd8ebd3-5491-448e-9346-85a8f3f75d9c
# ╠═6557e731-e074-4d42-b4df-7137be2ba358
# ╟─571c0b4d-c554-440c-927f-2996710e2246
# ╠═c7ba4361-5a9e-4961-87d8-ad934d0bd31e
# ╠═8999a1fd-d27c-4b60-acfb-927e0ff6d3f2
# ╠═6d946d69-33c0-40d2-8ddf-ebe9b7957ab6
# ╠═a0f45e02-644c-4b96-ba6e-f2e1c32e60a0
# ╠═3c0476c6-da9b-44e6-964b-37ea70854ba6
# ╠═817be26b-897e-4996-a452-5c55dd7b45bf
# ╠═9ca51bea-10ac-43e9-8da6-f3d960c634b8
# ╠═3b1d3646-eef1-4a2b-a6fa-e631307eca40
# ╠═a96245fd-4ced-4314-b48e-ae2eec93abc8
# ╠═de85876d-edc2-4a32-baf8-35804dbdc8dc
# ╟─534d40ce-70cc-4f87-ac48-bfe7330c5acb
# ╠═2f74a681-310c-4914-9b25-ee780b513417
# ╠═f887dc67-8095-442c-aa3b-9f57e5ec7362
# ╠═4a194d15-566d-404a-a022-5347447cf1b6
# ╠═3483e954-f3d0-419e-87c8-1aca0788f36c
# ╠═b98f48b8-9445-49b9-94f7-65cf4c72b2dc
# ╠═228a018c-dcae-4804-95d0-a1a6db9dcfd2
# ╠═74ee59b9-4f4d-4b9a-b0fa-cfaf1d61cd69
# ╠═db46f5fd-20d6-4177-9c4c-1db73fab1b94
# ╟─11c323c3-760e-49ea-bf5e-b291a4667075
# ╠═61f995bc-d2cb-4d93-a05c-b98041c3c564
# ╠═984e422b-cc31-44c0-8aec-e3b7ae40b0dc
# ╠═1081c803-2fef-4141-8540-5d2cdd70fd39
# ╠═db72332d-9c76-4611-a853-663cbb9a79af
# ╟─6976b4d9-3d80-429c-8b19-4323acf6ae63
# ╠═5c77e6ad-878c-4e90-a3bf-f5fc87b8e4fe
# ╠═42fa897b-36d6-4568-a9ae-98e5ac359f75
# ╠═3ca78d4a-52f2-4b1f-9c2c-d1f3cfbcd14d
# ╠═e0bd1a1e-6274-4062-aab3-886ae3146eb3
# ╠═8e112380-2e85-4c16-a7ff-417f4d2c9d0b
# ╠═0777a0d1-673a-4db3-8eb6-46d9b35e7701
# ╠═209747a6-5a77-431c-97ad-f1ad3b185a7f
# ╠═a6177074-45db-4a41-8946-30a4a434df22
# ╠═71e7af35-3070-4979-b416-aebf41c996ec
# ╠═527bf9ce-62d6-45f9-a833-9692c6fa36f8
# ╠═e662c24a-2e0b-4e64-b1b0-db0d2084b0bc
# ╠═da6100ac-e6b6-40a7-b34c-a6d120a7e417
# ╠═b0c0fa14-c303-4709-9c29-f1fad2a38758
# ╠═080f125c-1172-427e-8ce6-5c130f28343f
# ╠═f2320190-8caf-4539-8360-e2bff21adb06
# ╠═b967a7d2-83b8-4904-9eb8-2a524de44351
# ╠═755e3bb2-9fdc-4ab1-8206-a122fa0e050a
# ╠═e48da31e-b52a-4790-bb88-26142c3c32e1
# ╠═8259abcd-5eab-43de-a345-86624e79e696
# ╠═eedd1f33-a337-4239-8237-11d83a95fbde
# ╠═6464b1bb-6c3f-4a6a-af39-69b6718bda96
# ╠═b56f168a-992b-4043-a06a-8620819416a7
# ╠═ad7fbfdd-b1d5-4f7a-813e-0916b868cada
# ╠═9462c377-93c8-43bd-9b83-6aed6bdbf379
# ╠═a99e1894-f7c9-4132-9bf4-ddc91dde26f5
# ╠═c14ed4e6-d688-4451-9c80-de16d9e9e204
# ╠═f8fb5bd8-408f-4203-8eec-d94146a96dc3
# ╠═bc31b8c8-d5c6-4724-97b8-485075463549
# ╠═cc999b93-2a71-4f8e-bae9-4d5b0b6f616d
# ╠═40714bf8-b620-4229-a326-30b8c8a77551
# ╠═48c61758-bbe1-4832-9159-abe589a59894
# ╠═bdf0eb6f-3c18-408b-9151-89c725237a4f
# ╠═55e57114-347a-434b-be67-2ddb14d43067
# ╠═45a00153-71a1-4072-a548-507fcc241575
# ╠═1e6f38e5-f521-4467-9f8e-e4a62751f8bd
# ╠═d7f6f631-d03f-4c8e-a33e-31ad218fa43c
# ╠═472ba287-5f97-432d-bb3a-077f42365214
# ╠═5ee1d01d-54e3-40c7-a4dc-e3c40c625b07

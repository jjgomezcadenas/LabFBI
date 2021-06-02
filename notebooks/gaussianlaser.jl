### A Pluto.jl notebook ###
# v0.14.5

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

# ╔═╡ 2166631f-cd54-4dc7-a122-9616e6939345
import Pkg

# ╔═╡ eb573495-dfeb-471c-a0b2-11b2c76f6b25
Pkg.add("SpecialFunctions")

# ╔═╡ 25b471c6-e7b4-41a9-b7e6-bd316d78b0c7
Pkg.add("HCubature")

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
md"# Calibration of TOPATU with molecular solutions"

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
md"### Gaussian beam formulas

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

Since $\int_0^\infty e^-y dy = 1$, we find:

$P = I_0 \frac{\pi w_0^2}{2}$

or:

$I_0 = \frac{2 P}{\pi w_0^2}$ 

"

# ╔═╡ 406a52d5-f842-4304-9d9a-751b2afe50dd
md"### Define laser parameters at objective lens"

# ╔═╡ 05584d03-2987-45bf-b24d-a188a56fe8d5
md" ##### Set laser wavelength (in nm)"

# ╔═╡ e53540b8-93b6-489f-b533-c34bf680e4a1
@bind ll NumberField(10.0^2:10.0^3; default=405.0)

# ╔═╡ e43cefef-23da-408a-b205-a7ee4296fb58
md" #####  Set laser power (in mW)"

# ╔═╡ de34e9c1-48e9-4d30-9c69-f3e5358881ae
@bind lp NumberField(10.0^-1:10.0^3; default=0.1)

# ╔═╡ 324cab76-4809-4ccb-b44b-361ae009eb4e
md" ##### Set w0 at the objective lens or at the object plane (in mm)"

# ╔═╡ cb4e1c63-65c6-41f6-af91-319728046f85
@bind lw0 NumberField(10.0^-1:10.0^1; default=0.5)

# ╔═╡ 450184dd-8c86-47ca-b7b2-0ccdfa6e84ef
md" ##### Set thickness of the vase (in mm)"

# ╔═╡ 5b5ef870-3a41-48b1-99b8-6a95e90ad125
@bind lz NumberField(1.0:100.0; default=12.43)

# ╔═╡ 7c343c65-ac30-4662-b133-9055cb1ed8d5
md" ##### Set location of the CCD (in mm)"

# ╔═╡ f5f1f76b-9dc6-48c2-8e67-87a565cdf41b
@bind dccd NumberField(10.0:10.0^3; default=550.0)

# ╔═╡ 82c72876-4a85-403c-9642-96d2505aba32
md" ##### Set diameter of the CCD (in mm)"

# ╔═╡ ea8a0d64-4329-4080-80e0-0fa12ec1554b
@bind Dccd NumberField(10.0:20.0; default=13.3)

# ╔═╡ dac46998-e6a2-41fb-adf7-c1d20d7493ff
md""" #### Concentration in solution"""

# ╔═╡ c3a6a828-b18d-461c-bb17-2b13feaa593e
md"- Set units (in M)"

# ╔═╡ 0ac0eeb8-1bde-4a47-874d-f68649c19320
@bind csu NumberField(1.0:10.0; default=1.0)

# ╔═╡ 6d57e4ea-857a-414e-a848-dd5ef4541a91
md"- Set power (in M)"

# ╔═╡ f3c687c7-03a2-4502-bbab-d6366121d747
@bind csp NumberField(1e-5:1e-9; default=1e-5)

# ╔═╡ b0e199e1-c1b9-48ce-bd25-9b48cd5ef9d1
cs = csu * csp 

# ╔═╡ b91be939-f092-4bf3-914c-fbe1b02dc9fd
solfbi = lfi.LabFbi.Solution("FBI G1 C = 1E-5 M", cs*M)

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

# ╔═╡ b2eef01f-774e-43b7-88d7-1b9ddbe10355
gl = lfi.LabFbi.GLaser(λ, P, w0);

# ╔═╡ ba8e38da-c4ff-4638-a362-10ddad2e296e
md"""### Calibration setup

Consider now the problem of calibrating the TOPATU setup using a fluorescent solution of known fluorophores cross section (σ) and quantum yield (Q), with a given concentration (C). The solution is contained in a vial of quartz. The laser beam is collimated to a waist of $(gl.w0) mm and crosses the vial, which has a thickness $(lz) mm. A CCD camera of diameter $Dccd mm, located at a position $dccd mm, with respect to the vial, along the z axis, records the photons emitted by the solution. 

For a waist of $(gl.w0), zr = $(lti.LabTools.to_fstr(gl.zr/mm,"%5.2g")) and thus any point in the vial has a small z compared with zr. Therefore:

$w(z) = w_0$

$I(r, z) = I_0 e^{-2(r/w_0)^2}$

When the beam crosses the vial it iluminates all the fluorophores. The emitted rate of photons is found by integrating the irradiance and multiplying by the cross section, the quantum yield and the concentration:

$N = \sigma\, Q\, C \int_z \int_r I(r) 2 \pi r dr$

The integral in $z$ is trivial and yields the thickness of the vial, $t$. The integral in $r$ can be precisely approximated by the integral to all the space, since the dimensions in r of the vial are much larger than the beam spread. Thus:

$N = I_0\, \sigma\, Q\, C\, t\, 2 \pi \int e^{-2(r/w_0)^2} r dr$

$N =  \frac{ I_0\, \sigma\, Q\, C\, t\, \pi w_0^2}{2}$

For a gaussian beam > 99.9 % of the power is contained in a radius $\rho = 2 w(z)$. In our case $w(z) \sim w_0$, and thus all the beam power is contained in $2 w_0$. Thus:

- Beam power contained in $(2*gl.w0)

This means that the photons are emitted along a "wire" of $(2*gl.w0) mm diameter and  $lz length. The dimensions of the fluorescent volume are small compared with the solid angle, since the CCD is at a distance d = $dccd mm much larger than t = $lz mm. Therefore, we can approximate that the emission occurs from a single point located in the center of the vase. 

The geometrical acceptance is: 

$A = \frac{1}{2}(1 - \frac{d}{\sqrt{(d^2 + (D/2)^2}})$

where $D$ is the aperture of the CCD and $d$ is the distance to the interaction point.
"""

# ╔═╡ db72332d-9c76-4611-a853-663cbb9a79af
xgI = lfi.LabFbi.Iρ(gl);

# ╔═╡ 984e422b-cc31-44c0-8aec-e3b7ae40b0dc
xr = xgI(1.0mm, 0.0mm) / xgI(0.0mm, 0.0mm)

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

# ╔═╡ 48b0db18-f278-4ca1-8ec6-88ad3e0f6055
md"## FBI solution"

# ╔═╡ d0b047ca-9496-4bd6-8e7a-57c9159caec3
cdf = lfi.LabFbi.load_df_from_csv(datadir("fbi"),
	                          "molar_extinction_coefficient_G1G2.csv", 
	                          lfi.LabFbi.enG) 


# ╔═╡ c33167ea-2317-4c66-9a55-33abb9a572f6
fbiFluo = lfi.LabFbi.fbi_fluorophores(cdf, ["g1", "g2"], [325,405], [0.67,0.67])


# ╔═╡ 05a675da-a726-4f60-9af6-33292a0dd207
fbiFluo.fbi["l405g1"]

# ╔═╡ cc999b93-2a71-4f8e-bae9-4d5b0b6f616d
function fluorescence_rate(gl, fl, t, C)
	ir = lfi.LabFbi.Iρ(gl)
	i0 = ir(0.0mm, 0.0mm)
	return 0.5 * i0 * fl.σ * fl.Q * π * gl.w0^2 * t * C
end

# ╔═╡ 586f3ee4-d245-4899-b5cc-4b4bff69450e
fr = fluorescence_rate(gl, fbiFluo.fbi["l405g1"], lzr, Cmol_mm3)

# ╔═╡ a96439b7-7459-4865-b9fc-c20c3330d76c
ga = lfi.LabFbi.geometrical_acceptance(dccd, Dccd)

# ╔═╡ bf586f74-73b5-4926-8349-8b26a3994c81
rate_ccd = fr * ga

# ╔═╡ cdd8ebd3-5491-448e-9346-85a8f3f75d9c
md""" ## Predicted rate

- The integrated rate expected in the CCD is for FBI (G1) with a concentration of $cs M is $(lti.LabTools.to_fstr(rate_ccd/Hz,"%5.2g")) Hz
""" 

# ╔═╡ 28ea436f-2ae9-4604-90ef-a84ef5c49b82
md"## Focusing the GL with a thin lens"

# ╔═╡ 38e92a20-c881-421f-a55b-e63a3d9bc25f
begin
	of = 10.0mm # focal distance
	od = 6.0mm  # lens diameter
	oM = 25.0   # magnification
	obj = lfi.LabFbi.Slobjective(of,od,oM)
end

# ╔═╡ 27508f89-d4b3-4365-b4d4-a38782711ed1
md" #### Objective:

- f     = $(obj.f) 
- d     = $(obj.d)
- M     = $(obj.M)
- NA    = $(obj.NA)
"

# ╔═╡ 5544a3bb-24e1-4a3b-895d-824bd9d3bb44
md"### Focusing the laser
- We assume that the objective is equivalent to a thin lens with focal distance of $(obj.f). The gaussian beam will be focused at the focal, with a new waist computed below."

# ╔═╡ 889708cb-946f-4812-b661-1d59f708e3dc
fgl = lfi.LabFbi.focus_lens_at_beam_waist(gl, obj.f)

# ╔═╡ ed0618e2-4c09-458f-932c-45e06a486f4b
md"""- Notice that the divergence of the beam is much larger than before focusing. If we want to keep 99 % of the beam power we need to consider a cone of 2.0 * θr = $(lti.LabTools.to_fstr(fgl.θr*1e+3*2.0,"%5.2g")) mrad		
"""

# ╔═╡ 9e983917-94bf-4b16-a45c-8d429f74c226
gI = lfi.LabFbi.Iρ(fgl);

# ╔═╡ 9e9cac03-ad49-47c5-9425-c298aa808d3a
begin
	Zmm =collect(-500.0:1.0:500.0)*μm
	pr0z = lpi.LabPlots.plot_xy(Zmm/μm, gI.((0.0*mm),Zmm)/(Hz*cm^-2), "z (μm)", 
	                "I(0,z) (Hz*cm^-2) ",
	                "I(0,z) for beam")
end

# ╔═╡ 4404f872-e5cd-4db6-8f42-dcdc6a188319
ir0z0 = gI(0.0*mm,0.0mm)

# ╔═╡ 1fc2c3b8-46cb-45a7-b370-5697859b0580
ir0z100mu = gI(0.0*mm,200.0μm)

# ╔═╡ cd151e85-5a05-49f5-889b-b90e40d06362
gI(0.0*mm,500.0μm) / gI(0.0*mm,0.0mm)

# ╔═╡ a6c9caa6-6ef8-44ab-972b-0a026ff2c9fb
md"- The beam power at 500 μm is already 1% of the beam power at z=0."

# ╔═╡ f13edc7c-4ef2-4ce4-a1ea-31b82b9d901d
begin
	Rmm =collect(-50:0.1:50.0)*μm
	Rmm2 =collect(-100:0.1:100.0)*μm
	Rmm3 =collect(-1.0:0.01:1.0)*mm
	prz0 = lpi.LabPlots.plot_xy(Rmm/μm, gI.(Rmm,(0.0*μm),)/(Hz*cm^-2), "r (mm)", 
	                "I(r,0) (Hz*cm^-2) ",
	                "I(r,0)");
	prz1 = lpi.LabPlots.plot_xy(Rmm/μm, gI.(Rmm,(100.0*μm),)/(Hz*cm^-2), "r (μm)", 
	                "I(r,z) (Hz*cm^-2) ",
	                "I(r,100μm)");
	prz2 = lpi.LabPlots.plot_xy(Rmm/μm, gI.(Rmm,(250.0*μm),)/(Hz*cm^-2), "r (mm)", 
	                "I(r,z) (Hz*cm^-2) ",
	                "I(r,250μm)");
	prz3 = lpi.LabPlots.plot_xy(Rmm/μm, gI.(Rmm,(500.0*μm),)/(Hz*cm^-2), "r (mm)", 
	                "I(r,z) (Hz*cm^-2) ",
	                "I(r,500μm)");
	void()
end

# ╔═╡ c0b3db15-9baa-46a2-889f-dee45a3a12e2
plot(prz0,prz1,prz2,prz3, layout = (2, 2), legend=false, fmt = :png)

# ╔═╡ bf9118ac-3ee1-4adc-9980-957a769341b5
md"- As the beam moves in z, the spread in r increases, but the beam is always contained in less than 50 μm."

# ╔═╡ 5b139e57-9bcf-4abc-a4a5-420f9fc1048c
gI(0.0μm,500.0μm)

# ╔═╡ fb1ed0f5-2df8-4126-9dee-66d952f11aa0
gI(40.0μm,500.0μm) / gI(0.0μm,500.0μm)

# ╔═╡ f4bfb50a-6a7a-4156-98df-1b8a337959a5
gI(20.0μm,250.0μm) / gI(0.0μm,250.0μm)

# ╔═╡ feef3cda-c04b-4816-9cb2-e0af21faefff
gI(10.0μm,100.0μm) / gI(0.0μm,100.0μm)

# ╔═╡ da7ef3e7-ee1e-4886-a1c8-71f86493e76b
md"- In summary, when the beam is focalised, the beam power is contained in a cube of 1mm in z (500 μm after the focus and 500 μm before the focus) and at most  50 μm diameter. Thus, the acceptance can be computed as if all the photons were emitted from the focal point. On the other hand, the computation of the intensity is harder in this case, since the beam changes both in r and in z. Thus:

$N = \int_{f-z_c}^{f+z_c} \int_0^{r_c} I(r,z) \cdot \sigma \cdot Q dr dz$

Where $z_c \sim 500\mu m$ and $r_c \sim 50 \mu m$.
"

# ╔═╡ c14ed4e6-d688-4451-9c80-de16d9e9e204
md" ## Functions"

# ╔═╡ 3a7ce9a7-72be-4e44-a3ff-1b587e5c068f
md"- example of use of hcubature to integrate a 2d function (first value of the tuple is the integral, second value the error)"

# ╔═╡ 45edc402-7965-49d5-b95c-b1848cbf17c9
hcubature(x -> cos(x[1])*cos(x[2]), [0,0], [1,1])

# ╔═╡ 6202eaa5-4e66-4148-bc9b-356148e39a94
md"### Defining the 2d integral for cubature"

# ╔═╡ df4daa39-e139-4636-a49d-3a18ac87e7ef
function cf2x(x)
	return gI(x[1]*μm,x[2]*μm) /(Hz*cm^-2)
end

# ╔═╡ 57cf38f3-07c8-4a5a-99f9-996c9be53a0a
hcubature(cf2x, (0.0, -500.0,), (50.0, 500.0))

# ╔═╡ 44cd861e-e7c7-4645-b28f-0bc27f238504
function focused_beam_photon_density(r::Float64, z::Float64)
	return hcubature(cf2x, (0, -z), (r, z))[1] * Hz*cm^-2
end

# ╔═╡ 7f4c9b6c-c069-4d26-ad5b-c7c518e970cb
focused_beam_photon_density(50.0, 1000.0)

# ╔═╡ 69e3d574-7a9e-4609-b988-5fb352a0a753
focused_beam_photon_density(500.0, 500.0)

# ╔═╡ 23911f69-b3e9-438d-91f4-2de3a7ab283e
focused_beam_photon_density(500.0, 10000.0)

# ╔═╡ 0112c05c-929f-4b46-bb00-42718e9900cb
focused_beam_photon_density(5.0, 500.0) / focused_beam_photon_density(5.0, 10000.0) 

# ╔═╡ 5271b4c3-2125-41db-8c86-741dfc6096dd
function photon_rate(f::lfi.LabFbi.Fluorophore, 
		             g::lfi.LabFbi.GLaser, 
		             t::typeof(1.0mm))
	
	c = 0.5 * sqrt(π/2.0)
	return c * f.σ * f.Q * t * g.ρ0 * g.w0 
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

# ╔═╡ 35b71070-f4bf-488b-b750-c4331d060caa

rθ = r_θ(fgl.θr, lzr) 

# ╔═╡ 45a00153-71a1-4072-a548-507fcc241575
function glaser_md(gl)
	md = " ### Laser \n"
	ll = lti.LabTools.to_fstr(gl.λ/nm,"%5.2g")
	lp = lti.LabTools.to_fstr(gl.P/mW,"%5.2g")
	lw0 = lti.LabTools.to_fstr(gl.w0/mm,"%5.2g")
	lw0 = lti.LabTools.to_fstr(gl.w0/mm,"%5.2g")
	lzr = lti.LabTools.to_fstr(gl.zr/mm,"%5.2g")
	lθ  = lti.LabTools.to_fstr(gl.θr*1e+3,"%5.2g")
	rho = lti.LabTools.to_fstr(gl.ρ0/(Hz* cm^-2),"%5.2g")
	md = string(md, "  - laser λ     = $ll  nm \n")
	md = string(md, "  - laser P     = $lp  mW \n")
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

# ╔═╡ 20730e56-145d-466c-a82e-f5c0cc685955
begin
	fginfo = glaser_md(fgl)
	Markdown.parse("$fginfo")
end

# ╔═╡ Cell order:
# ╠═50df4350-0159-4a32-a4a3-2b196734853c
# ╠═2166631f-cd54-4dc7-a122-9616e6939345
# ╠═eb573495-dfeb-471c-a0b2-11b2c76f6b25
# ╠═25b471c6-e7b4-41a9-b7e6-bd316d78b0c7
# ╠═f5a9c61a-c216-11eb-1ddc-6f24685e81bc
# ╠═7c296957-7de3-411e-be97-bb3272c29102
# ╠═65c05bfa-a565-41c1-aff7-71f96c7b9cbd
# ╠═189dd02c-5243-49c3-88e4-fdb8f79380c9
# ╠═ac29bcd2-37be-40f8-bbcf-095e3af08a7b
# ╠═a0c6b3d3-a43d-45c3-b511-9a21ce8e2390
# ╠═8c9a2ea3-6a7e-45ae-99a1-6db448064e35
# ╠═223d92a7-6db0-4383-848e-72f592c8627e
# ╠═160eb718-f450-494a-a35b-1d3771d5b35e
# ╠═27c53ff7-e67e-4ce8-bcc1-adaeef8c2129
# ╠═5869c147-6e13-4269-942a-fd84340011ba
# ╠═8077956e-84c6-4565-8647-c8bb197ae65e
# ╠═c449a041-0106-4cb3-a3b0-4b5689c5391d
# ╠═189b1cae-a002-4b1e-9ca2-fd8e4ae00826
# ╠═ba8e38da-c4ff-4638-a362-10ddad2e296e
# ╠═406a52d5-f842-4304-9d9a-751b2afe50dd
# ╠═05584d03-2987-45bf-b24d-a188a56fe8d5
# ╠═e53540b8-93b6-489f-b533-c34bf680e4a1
# ╠═e43cefef-23da-408a-b205-a7ee4296fb58
# ╠═de34e9c1-48e9-4d30-9c69-f3e5358881ae
# ╠═324cab76-4809-4ccb-b44b-361ae009eb4e
# ╠═cb4e1c63-65c6-41f6-af91-319728046f85
# ╠═450184dd-8c86-47ca-b7b2-0ccdfa6e84ef
# ╠═5b5ef870-3a41-48b1-99b8-6a95e90ad125
# ╠═7c343c65-ac30-4662-b133-9055cb1ed8d5
# ╠═f5f1f76b-9dc6-48c2-8e67-87a565cdf41b
# ╠═82c72876-4a85-403c-9642-96d2505aba32
# ╠═ea8a0d64-4329-4080-80e0-0fa12ec1554b
# ╠═dac46998-e6a2-41fb-adf7-c1d20d7493ff
# ╠═c3a6a828-b18d-461c-bb17-2b13feaa593e
# ╠═0ac0eeb8-1bde-4a47-874d-f68649c19320
# ╠═6d57e4ea-857a-414e-a848-dd5ef4541a91
# ╠═f3c687c7-03a2-4502-bbab-d6366121d747
# ╠═b0e199e1-c1b9-48ce-bd25-9b48cd5ef9d1
# ╠═b91be939-f092-4bf3-914c-fbe1b02dc9fd
# ╠═a969225a-1254-48ef-9775-8189ccfb9110
# ╠═b2eef01f-774e-43b7-88d7-1b9ddbe10355
# ╠═371dfb8c-4b37-427a-8be5-9448262ef65d
# ╠═b3370280-5321-4a05-816f-ce41acb09952
# ╠═11c323c3-760e-49ea-bf5e-b291a4667075
# ╠═61f995bc-d2cb-4d93-a05c-b98041c3c564
# ╠═984e422b-cc31-44c0-8aec-e3b7ae40b0dc
# ╠═db72332d-9c76-4611-a853-663cbb9a79af
# ╠═6976b4d9-3d80-429c-8b19-4323acf6ae63
# ╠═5c77e6ad-878c-4e90-a3bf-f5fc87b8e4fe
# ╠═48b0db18-f278-4ca1-8ec6-88ad3e0f6055
# ╠═d0b047ca-9496-4bd6-8e7a-57c9159caec3
# ╠═c33167ea-2317-4c66-9a55-33abb9a572f6
# ╠═05a675da-a726-4f60-9af6-33292a0dd207
# ╠═cc999b93-2a71-4f8e-bae9-4d5b0b6f616d
# ╠═586f3ee4-d245-4899-b5cc-4b4bff69450e
# ╠═a96439b7-7459-4865-b9fc-c20c3330d76c
# ╠═bf586f74-73b5-4926-8349-8b26a3994c81
# ╠═cdd8ebd3-5491-448e-9346-85a8f3f75d9c
# ╠═28ea436f-2ae9-4604-90ef-a84ef5c49b82
# ╠═38e92a20-c881-421f-a55b-e63a3d9bc25f
# ╠═27508f89-d4b3-4365-b4d4-a38782711ed1
# ╠═5544a3bb-24e1-4a3b-895d-824bd9d3bb44
# ╠═889708cb-946f-4812-b661-1d59f708e3dc
# ╠═20730e56-145d-466c-a82e-f5c0cc685955
# ╠═ed0618e2-4c09-458f-932c-45e06a486f4b
# ╠═35b71070-f4bf-488b-b750-c4331d060caa
# ╠═9e983917-94bf-4b16-a45c-8d429f74c226
# ╠═9e9cac03-ad49-47c5-9425-c298aa808d3a
# ╠═4404f872-e5cd-4db6-8f42-dcdc6a188319
# ╠═1fc2c3b8-46cb-45a7-b370-5697859b0580
# ╠═cd151e85-5a05-49f5-889b-b90e40d06362
# ╠═a6c9caa6-6ef8-44ab-972b-0a026ff2c9fb
# ╠═f13edc7c-4ef2-4ce4-a1ea-31b82b9d901d
# ╠═c0b3db15-9baa-46a2-889f-dee45a3a12e2
# ╠═bf9118ac-3ee1-4adc-9980-957a769341b5
# ╠═5b139e57-9bcf-4abc-a4a5-420f9fc1048c
# ╠═fb1ed0f5-2df8-4126-9dee-66d952f11aa0
# ╠═f4bfb50a-6a7a-4156-98df-1b8a337959a5
# ╠═feef3cda-c04b-4816-9cb2-e0af21faefff
# ╠═da7ef3e7-ee1e-4886-a1c8-71f86493e76b
# ╠═c14ed4e6-d688-4451-9c80-de16d9e9e204
# ╠═3a7ce9a7-72be-4e44-a3ff-1b587e5c068f
# ╠═45edc402-7965-49d5-b95c-b1848cbf17c9
# ╠═6202eaa5-4e66-4148-bc9b-356148e39a94
# ╠═df4daa39-e139-4636-a49d-3a18ac87e7ef
# ╠═57cf38f3-07c8-4a5a-99f9-996c9be53a0a
# ╠═44cd861e-e7c7-4645-b28f-0bc27f238504
# ╠═7f4c9b6c-c069-4d26-ad5b-c7c518e970cb
# ╠═69e3d574-7a9e-4609-b988-5fb352a0a753
# ╠═23911f69-b3e9-438d-91f4-2de3a7ab283e
# ╠═0112c05c-929f-4b46-bb00-42718e9900cb
# ╠═5271b4c3-2125-41db-8c86-741dfc6096dd
# ╠═55e57114-347a-434b-be67-2ddb14d43067
# ╠═45a00153-71a1-4072-a548-507fcc241575

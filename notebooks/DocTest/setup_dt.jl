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
end

# ╔═╡ 955d45a1-9d25-4a23-86c2-cc279d14e710
using DrWatson

# ╔═╡ c60716c4-108c-4c1e-85c1-024520065952


# ╔═╡ 0447a463-0ceb-49e9-a565-a47d9e881455
md"# Test and documentation of functions in module `setup.jl`"


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

# ╔═╡ 7d35cba2-5ea0-4dd4-9fc9-2e4e7d1f8699
setup = ingredients(srcdir("setup.jl"))

# ╔═╡ e96c8a7c-8f79-4608-ab29-444950b65f12
plt   = ingredients(srcdir("plotDF.jl"))

# ╔═╡ 583c4139-53a9-4c85-a2b0-c09285e57da9
dff   = ingredients(srcdir("dffunctions.jl"))

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

# ╔═╡ 4c9c1952-9a1b-4bfa-b665-a05ad399e112
md"## Notebook"

# ╔═╡ 697c1d8a-2812-4e67-9954-71d8739f9186
md"""
	struct Fov

Represent a field of view

\# Fields
- `d::Unitful.Length`  : diameter of Fov
- `z::Unitful.Length`  : thickness
- `a::Unitful.Area`    : area (computed)
- `v::Unitful.Volume`  : volume (computed)

"""


# ╔═╡ 3bf2b831-219e-48f4-aed9-903c44a8de0e
fov = setup.Fov(1mm, 1mm)

# ╔═╡ e310b391-2696-4593-bef4-68ba63663373
@test fov.a ≈ π * (fov.d/2.0)^2

# ╔═╡ 4b691190-4380-4bb7-95b5-34e747e97f75
@test fov.v ≈ fov.a * fov.z

# ╔═╡ f3905891-0212-4aaa-bf38-f0f5334fbb9e
md"""
	struct Laser

Simplest representation of a laser

\# Fields
- `λ::typeof(1.0nm)`  : Laser wavelength
- `P::typeof(1.0mW)`  : Power

"""

# ╔═╡ f5c211f0-fec6-4504-b35e-67926207676d
l405 = setup.Laser(405nm, 1mW)

# ╔═╡ 9bc6a36b-81f3-48b4-956e-e6c2b5614a26
md"""
	struct Objective

Simple representation of a microscope objective

\# Fields
- `name::String` : identifies the objective
- `NA::Float64`  : Numerical aperture
- `M::Float64`   : Magnification

"""

# ╔═╡ 91586ae0-fc7c-432a-8339-576be9a36267
obj = setup.Objective("high NA objective", 0.9, 100.0)

# ╔═╡ 8a860981-f307-459e-8f71-f68e86a178ad
md"""
	struct GaussianLaser

Representation of a Gaussian laser

\# Fields
- `laser::Laser`       : A laser type
- `obj::Objective`     : An objective type
- `w0::typeof(1.0nm)`  : Waist of laser at focusing point
- `zr::typeof(1.0nm)`  : z of laser at focusing point

"""

# ╔═╡ b36913a0-f326-4c11-be0e-c6aab0c929ad
gl = setup.GaussianLaser(l405,obj)

# ╔═╡ 7a46771e-56f2-49f2-bbb4-d50ece6d9396
@test gl.w0/gl.zr ≈ gl.obj.NA  

# ╔═╡ ed74c7b1-0dc3-4017-91e0-2e7284a6e670
md"""
	diffraction_limit(l::Laser, obj:: Objective)

Return the diameter of diffractive spot for a laser l
focused with an objective obj

\# Fields

- `l::Laser`         : A laser
- `obj:: Objective`  : An objective
"""

# ╔═╡ 66c91053-5c88-4de7-ad51-38cbfd80cf51
dl = setup.diffraction_limit(l405, obj)

# ╔═╡ db83f482-54d4-4bc1-9f87-dc92c9fa0af5
md"""
	diffraction_limit(gl::GaussianLaser)

Return the diameter of diffractive spot for a gaussian laser

\# Fields

- `gl::GaussianLaser`  : A gaussian laser
"""

# ╔═╡ 69e7a83b-a233-4e0f-a1d5-024d879734ea
dlgl = setup.diffraction_limit(gl)

# ╔═╡ 78619a76-6875-4db7-b8d8-d39d63c8fae4
@test setup.diffraction_limit(l405, obj) ≈ setup.diffraction_limit(gl)

# ╔═╡ a23afa19-251e-475f-9858-7a4cf4edce7a
md"""
	photon_energy(λ::Unitful.Length)

Given wavelength of photon return its energy.
\# Fields

- `λ::Unitful.Length`  : Photon wavelength

"""

# ╔═╡ ace197ec-f75f-4a61-af07-c4fdff2a66a3
md"""
	delivered_energy(laser::Laser, t::Unitful.Time)

Delivered energy of a laser in a given time.
\# Fields

- `laser::Laser`     : Laser
- `t::Unitful.Time`  : Time in which target is illuminated

"""

# ╔═╡ 553fe29c-9d4d-4c17-8d96-c034ffadb527
pe405nm = setup.photon_energy(405nm)

# ╔═╡ a7924698-401d-458b-a93f-8db833cbef6e
begin
	Lx = collect(250:800) * nm
	Ex = setup.photon_energy.(Lx)
end

# ╔═╡ e5a86e45-bb67-4fc7-8d18-c5031100976d
plt.plot_xy(Lx/nm, Ex/eV, 
	        "λ (nm)", "E (eV)", "Photon energy as a function of wavelength")

# ╔═╡ adf5d46f-8c72-4d4e-a02e-f75c9041df36
@test uconvert(mJ, setup.delivered_energy(l405, 1s)) ≈ 1mJ

# ╔═╡ 56dae56b-7cd6-4714-b7e7-0ae8c4102446
md"""
	n_photons(laser::Laser)

Rate of photons (number of photons per unit time) produced by a laser
\# Fields

- `laser::Laser`     : Laser
- `t::Unitful.Time`  : Time in which target is illuminated

"""

# ╔═╡ c3f221f6-bfc7-433f-b1b7-b4d7ad8d77cc
md"""
	n_photons(λ::Unitful.Length, p::Unitful.Power)

Rate of photons (number of photons per unit time) corresponding to a wavelength
λ and a power P

\# Fields

- `λ::Unitful.Length` : photon wavelength
- `p::Unitful.Power`  : Power

"""

# ╔═╡ e2a45f4d-2dc4-4b48-8592-ba5eb55d1a4e
@test setup.n_photons(405nm, 1mW) ≈ setup.n_photons(l405)

# ╔═╡ ac3cb632-15f8-470e-895c-30cfa6ff8937
Nl = setup.n_photons.(Lx, (1mW,))

# ╔═╡ eff117b4-c2d1-4bd5-ab2c-80a9d0a4fb4f
plt.plot_xy(Lx/nm, Nl/Hz, 
	        "λ (nm)", "N", "Number of photons as a function of wavelength")

# ╔═╡ 4d1a4f4e-d43c-409b-bc2a-80765a7a59fe
logrange(x1, x2, n) = (10^y for y in range(log10(x1), log10(x2), length=n))


# ╔═╡ 8af87fd1-f8e5-4788-ab1a-5f677904affc
begin
	Px = collect(logrange(1, 10^6, 100)) *μW
	Np = setup.n_photons.((405nm,), Px)
end

# ╔═╡ 07c584b4-771a-4692-9330-e3e219845890
plt.plot_xy(Px/μW, Np/Hz, 
	        "P (μW)", "N", "Number of photons as a function of Power")

# ╔═╡ 5819a890-4011-4da9-bd99-fec65761f2e9
md"""
	photon_density(λ::Unitful.Length, p::Unitful.Power, a::Unitful.Area)

number of photons per unit time per unit area

\# Fields

- `λ::Unitful.Length` : photon wavelength
- `p::Unitful.Power`  : Power
- `a::Unitful.Area`   : Area

"""

# ╔═╡ b08c5765-2717-4ba4-9a99-eb2489204ded
md"""
	photon_density(l::Laser, fov::FoV)

number of photons per unit time per unit area, in a Fov illuminated by a laser

\# Fields

- `laser::Laser` : Laser
- `fov::Fov`     : Field of view

"""

# ╔═╡ 5f5e91cf-f83e-4dbe-8381-8e482e95cf17
begin
	Ax = collect(logrange(1, 10^-6, 100)) *mm^2
	Pd = setup.photon_density.((405nm,), (1mW,), Ax)
end

# ╔═╡ ec39575c-8c2c-43b4-bc19-d3ef4c8d1936
plt.loglog_xy(Ax/mm^2, Pd*mm^2/Hz, 
	        L"A(mm^2)", L"I(Hz/mm^2)", "Photon density as a function of area")

# ╔═╡ 0aacefd0-33b2-4546-b946-33d7d42a0ec1
@test setup.photon_density(l405, fov) ≈ setup.photon_density(405nm, 1mW, fov.a)

# ╔═╡ db23c467-4abb-4e33-8750-5be31d2f7251
md"""
	function w(gl::GaussianLaser, z::Unitful.Length)

Waist of a laser at length z

\# Fields

- `gl::GaussianLaser`  : A gaussian laser
- `z::Unitful.Length`  : z distance from focusing point
"""

# ╔═╡ f6f2ad73-e862-4295-9b90-a78e33e4e131
begin
	Z = collect(logrange(1, 10^2, 100)) * μm
	Wz = setup.w.((gl,), Z)
end

# ╔═╡ 15fac2a5-761d-4a33-b714-010e3fdd5880
plt.plot_xy(Z/μm, uconvert.(μm, Wz)/μm, "z(μm)", "I(μm)", "wz vs z for a gaussian beam")

# ╔═╡ e5095cff-55eb-4433-a05d-3565c89b7882
md"""
	function gf(gl::GaussianLaser, z::Unitful.Length, r::Unitful.Length)

Gaussian distribution of a gaussian beam

\# Fields

- `gl::GaussianLaser`  : A gaussian laser
- `z::Unitful.Length`  : z distance from focusing point
- `r::Unitful.Length`  : r distance from focusing point
"""

# ╔═╡ a78687d6-b027-4c20-97cf-00c95ae8399a
begin
	Zr = collect(logrange(1, 500, 100)) 
	R = Zr
end

# ╔═╡ f5d23440-d4e5-4f93-be19-39576ab66c4b
fgy(x,y) = setup.gf(gl, x*nm, y*nm)

# ╔═╡ 61c2b72b-1baf-4ceb-a9a7-3db0cdff9909
contour(Zr,R,fgy,fill=true)

# ╔═╡ a52c7054-75f6-43c2-a47c-481c6c73ad64
md"""
	function gI(gl::GaussianLaser, z::Unitful.Length, r::Unitful.Length)

Intensity of a gaussian beam

\# Fields

- `gl::GaussianLaser`  : A gaussian laser
- `z::Unitful.Length`  : z distance from focusing point
- `r::Unitful.Length`  : r distance from focusing point
"""

# ╔═╡ e45153e3-47dc-4865-a983-ede17ba5bea6
md"### Gaussian beam formulas

$g(z,r) = e^{-2(r/w(z)^2}$

$I(z,r) = (\frac{w_0}{w(z)})^2 g(z,r)$

$w(z) = w_0 \sqrt{cw(z)}$

$cw(z) =1 + (\frac{z}{z_r})^2$

The formulas verify the following limits:

$I(0, r) = g(0,r)$

$I(z, 0) = g(z,0)^2$

Both conditions are used to provide tests
"

# ╔═╡ 2bf791dc-2b6f-46cb-aee2-10536cc1ff3e
gz0(r) = setup.gf(gl, 0*nm, r*nm)

# ╔═╡ 52d80fab-06a6-42f0-bcd9-89f611edb7d9
gr0(z) = setup.gf(gl, z*nm, 0*nm)

# ╔═╡ c9eba461-3b98-4f3b-a9bd-fe99cd9066f9
pgz0 = plt.plot_xy(R, gz0.(R), "R(nm)", "GR(z0)", "Beam profile in R at Z= 0");

# ╔═╡ a8478687-bb69-4df6-8452-60ac1819c327
pgr0 = plt.plot_xy(R, gr0.(R), "Z(nm)", "GZ(r0)", "Beam profile in Z at R= 0");

# ╔═╡ 7387cc0c-d049-4acd-9b26-701dfeac00a2
plot(pgz0,pgr0, layout = (1, 2), legend=false, fmt = :png)

# ╔═╡ cda298c3-1471-409e-8b83-7f239bbdb52d
gir0(z) = setup.gI(gl, z*nm, 0*nm)

# ╔═╡ d0edb6d9-37ea-43a8-841b-73db6b948d4c
giz0(r) = setup.gI(gl, 0*nm, r*nm)

# ╔═╡ bb93ca5a-0032-40e0-842e-b85964af6633
piz0 = plt.plot_xy(R, giz0.(R), "R(nm)", "IG(z0)", "Beam profile in R at Z= 0");

# ╔═╡ 931ce485-a0f2-4cca-8fc7-8103962256fe
pir0 = plt.plot_xy(R, gir0.(R), "Z(nm)", "IG(r0)", "Beam profile in Z at R= 0");

# ╔═╡ b32efcc4-f0dc-402c-9f90-351fe5ff6031
plot(piz0,pir0, layout = (1, 2), legend=false, fmt = :png)

# ╔═╡ c4abbd03-670b-429f-8903-6b5802e2aeff
@test dff.qpdf(giz0, 0.0, 100.0) ≈ dff.qpdf(gz0, 0.0, 100.0)

# ╔═╡ 84a447c0-7155-481f-8744-2e033270e0ac
rg0(x) = gir0(x) / gr0(x)

# ╔═╡ 6441c23c-e07a-4044-8c5c-b44bf1cd4fc5
@test dff.qpdf(rg0, 0.0, 100.0) ≈ dff.qpdf(gr0, 0.0, 100.0)

# ╔═╡ 92c4c3c5-1a64-44c6-9429-142445a2ccbd
pr0 = plt.plot_xy(R, rg0.(R), "Z(nm)", "I/G", "Ratio /IG")

# ╔═╡ ee05ed95-a708-487a-8cd1-92ea593c7c7f
plt.merge_plots!(pgr0, pr0) 

# ╔═╡ 46a6cf57-0c1f-4287-82bb-fda7b4923cc5
md"""
	geometrical_acceptance(d::Float64, D::Float64)

Compute the fraction of photons that make it through an iris
of diameter D located at a distance d from the emission point.

\# Fields

- `d::Float64`   : distance between emission point and iris
- `D::Float64`   : Diameter of iris

"""

# ╔═╡ 44aa9670-0fcb-4825-9d7e-5d3d336ccad6
begin
	Dx = collect(logrange(1, 10^2, 100))
	GA = setup.geometrical_acceptance.((1.0,), Dx)
end

# ╔═╡ 65ced90b-ca33-4f1b-953c-7d313b777e3c
plt.plot_xy(Dx, GA, 
	        "D (in units of d)", "Acceptance", "Geometrical acceptance D/d")

# ╔═╡ 9fa534e0-a798-4553-9253-7b95bc2f3310
begin
	dx = collect(logrange(1, 10^2, 100))
	GAd = setup.geometrical_acceptance.(Dx, (1.0,))
end

# ╔═╡ 2ce5076c-270c-43c2-a867-89fa15ed0119
plt.loglog_xy(Dx, GAd, 
	        "d (in units of D)", "Acceptance", "Geometrical acceptance D/d")

# ╔═╡ 9f3fdaec-f032-4bfd-b18a-e001eefcee80
md"""
	struct Objective

Simple representation of a microscope objective

\# Fields
- `name::String` : identifies the objective
- `NA::Float64`  : Numerical aperture
- `M::Float64`   : Magnification

"""

# ╔═╡ 519faa69-72aa-4bf3-a03c-e8215c099992
ohna = setup.Objective("High NA", 0.9, 100.0)

# ╔═╡ 4ebb0e8c-9b1d-43c2-a3f1-eadee4b72acb
olna = setup.Objective("Low NA", 0.5, 100.0)

# ╔═╡ f8987592-282e-4b54-af72-a3e36ea1692b
md"""
	transmission(objective::Objective)

Compute the transmission of an objective (depends only of NA).

\# Fields

- `objective::Objective` : Objective

"""

# ╔═╡ 116c46d9-f02d-4267-a643-14750b121029
md"""
	transmission(A::Float64)

Compute the transmission as a function of NA.

\# Fields

- `A::Float64` : Numerical acceptance (NA)

"""

# ╔═╡ 939d8509-b57a-4e62-b81e-efeef137ef4c
setup.transmission(ohna)

# ╔═╡ e005d049-4750-46f2-9f11-79e8e15a1a31
setup.transmission(olna)

# ╔═╡ 05ab7532-6efe-403a-8a22-0cbe93035a32
@test setup.transmission(0.9) ≈ setup.transmission(ohna)

# ╔═╡ eca4fcfd-0997-4bef-bdea-bf4b4844bdd2
begin
	Nx = collect(0.1:0.005:1.0)
	Tx = setup.transmission.(Nx)
end

# ╔═╡ 6d39c50b-e465-4ac3-8903-10fc1207c8a8
plt.plot_xy(Nx, Tx, 
	        "NA", "T", "Transmission as a function of Numerical Aperture")

# ╔═╡ 14c6aec0-37c8-417b-a874-5fbe9e9d291f
md"""
	ccd(lmin::Float64=350.0, lmax::Float64=1000.0)

Return the efficiency of a CCD as a function of wavelength.

\# Fields

- `lmin::Float64=350.0` : Minimum wavelength for which efficiency is defined
- `lmax::Float64=350.0` : Maximum wavelength for which efficiency is defined

"""

# ╔═╡ 61137405-f013-4807-87d5-b10c1c3ef99f
orca = setup.ccd()

# ╔═╡ c744ee52-0b9d-44ba-9077-fcacd32a1ed9
Ly = collect(350.0:1000.0)

# ╔═╡ 465a72e1-0494-492d-b093-d617ea48a12e
Eo = orca.(Ly)

# ╔═╡ 99b59d73-2303-4f53-8b12-c40e6726a302
plt.plot_xy(Ly, Eo, 
	        "λ(nm)", "ϵ", "CCD efficiency as a function of λ")

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
# ╠═7d35cba2-5ea0-4dd4-9fc9-2e4e7d1f8699
# ╠═e96c8a7c-8f79-4608-ab29-444950b65f12
# ╠═583c4139-53a9-4c85-a2b0-c09285e57da9
# ╠═14f34300-f974-432d-8a53-91da2a0dc298
# ╠═4c9c1952-9a1b-4bfa-b665-a05ad399e112
# ╟─697c1d8a-2812-4e67-9954-71d8739f9186
# ╠═3bf2b831-219e-48f4-aed9-903c44a8de0e
# ╠═e310b391-2696-4593-bef4-68ba63663373
# ╠═4b691190-4380-4bb7-95b5-34e747e97f75
# ╟─f3905891-0212-4aaa-bf38-f0f5334fbb9e
# ╠═f5c211f0-fec6-4504-b35e-67926207676d
# ╟─9bc6a36b-81f3-48b4-956e-e6c2b5614a26
# ╠═91586ae0-fc7c-432a-8339-576be9a36267
# ╟─8a860981-f307-459e-8f71-f68e86a178ad
# ╠═b36913a0-f326-4c11-be0e-c6aab0c929ad
# ╠═7a46771e-56f2-49f2-bbb4-d50ece6d9396
# ╠═ed74c7b1-0dc3-4017-91e0-2e7284a6e670
# ╠═66c91053-5c88-4de7-ad51-38cbfd80cf51
# ╠═db83f482-54d4-4bc1-9f87-dc92c9fa0af5
# ╠═69e7a83b-a233-4e0f-a1d5-024d879734ea
# ╠═78619a76-6875-4db7-b8d8-d39d63c8fae4
# ╟─a23afa19-251e-475f-9858-7a4cf4edce7a
# ╟─ace197ec-f75f-4a61-af07-c4fdff2a66a3
# ╠═553fe29c-9d4d-4c17-8d96-c034ffadb527
# ╠═a7924698-401d-458b-a93f-8db833cbef6e
# ╠═e5a86e45-bb67-4fc7-8d18-c5031100976d
# ╠═adf5d46f-8c72-4d4e-a02e-f75c9041df36
# ╟─56dae56b-7cd6-4714-b7e7-0ae8c4102446
# ╟─c3f221f6-bfc7-433f-b1b7-b4d7ad8d77cc
# ╠═e2a45f4d-2dc4-4b48-8592-ba5eb55d1a4e
# ╠═ac3cb632-15f8-470e-895c-30cfa6ff8937
# ╠═eff117b4-c2d1-4bd5-ab2c-80a9d0a4fb4f
# ╠═4d1a4f4e-d43c-409b-bc2a-80765a7a59fe
# ╠═8af87fd1-f8e5-4788-ab1a-5f677904affc
# ╠═07c584b4-771a-4692-9330-e3e219845890
# ╟─5819a890-4011-4da9-bd99-fec65761f2e9
# ╟─b08c5765-2717-4ba4-9a99-eb2489204ded
# ╠═5f5e91cf-f83e-4dbe-8381-8e482e95cf17
# ╠═ec39575c-8c2c-43b4-bc19-d3ef4c8d1936
# ╠═0aacefd0-33b2-4546-b946-33d7d42a0ec1
# ╟─db23c467-4abb-4e33-8750-5be31d2f7251
# ╠═f6f2ad73-e862-4295-9b90-a78e33e4e131
# ╠═15fac2a5-761d-4a33-b714-010e3fdd5880
# ╟─e5095cff-55eb-4433-a05d-3565c89b7882
# ╠═a78687d6-b027-4c20-97cf-00c95ae8399a
# ╠═f5d23440-d4e5-4f93-be19-39576ab66c4b
# ╠═61c2b72b-1baf-4ceb-a9a7-3db0cdff9909
# ╟─a52c7054-75f6-43c2-a47c-481c6c73ad64
# ╠═e45153e3-47dc-4865-a983-ede17ba5bea6
# ╠═2bf791dc-2b6f-46cb-aee2-10536cc1ff3e
# ╠═52d80fab-06a6-42f0-bcd9-89f611edb7d9
# ╠═c9eba461-3b98-4f3b-a9bd-fe99cd9066f9
# ╠═a8478687-bb69-4df6-8452-60ac1819c327
# ╠═7387cc0c-d049-4acd-9b26-701dfeac00a2
# ╠═cda298c3-1471-409e-8b83-7f239bbdb52d
# ╠═d0edb6d9-37ea-43a8-841b-73db6b948d4c
# ╠═bb93ca5a-0032-40e0-842e-b85964af6633
# ╠═931ce485-a0f2-4cca-8fc7-8103962256fe
# ╠═b32efcc4-f0dc-402c-9f90-351fe5ff6031
# ╠═c4abbd03-670b-429f-8903-6b5802e2aeff
# ╠═84a447c0-7155-481f-8744-2e033270e0ac
# ╠═6441c23c-e07a-4044-8c5c-b44bf1cd4fc5
# ╠═92c4c3c5-1a64-44c6-9429-142445a2ccbd
# ╠═ee05ed95-a708-487a-8cd1-92ea593c7c7f
# ╟─46a6cf57-0c1f-4287-82bb-fda7b4923cc5
# ╠═44aa9670-0fcb-4825-9d7e-5d3d336ccad6
# ╠═65ced90b-ca33-4f1b-953c-7d313b777e3c
# ╠═9fa534e0-a798-4553-9253-7b95bc2f3310
# ╠═2ce5076c-270c-43c2-a867-89fa15ed0119
# ╟─9f3fdaec-f032-4bfd-b18a-e001eefcee80
# ╠═519faa69-72aa-4bf3-a03c-e8215c099992
# ╠═4ebb0e8c-9b1d-43c2-a3f1-eadee4b72acb
# ╟─f8987592-282e-4b54-af72-a3e36ea1692b
# ╟─116c46d9-f02d-4267-a643-14750b121029
# ╠═939d8509-b57a-4e62-b81e-efeef137ef4c
# ╠═e005d049-4750-46f2-9f11-79e8e15a1a31
# ╠═05ab7532-6efe-403a-8a22-0cbe93035a32
# ╠═eca4fcfd-0997-4bef-bdea-bf4b4844bdd2
# ╠═6d39c50b-e465-4ac3-8903-10fc1207c8a8
# ╟─14c6aec0-37c8-417b-a874-5fbe9e9d291f
# ╠═61137405-f013-4807-87d5-b10c1c3ef99f
# ╠═c744ee52-0b9d-44ba-9077-fcacd32a1ed9
# ╠═465a72e1-0494-492d-b093-d617ea48a12e
# ╠═99b59d73-2303-4f53-8b12-c40e6726a302

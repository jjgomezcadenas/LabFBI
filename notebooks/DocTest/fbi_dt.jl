### A Pluto.jl notebook ###
# v0.14.3

using Markdown
using InteractiveUtils

# ╔═╡ 9e457f71-8d78-442d-9182-9566e8e55cec
begin
	using CSV
	using DataFrames
	using PlutoUI
	using Shapefile
	using ZipFile
	using LsqFit
	using Plots
	using Statistics
	using Dates
	using Interpolations
	using QuadGK
	using Test
end

# ╔═╡ 3501b7fc-f721-47e6-836c-058adc261d07
using DrWatson

# ╔═╡ 9294d422-8b57-4081-8895-00851fae187b
md"# Test and documentation of functions in module `fbi.jl`"

# ╔═╡ 992dfd75-ea4b-4dc0-99ee-7c6e27ab1aee
@quickactivate "LabFBI"

# ╔═╡ 269e2e51-848a-4638-ad11-191b8abf34cd
projectdir()

# ╔═╡ 9938294d-f466-4a74-b816-8d8607da7a64
datadir()

# ╔═╡ 6d6f0cdc-a84a-41dc-8b86-ebd4c6440bcb
srcdir()

# ╔═╡ def7471f-ab92-4958-9d47-f8ef32ec8d38
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

# ╔═╡ 29e6e72d-2921-4364-938d-c0320085fcab
fbi = ingredients(srcdir("fbi.jl"))

# ╔═╡ d58cabdb-703a-4a40-b8fd-80445f01a736
plt = ingredients(srcdir("plotDF.jl"))

# ╔═╡ 0a19aa18-ea96-4836-b1b2-5fe3f0a62251
md"## Documentation"

# ╔═╡ dd1f586d-b726-4753-9f66-84a56a47c1f5
md"""
	struct Φ
		Φu::Float64
		Φc::Float64
		ϵu::Float64
		ϵc::Float64
		Φuc::Float64
		fuc::Float64
	end
Characterizes the separation between free (u) and chelated (c) molecules.


$\Phi_u = \int_{\lambda_{min}}{\lambda_{max}} f_u(\lambda) d\lambda$

$\Phi_c = \int_{\lambda_{min}}{\lambda_{max}} f_c(\lambda) d\lambda$

$\epsilon_u = \int_{\lambda_{min}}{\lambda_{max}} pdf_u(\lambda) d\lambda$

$\epsilon_c = \int_{\lambda_{min}}{\lambda_{max}} pdf_c(\lambda) d\lambda$

$\Phi_{uc} = \Phi_c / \Phi_u$

$f_{uc} = \Phi_c / \sqrt{Phi_u}$


"""

# ╔═╡ 894b13cb-0237-4f28-a8fe-905ec7c8db72
md"""
	stn_band(fs::Function, fn::Function, λmin::Float64, λmax::Float64)

Compute the signal to noise ratio between fs (signal) and fn (noise) in the band
defined by (λmin, λmax).
where:

$stn = S /\sqrt{N}$

$S = \int_{\lambda_{min}}^{\lambda_{max}} fs(\lambda) d\lambda$

$N = \int_{\lambda_{min}}^{\lambda_{max}} fn(\lambda) d\lambda$

**Arguments**

`fs::Function`   : input signal function 
`fn::Function`   : input noise function

`λmin::Float64`  : minimum of range

`λmax::Float64`  : maximum of range
"""

# ╔═╡ e7ffc713-65ee-4e10-8c33-4d9b5efbe208
md"
	double_band_ratio(fbp::Function, flp::Function,
			     	  λbpmin::Float64, λbpmax::Float64,
				      λlpmin::Float64, λlpmax::Float64)

Computes the ratio:

$r = BP/LP$

$BP = \int_{\lambda_{bpmin}}^{\lambda_{bpmax}} fbp(\lambda) d\lambda$

$LP = \int_{\lambda_{lpmin}}^{\lambda_{lpmax}} flp(\lambda) d\lambda$

**Arguments**

`fbp::Function`   : input fbp function (e.g, a function that peaks in the interval BP) 

`fl`::Function`   : input flp function (e.g, a function that peaks in the interval LP) 

`λbpmin::Float64`  : minimum of range BP

`λbpmax::Float64`  : maximum of range BP

`λlpmin::Float64`  : minimum of range LP

`λlpmax::Float64`  : maximum of range LP
"

# ╔═╡ 6b038754-90ad-4e46-82a5-72bde9d91f27
md"""
	fom(fbi::Gf, fbiba::Gf, λmin::Float64, λmax::Float64)

Return figure of merit Φ (see definition of type)

**Arguments**

`fbi::Gf`   : Generalized function defining the fbi spectrum

`fbiba::Gf` : Generalized function defining the fbiba spectrum

`λmin::Float64`  : minimum of range 

`λmax::Float64`  : maximum of range 
"""

# ╔═╡ 6a637888-2cec-4a88-a237-e9f3c123af92
md"## Tests"

# ╔═╡ 5527eeff-471b-4cfb-a499-12e4e24bba28
md"
- Test eps band definition using interpolated and generalized functions
- Test double band ratio definition using interpolated and generalized functions
"

# ╔═╡ 9a048427-a11f-4e20-bfce-4c4549998e4f
begin
    wf = 350.
    we = 800.
    ws = 2.
    wl = wf:ws:we
end


# ╔═╡ 4334b44e-71e5-4578-8488-f4cea53c8dbf
begin
	λminbp = 425.0  # in nm
	λmaxbp = 435.0  # in nm
	λminlp = 450.0  # in nm
	λmaxlp = 800.0  # in nm
end

# ╔═╡ ae514258-e949-46dc-9aee-7c9f5f71b875
fbidf = fbi.load_df_from_csv(datadir("fluorimeter/325nm"),
	                                  "FBI_G1_vs_FBI_G2_em325nm.csv",  fbi.spG)


# ╔═╡ 51969781-4390-4478-9aba-b78dd50e7d92
begin
	fbig1 = fbi.dftof(wl, fbidf, "FBIG1")
	fbibag1 = fbi.dftof(wl, fbidf, "FBIBaG1")
	fbig2 = fbi.dftof(wl, fbidf, "FBIG2")
	fbibag2 = fbi.dftof(wl, fbidf, "FBIBaG2")
end

# ╔═╡ 4063ab15-3818-43e9-98f4-df1e24598df6
begin
	fbigfg1   = fbi.dftogf(wl, fbidf, "FBIG1")
	fbibagfg1 = fbi.dftogf(wl, fbidf, "FBIBaG1")
	fbigfg2   = fbi.dftogf(wl, fbidf, "FBIG2")
	fbibagfg2 = fbi.dftogf(wl, fbidf, "FBIBaG2")
end

# ╔═╡ 8a3642e1-f9bb-4ac7-ac45-f6b8ed4fc3db
@test fbi.stn_band(fbibagfg1.f, fbigfg1.f, λminbp, λmaxbp) ≈ fbi.qpdf(fbibag1, λminbp, λmaxbp) / sqrt(fbi.qpdf(fbig1, λminbp, λmaxbp))


# ╔═╡ dac05e87-c27c-48fd-9334-c08d32e60c53
@test fbi.stn_band(fbibagfg2.f, fbigfg2.f, λminbp, λmaxbp) ≈ fbi.qpdf(fbibag2, λminbp, λmaxbp) / sqrt(fbi.qpdf(fbig2, λminbp, λmaxbp))


# ╔═╡ fa6341f1-7baf-4772-8052-af2b8c75af2a
@test fbi.double_band_ratio(fbibagfg1.f, fbigfg1.f,
	                       λminbp, λmaxbp,
				           λminbp, λmaxbp) ≈ fbi.qpdf(fbibag1, λminbp, λmaxbp) / fbi.qpdf(fbig1, λminbp, λmaxbp)


# ╔═╡ fd353197-ac4f-4713-81b3-dac0c4f395d1
Phi = fbi.fom(fbigfg1, fbibagfg1, λminbp, λmaxbp)

# ╔═╡ 484be2b8-c10e-4a2c-a8a9-c2a1c7dae498
@test Phi.Φu ≈fbi.qpdf(fbig1, λminbp, λmaxbp)

# ╔═╡ 5d49be14-df4e-4d95-88cf-82ed0c07bb0d
@test Phi.ϵu ≈fbi.qpdf(fbig1, λminbp, λmaxbp) / fbi.qpdf(fbig1, 0.0, 800.0)

# ╔═╡ 02af59d6-7ad7-49d6-bdc0-d46fba086750
@test Phi.Φc ≈fbi.qpdf(fbibag1, λminbp, λmaxbp)

# ╔═╡ 888f8987-5ccc-48dc-9424-b7974db2f421
@test Phi.ϵc ≈fbi.qpdf(fbibag1, λminbp, λmaxbp) / fbi.qpdf(fbibag1, 0.0, 800.0)

# ╔═╡ 1b54c464-21eb-460f-8ae4-90f85ff4a2ff
@test Phi.Φuc ≈fbi.double_band_ratio(fbibagfg1.f, fbigfg1.f,
	                       λminbp, λmaxbp,
				           λminbp, λmaxbp)


# ╔═╡ 4d29e1b7-9bbb-4a26-9a32-cb954a47438b
@test Phi.fuc ≈fbi.stn_band(fbibagfg1.f, fbigfg1.f, λminbp, λmaxbp)

# ╔═╡ Cell order:
# ╠═9e457f71-8d78-442d-9182-9566e8e55cec
# ╠═3501b7fc-f721-47e6-836c-058adc261d07
# ╠═9294d422-8b57-4081-8895-00851fae187b
# ╠═992dfd75-ea4b-4dc0-99ee-7c6e27ab1aee
# ╠═269e2e51-848a-4638-ad11-191b8abf34cd
# ╠═9938294d-f466-4a74-b816-8d8607da7a64
# ╠═6d6f0cdc-a84a-41dc-8b86-ebd4c6440bcb
# ╟─def7471f-ab92-4958-9d47-f8ef32ec8d38
# ╠═29e6e72d-2921-4364-938d-c0320085fcab
# ╠═d58cabdb-703a-4a40-b8fd-80445f01a736
# ╟─0a19aa18-ea96-4836-b1b2-5fe3f0a62251
# ╟─dd1f586d-b726-4753-9f66-84a56a47c1f5
# ╟─894b13cb-0237-4f28-a8fe-905ec7c8db72
# ╟─e7ffc713-65ee-4e10-8c33-4d9b5efbe208
# ╟─6b038754-90ad-4e46-82a5-72bde9d91f27
# ╠═6a637888-2cec-4a88-a237-e9f3c123af92
# ╠═5527eeff-471b-4cfb-a499-12e4e24bba28
# ╠═9a048427-a11f-4e20-bfce-4c4549998e4f
# ╠═4334b44e-71e5-4578-8488-f4cea53c8dbf
# ╠═ae514258-e949-46dc-9aee-7c9f5f71b875
# ╠═51969781-4390-4478-9aba-b78dd50e7d92
# ╠═4063ab15-3818-43e9-98f4-df1e24598df6
# ╠═8a3642e1-f9bb-4ac7-ac45-f6b8ed4fc3db
# ╠═dac05e87-c27c-48fd-9334-c08d32e60c53
# ╠═fa6341f1-7baf-4772-8052-af2b8c75af2a
# ╠═fd353197-ac4f-4713-81b3-dac0c4f395d1
# ╠═484be2b8-c10e-4a2c-a8a9-c2a1c7dae498
# ╠═5d49be14-df4e-4d95-88cf-82ed0c07bb0d
# ╠═02af59d6-7ad7-49d6-bdc0-d46fba086750
# ╠═888f8987-5ccc-48dc-9424-b7974db2f421
# ╠═1b54c464-21eb-460f-8ae4-90f85ff4a2ff
# ╠═4d29e1b7-9bbb-4a26-9a32-cb954a47438b

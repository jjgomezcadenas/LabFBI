### A Pluto.jl notebook ###
# v0.14.3

using Markdown
using InteractiveUtils

# ╔═╡ 5d5e5072-7ad8-11eb-1559-9763c2eff6a0
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

# ╔═╡ 94f44f29-d95a-40bb-b3b4-ba217ff3cb21
using DrWatson

# ╔═╡ f58fd9ea-7b52-11eb-2ad5-c5e14f503d9d
md"# Test and documentation of functions in module `dfftest.jl`"

# ╔═╡ d4f547ff-efcf-4c70-ab38-9ed2abf52296
#begin
#	using Pkg
#	Pkg.add.(["CSV", "DataFrames", "PlutoUI", "Shapefile", "ZipFile", "LsqFit", 	"Plots","Interpolations","QuadGK","Test"])
#end

# ╔═╡ 914b98d3-eca5-411b-9ddb-032540e51ab7
@quickactivate "LabFBI"

# ╔═╡ d375839e-69cc-4ce4-abaf-23d71c837950
projectdir()

# ╔═╡ cbe4f338-3c9a-49a9-8efc-0f31f9ed6c29
datadir()

# ╔═╡ ffcb23c5-26f2-4f27-939d-7d855af1b142
srcdir()

# ╔═╡ e4642e06-4fac-4226-a7d6-030cb5709f95
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

# ╔═╡ f6d37f9c-5647-4b65-af0d-c62bd9e37dc9
fbi = ingredients(srcdir("fbi.jl"))

# ╔═╡ 0f3dc8ca-6917-459e-8c77-cc9e23286d99
plt = ingredients(srcdir("plotDF.jl"))

# ╔═╡ baebd1fc-1f08-4c44-b63f-aad0638a5738
md"## `dffunctions.jl`"

# ╔═╡ a767cf02-d743-49b8-a266-019fb7c3fc48
md"""
	struct CsvG  # stands for CSV Grammar
		delim::Char
		decimal::Char
	end
Describe a simple grammar needed to read CSV files

**Fields**
- `delim::Char`  : A character defining the delimiter (eg. ',' ';' etc)
- `decimal::Char`: A character defining the decimal (eg. ',' '.')
- `pdf::Function`: Normalized to area 1 (PDF).
	
"""

# ╔═╡ 0ebe3fa5-46e5-41fc-bc00-3b6e4d9d0b5d
md"""
		dftof(wl::AbstractRange, df::DataFrame, cname::String,
		  bkgnd::Float64=0.0)

Return an interpolated function, valid in range wl.

**Arguments**
- `wl::Range`: interpolation range.
- `df::DataFrame`: data frame holding the data.
- `cname::String`: name of the column holding data to be interplated.
- `bkgnd::Float64`: value of the data outside interpolation range
                    (zero by default).
"""

# ╔═╡ a7c7832e-0b2d-41d3-9306-e4d900be4fc5
md"""
	qpdf(f::Function, λmin::Number, λmax::Number)

Return the integral of f in the interval (λmin, λmax)
(Syntactic sugar for quadgk)

**Arguments**
- `f::Function`: Function to be integrated.
- `λmin::Number`: lower bound of range.
- `λmax::Number`: upper bound of range.
"""

# ╔═╡ 622edbc0-bae3-421c-8324-17dd82ecb178
md"""
	ftopdf(wl::AbstractRange, f::Function)

Compute the PDF of function f in range wl.

**Arguments**
- `wl::AbstractRange`: Range of application.
- `f::Function`: Input function.
"""

# ╔═╡ 4607ee65-f8c3-4f86-a89c-436ab40560d8
md"""
	struct Gf
		N  ::Number
		f  ::Function
		pdf::Function
	end

Represents a generalized function:

**Fields**
- `N::Number`: Normalization constant: pdf(x) = f(x)/N.
- `f::Function`: Any integrable function.
- `pdf::Function`: Normalized to area 1 (PDF).
	
"""

# ╔═╡ f2731e25-581c-454b-83a9-15d08d3c45f3
md"""
	dftogf(wl::AbstractRange, df::DataFrame, cname::String, 
	       bkgnd::Float64=0.0)
Return a generalized function from dataframe, column and range

**Arguments**
- `wl::AbstractRange`: Range of application.
- `df::DataFrame`: data frame holding the data.
- `f::Function`: Input function.
- `cname::String`: name of the column holding data to be interplated.
- `bkgnd::Float64`: value of the data outside interpolation range
                    (zero by default).
"""

# ╔═╡ bba127ba-2f4b-44d3-a0e9-3d8f54edb84c
md"""
	find_max_xy(df, xc, yc)

Return ymax and x such that ymax = f(x).

Description:
In a DataFrame one has often "XY" variables, that is, a pair of columns
"X" and "Y" which represent correlated variables (e.g intensity and wavelength).
In such cases one often wants the XY maximum, that is, finding the maximum
of "Y" (call it ymax) and the corresponding value in "X"
(x for which y is ymax). This corresponds, for example, to the wavelength at
which the intensity is maximal.

**Arguments**
- `df::DataFrame`: data frame holding the data.
- `xc::String`: the X column.
- `yc::String`: the Y column.

"""

# ╔═╡ 8f50282d-db25-4658-bddf-ba8d2d4eaae7
md"## Tests `dffunctions.jl`
- Test that the PDF (normalize function) is equal to the interpolated function divided by the normalization constant
- Test that the interpolation function equals the DF in all nodes
"

# ╔═╡ f50f9466-e52d-4ec8-9c68-cfc8cd904b4e
fbidf = fbi.load_df_from_csv(datadir("fluorimeter/325nm"),
	                                  "FBI_G1_vs_FBI_G2_em325nm.csv",  fbi.spG) 

# ╔═╡ fee8f210-d8a0-4128-accb-846a5fdcc419
begin
    wf = 350.
    we = 800.
    ws = 2.
    wl = wf:ws:we
end

# ╔═╡ e71bd44f-b3c3-4ea0-9c11-a3a08afe153c
begin
	fbig1 = fbi.dftof(wl, fbidf, "FBIG1")
	fbibag1 = fbi.dftof(wl, fbidf, "FBIBaG1")
	fbig2 = fbi.dftof(wl, fbidf, "FBIG2")
	fbibag2 = fbi.dftof(wl, fbidf, "FBIBaG2")
end

# ╔═╡ de72040c-64ee-46fc-9f9c-9f1eafd5dd38
begin
	fbigfg1   = fbi.dftogf(wl, fbidf, "FBIG1")
	fbibagfg1 = fbi.dftogf(wl, fbidf, "FBIBaG1")
	fbigfg2   = fbi.dftogf(wl, fbidf, "FBIG2")
	fbibagfg2 = fbi.dftogf(wl, fbidf, "FBIBaG2")
end

# ╔═╡ 9c4311ef-fc12-4b2d-b8e7-c60190ab7e00
function test1_df_functions(w, fbis)
	L = collect(wl)
	test1 = [(F.f).(L) / (F.N) ≈ (F.pdf).(L) for F in fbis]
	return all(test1)
end

# ╔═╡ 844c7591-dbd2-4ffa-a425-c728071c18b0
function test2_df_functions(w, fbis)
	L = collect(wl)
	LBL =["FBIG1", "FBIBaG1", "FBIG2", "FBIBaG2"]
	test2 =[(F.f).(L) ≈ fbidf[!,LBL[i]] for (i,F) in enumerate(fbis)]
	return all(test2)
end

# ╔═╡ cb0f28ed-b34d-4bf5-a055-ad3c8dbea5c1
function test_fbi_pdfs(wr, fbis)
	Ll = collect(wr)
	test = [fbi.qpdf(f.pdf, wr[1], wr[end]) ≈ 1 for f in fbis]
	return all(test)
end


# ╔═╡ 965e86ea-bbe4-4e46-af0b-c972c0a533cf
function test_maximum_fbi(wl, fbis)
	L = collect(wl)
	LBL =["FBIG1", "FBIBaG1", "FBIG2", "FBIBaG2"]
	test =[maximum(F.(L)) ≈ fbi.find_max_xy(fbidf, "λ", LBL[i])[1] for (i,F) in enumerate(fbis)]
	return all(test)
end

# ╔═╡ 0e2a18d4-9c2f-495e-8d43-34f824b71b91
@test test_maximum_fbi(wl, [fbig1,fbibag1,fbig2,fbibag2]) 

# ╔═╡ ca5eb8e7-b246-4737-8ce4-16b86da3a48d
@test test1_df_functions(wl, [fbigfg1,fbibagfg1,fbigfg2,fbibagfg2]) 

# ╔═╡ 13f89e00-0a32-47ed-84e6-dc0b8074ffbe
@test test2_df_functions(wl, [fbigfg1,fbibagfg1,fbigfg2,fbibagfg2])

# ╔═╡ 24ec362f-9038-4b8c-85aa-4845d3fbb936
@test test_fbi_pdfs(wl, [fbigfg1,fbibagfg1,fbigfg2,fbibagfg2])

# ╔═╡ 70ded0d4-0052-4fad-920a-6dc985ae41ee
md"### Plot functions"

# ╔═╡ 3bb51375-8f81-4f64-96c2-147f6c88b23e
md"""
	plotdf_xy(df::DataFrame, dfx::String, dfy::String,
		 	  xlabel, ylabel;
		 	  label="DF data", shape = :circle, color = :black,
		 	  markersize = 3, legend=false)

Plot columns (dx,dy) of dataframe df

"""

# ╔═╡ 12b7c030-f194-49e8-9801-00d92d35655a
md"""
	plotf(f, X, xlabel, ylabel, title;
		  label="DF", lw = 2, color = :black, legend=false)

Plot f.(X) where X is a vector
"""

# ╔═╡ 5fb83f9d-c4d4-494b-9339-2c4763a71926
md"""
	merge_plots!(plt, plts...)

Merge plts plts with plt.
"""

# ╔═╡ 6a167929-257c-4036-b6a8-726ae2cc9c48
pfbidf = plt.plotdf_xy(fbidf, "λ", "FBIG1", "λ (nm)", "I (au)", color=:green, label="Fbi-325nm");

# ╔═╡ 3aa65b41-a6bc-4f10-aea6-3b30571c4f12
pfbibadf = plt.plotdf_xy(fbidf, "λ", "FBIBaG1",  "λ (nm)", "I (au)", color=:Blue, label="FbiBa-325nm");

# ╔═╡ 164d1f4c-0f3f-4872-8de0-6eafa2307ccc
pfbif = plt.plotf(fbig1, collect(wl), "λ (nm)", "I (au)", "FBI interpol", color=:black, label="Fbi-325nm-i", legend=true);

# ╔═╡ f98e679b-1f4e-460f-a6be-1dacde7cc09e
pfbi= plt.merge_plots!(pfbidf, pfbif);

# ╔═╡ 601cb98a-e82f-4361-a6f3-ea0f37b96b22
pfbibaf = plt.plotf(fbibag1, collect(wl), "λ (nm)", "I (au)", "FBIBa interpol", color=:black, label="FbiBa-325nm-i", legend=true);

# ╔═╡ e289b9e2-9b30-4068-ad44-60665a5f5146
pfbiba= plt.merge_plots!(pfbibadf, pfbibaf);

# ╔═╡ 62fa1442-945b-4899-8402-67a1ff6d321e
plot(pfbi, pfbiba, layout = (1, 2), legend = true, fmt = :png)

# ╔═╡ 2137d817-d5cd-4b06-bb92-6772571d843c
pF= plt.merge_plots!(pfbi, pfbiba)

# ╔═╡ 28db0899-2bd2-4e6d-b7fd-29f5670e8950
pfbigfg1 = plt.plotf(fbigfg1.pdf, collect(wl), "λ (nm)", "I (au)", "FBI PDF", color=:green, label="Fbi G1", legend=true);

# ╔═╡ 9e0fd9c7-5c9a-49d7-8f86-10e0e6753f02
pfbibagfg1 = plt.plotf(fbibagfg1.pdf, collect(wl), "λ (nm)", "I (au)", "FBIBa G1 PDF", color=:blue, label="FbiBa G1", legend=true);

# ╔═╡ 0b75077f-6141-4b51-ae24-5fd9d96c39a9
pG1= plt.merge_plots!(pfbigfg1, pfbibagfg1);

# ╔═╡ 8de65552-7fc8-4b51-b721-82a854af4a9e
pfbigfg2 = plt.plotf(fbigfg2.pdf, collect(wl), "λ (nm)", "I (au)", "FBI G2 PDF", color=:red, label="Fbi G2", legend=true);

# ╔═╡ 97861dcb-e26c-44fd-8234-6a129f950dab
pfbibagfg2 = plt.plotf(fbibagfg2.pdf, collect(wl), "λ (nm)", "I (au)", "FBIBa G2 PDF", color=:orange, label="FbiBa G2", legend=true);

# ╔═╡ de3e2e3f-17d3-4a95-8ea8-75a99c7fffec
pG2= plt.merge_plots!(pfbigfg2, pfbibagfg2)

# ╔═╡ Cell order:
# ╠═f58fd9ea-7b52-11eb-2ad5-c5e14f503d9d
# ╠═d4f547ff-efcf-4c70-ab38-9ed2abf52296
# ╠═5d5e5072-7ad8-11eb-1559-9763c2eff6a0
# ╠═94f44f29-d95a-40bb-b3b4-ba217ff3cb21
# ╠═914b98d3-eca5-411b-9ddb-032540e51ab7
# ╠═d375839e-69cc-4ce4-abaf-23d71c837950
# ╠═cbe4f338-3c9a-49a9-8efc-0f31f9ed6c29
# ╠═ffcb23c5-26f2-4f27-939d-7d855af1b142
# ╠═e4642e06-4fac-4226-a7d6-030cb5709f95
# ╠═f6d37f9c-5647-4b65-af0d-c62bd9e37dc9
# ╠═0f3dc8ca-6917-459e-8c77-cc9e23286d99
# ╠═baebd1fc-1f08-4c44-b63f-aad0638a5738
# ╟─a767cf02-d743-49b8-a266-019fb7c3fc48
# ╟─0ebe3fa5-46e5-41fc-bc00-3b6e4d9d0b5d
# ╟─a7c7832e-0b2d-41d3-9306-e4d900be4fc5
# ╟─622edbc0-bae3-421c-8324-17dd82ecb178
# ╟─4607ee65-f8c3-4f86-a89c-436ab40560d8
# ╟─f2731e25-581c-454b-83a9-15d08d3c45f3
# ╟─bba127ba-2f4b-44d3-a0e9-3d8f54edb84c
# ╟─8f50282d-db25-4658-bddf-ba8d2d4eaae7
# ╠═f50f9466-e52d-4ec8-9c68-cfc8cd904b4e
# ╠═fee8f210-d8a0-4128-accb-846a5fdcc419
# ╠═e71bd44f-b3c3-4ea0-9c11-a3a08afe153c
# ╠═de72040c-64ee-46fc-9f9c-9f1eafd5dd38
# ╠═9c4311ef-fc12-4b2d-b8e7-c60190ab7e00
# ╠═844c7591-dbd2-4ffa-a425-c728071c18b0
# ╠═cb0f28ed-b34d-4bf5-a055-ad3c8dbea5c1
# ╠═965e86ea-bbe4-4e46-af0b-c972c0a533cf
# ╠═0e2a18d4-9c2f-495e-8d43-34f824b71b91
# ╠═ca5eb8e7-b246-4737-8ce4-16b86da3a48d
# ╠═13f89e00-0a32-47ed-84e6-dc0b8074ffbe
# ╠═24ec362f-9038-4b8c-85aa-4845d3fbb936
# ╠═70ded0d4-0052-4fad-920a-6dc985ae41ee
# ╟─3bb51375-8f81-4f64-96c2-147f6c88b23e
# ╟─12b7c030-f194-49e8-9801-00d92d35655a
# ╟─5fb83f9d-c4d4-494b-9339-2c4763a71926
# ╠═6a167929-257c-4036-b6a8-726ae2cc9c48
# ╠═3aa65b41-a6bc-4f10-aea6-3b30571c4f12
# ╠═164d1f4c-0f3f-4872-8de0-6eafa2307ccc
# ╠═f98e679b-1f4e-460f-a6be-1dacde7cc09e
# ╠═601cb98a-e82f-4361-a6f3-ea0f37b96b22
# ╠═e289b9e2-9b30-4068-ad44-60665a5f5146
# ╠═62fa1442-945b-4899-8402-67a1ff6d321e
# ╠═2137d817-d5cd-4b06-bb92-6772571d843c
# ╠═28db0899-2bd2-4e6d-b7fd-29f5670e8950
# ╠═9e0fd9c7-5c9a-49d7-8f86-10e0e6753f02
# ╠═0b75077f-6141-4b51-ae24-5fd9d96c39a9
# ╠═8de65552-7fc8-4b51-b721-82a854af4a9e
# ╠═97861dcb-e26c-44fd-8234-6a129f950dab
# ╠═de3e2e3f-17d3-4a95-8ea8-75a99c7fffec

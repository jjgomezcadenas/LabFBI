### A Pluto.jl notebook ###
# v0.14.4

using Markdown
using InteractiveUtils

# ╔═╡ 167d8fa6-0297-49ab-83f7-49113da458a3
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
	using DrWatson
	using GLM
end

# ╔═╡ e790f112-a73c-11eb-0997-47b38885410b
md"# FBI Fluorescence"

# ╔═╡ e3e38f9c-668a-453a-9ce1-e4b8a98ea931
@quickactivate "LabFBI"

# ╔═╡ 5c75d2f1-75a1-43d1-b604-9c3b0d116c67
projectdir()

# ╔═╡ 3860a099-3acb-4f2c-9926-8f0221500435
datadir()

# ╔═╡ ac6e7ea2-0d99-41ce-a600-b5cdffb21b8d
srcdir()

# ╔═╡ f5b87a7b-8345-4c05-a99f-48f87ae02f6a
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

# ╔═╡ ac21609b-e68a-418f-a18a-857492d6533c
dff = ingredients(srcdir("fbi.jl"))

# ╔═╡ 69dff321-a39c-4322-ba3f-4372aef840a5
dfp = ingredients(srcdir("plotDF.jl"))

# ╔═╡ 8eb858ce-83f1-490a-b85e-6edb6a4eb193
utils = ingredients(srcdir("utils.jl"))

# ╔═╡ 0540fcda-f003-428e-a2da-46764339059f
math = ingredients(srcdir("math.jl"))

# ╔═╡ 405ad1d1-aaf0-4dae-9f2a-fac0e60fad6d
colors = [:green :orange :black :purple :red  :yellow :brown :white]

# ╔═╡ f4a15e73-aeb1-4c90-a96c-20745d4d84ff
md"## Notebook"

# ╔═╡ 0f7e6f41-d418-4388-a947-9a1110bd7148
md"### Fluorescence 325 & 405 nm: FBI G1"

# ╔═╡ d09fb693-58e4-40f3-b3f2-0bda0550dffe
md"""
#### Main results:

- Fluorescence is not lineal with concentration for concentrations above ``10^{-6}``M. 
- No sensitivity below ``5 \times 10^{-8}``M.
- Free molecule shows a shoulder in blue region (circa 425 nm) for excitation at 325 nm not present at 405 nm
"""

# ╔═╡ 30c959a8-c36a-4088-95c6-8cadd5e28e17
g1em325df = dff.load_df_from_csv(datadir("fluorimeter/fluorescence"),
	                                  "FbiG1Em325.csv",  dff.spG) ;


# ╔═╡ 368e8ac9-ffd2-442e-a38c-a61f6e1c4356
g1em405df = dff.load_df_from_csv(datadir("fluorimeter/fluorescence"),
	                                  "FbiG1Em405.csv",  dff.spG) ;

# ╔═╡ 8586f28f-93e6-4508-8818-41a7004d4f3e
begin
	columns = names(g1em325df)[2:end]
	cs = [parse.(Float64, split(column, "_")[2]) for column in columns]
	col405 = names(g1em405df)[2:end]
	cs405 = [parse.(Float64, split(column, "_")[2]) for column in col405]
end

# ╔═╡ 0a5d227a-f2d7-45ae-a836-3eacfc313a98
pf325 = dfp.plotdf_xys(g1em325df, "λ", columns,
	                columns, colors,
	                "λ (nm)", "I (AU)", "Fbi G1; Em: 325 nm",
		            legend=:topright);

# ╔═╡ c8325077-fe1e-439b-b16c-6b795eb92245
pf405 = dfp.plotdf_xys(g1em405df, "λ", col405,
	                col405, colors,
	                "λ (nm)", "I (AU)", "Fbi G1; Em: 405 nm",
		            legend=:topright);

# ╔═╡ 88214d2e-2706-4ab4-9f11-1bba0bb421e1
plot(pf325, pf405, layout = (1, 2), legend=:topright, fmt = :png)


# ╔═╡ a180ec0f-c47d-47cb-82eb-b120a32a7ed6
md"### Fluorescence 325 & 405 nm: FbiBa G1"

# ╔═╡ 6b589891-16ff-42b8-95b0-a66e2f94c62b
g1baem325df = dff.load_df_from_csv(datadir("fluorimeter/fluorescence"),
	                                  "FbiBaG1Em325.csv",  dff.spG) ;

# ╔═╡ 44469ce2-173d-48c6-bb00-fea6211b4b21
g1baem405df = dff.load_df_from_csv(datadir("fluorimeter/fluorescence"),
	                                  "FbiBaG1Em405.csv",  dff.spG) ;

# ╔═╡ e23bd3ad-f76d-43c2-8a91-3b36a65db389
begin
	colba = names(g1baem325df)[2:end]
	csba = [parse.(Float64, split(column, "_")[3]) for column in colba]
	colba405 = names(g1baem405df)[2:end]
end

# ╔═╡ d21a4e66-6621-4c2d-97ff-85f8a6a09c34
pfba325 = dfp.plotdf_xys(g1baem325df, "λ", colba,
	                colba, colors,
	                "λ (nm)", "I (AU)", "Fbi G1; Em: 325 nm",
		            legend=:topright);

# ╔═╡ 66278a87-23c2-4354-9239-2c9a577ea68f
pfba405 = dfp.plotdf_xys(g1baem405df, "λ", colba405,
	                colba405, colors,
	                "λ (nm)", "I (AU)", "Fbi G1; Em: 405 nm",
		            legend=:topright);

# ╔═╡ 298e45c4-dc70-4088-a523-98d92f496e9b
plot(pfba325, pfba405, layout = (1, 2), legend=:topright, fmt = :png)


# ╔═╡ 6bd6d4e9-0e13-420a-ac5c-cda420c392a4
md"### Fluorescence 325: FbiBa G2

#### Main results:

- Fluorescence is not lineal with concentration for concentrations above ``10^{-6}``M. 
- No sensitivity below ``5 \times 10^{-8}``M.
- Peak of free molecule around 620 nm
"

# ╔═╡ 0358b5d4-4574-444e-95de-dcafdfc84bd3
g2em325df = dff.load_df_from_csv(datadir("fluorimeter/fluorescence"),
	                                  "FbiG2Em325.csv",  dff.spG) 

# ╔═╡ c431ca45-3c51-410f-aeff-20d8e5f3c725
begin
	g2col325 = names(g2em325df)[2:end]
	g2cs = [parse.(Float64, split(column, "_")[2]) for column in columns]
end

# ╔═╡ a144b6bf-9986-4246-9a24-aeaa69b82b80
g2p325 = dfp.plotdf_xys(g2em325df, "λ", g2col325,
	                g2col325, colors,
	                "λ (nm)", "I (AU)", "Fbi G2; Em: 325 nm",
		            legend=:topright)

# ╔═╡ 341c8e71-b73d-451d-bba9-aae64269d01d
md"### Shape invariance

- Unlike G1, G2 does not change shape with concentration 
"

# ╔═╡ 8ff15860-93a5-4669-9349-f6dd10f97cb7
wl = 350.0:2.0:900.0

# ╔═╡ e616c755-4d69-422f-b694-9c76dcaf0334
gfs = dff.dftogf.((wl,), (g2em325df,), g2col325)

# ╔═╡ 7b959d38-ae82-4355-87ee-5c2f4ff937c2
pgdfs = dfp.plotdf_gfs(gfs, wl,1:5,
	                g2col325, colors,
	                "λ (nm)", "I (pdf)", "Fbi G2; Em: 325 nm",
		            legend=:topright)

# ╔═╡ b9a30be5-17d7-4137-af7c-4935c567e0f8
md"### FbiBa G2 "

# ╔═╡ 3373fef3-08c9-4bb7-9d33-814ce4eedf24
g2ba325df = dff.load_df_from_csv(datadir("fluorimeter/fluorescence"),
	                                  "FbiBaG2Em325.csv",  dff.spG) 

# ╔═╡ 81a28886-603a-4618-a98c-d6562fc8f5a1
g2bacol325 = names(g2ba325df)[2:end]

# ╔═╡ 53e4fc32-c61a-4b17-90de-45fe0bb1c233
g2bap325 = dfp.plotdf_xys(g2ba325df, "λ", g2bacol325,
	                g2bacol325, colors,
	                "λ (nm)", "I (AU)", "FbiBa G2; Em: 325 nm",
		            legend=:topright)

# ╔═╡ b345f6ed-9224-418c-8076-d69e7e48f01e
gfsba = dff.dftogf.((wl,), (g2ba325df,), g2bacol325)

# ╔═╡ da49a021-4bfc-43c8-91f5-a93f9ffc4270
pbagdfs = dfp.plotdf_gfs(gfsba, wl,1:5,
	                g2bacol325, colors,
	                "λ (nm)", "I (pdf)", "FbiBa G2; Em: 325 nm",
		            legend=:topright)

# ╔═╡ dd00edef-1382-4f05-a632-18bcbd6d3872
plot(pgdfs, pbagdfs, layout = (1, 2), legend=false, fmt = :png)


# ╔═╡ f2d4aa10-fda1-4fdb-b134-6526b225c726
dfp.merge_plots!(pgdfs, pbagdfs)

# ╔═╡ 07aafc3c-a0c5-4f50-8f29-8e686a28b828
md"### Fbi G2 metales"

# ╔═╡ ae1e97ca-bcc0-4ba9-b878-1f38eafc714f
g2mtem325df = dff.load_df_from_csv(datadir("fluorimeter/fluorescence"),
	                                  "FbiMetalesG2Em325.csv",  dff.spG) 

# ╔═╡ 0488a732-a7c0-4531-aeab-d7a0b41701b7
g2mtcol325 = names(g2mtem325df)[2:end]

# ╔═╡ 02776897-4a63-4066-aff9-e308151f715a
g2mtp325 = dfp.plotdf_xys(g2mtem325df, "λ", g2mtcol325,
	                g2mtcol325, colors,
	                "λ (nm)", "I (AU)", "Fbi G2 Metals; Em: 325 nm",
		            legend=:topright)

# ╔═╡ 4f77b21b-b16a-41f4-bcc0-459e79db8896
md"### Functions"

# ╔═╡ a47c2d44-b11f-4d5e-810a-dfa83957a683
function maxfluo(df::DataFrame, columns::Array{String})
	maxXY = dff.find_max_xy.((df,), ("λ",), columns)
	fm = collect(zip(maxXY...))
	return collect(fm[1]), collect(fm[2])
end

# ╔═╡ 5b3d0a71-00ff-44fe-b004-9d15c92a1c2a
fmaxg2, lmaxg2 = maxfluo(g2em325df, g2col325)

# ╔═╡ c2f95a5d-3711-45af-835f-3770d1612954
g2peak = dfp.plot_xy(cs, fmaxg2, "C (M)", "I (AU)", "I (peak) vs C (G2 325 nm)")	

# ╔═╡ 01e7991c-91ce-48b9-969e-5a71bf055f28
function afluo(df::DataFrame, columns::Array{String})
	return [sum(df[!,col]) for col in columns]
end

# ╔═╡ 5c6c34df-c47b-474d-8762-b1b03da9451c
begin
	fmax, lmax = maxfluo(g1em325df, columns);
	Is = afluo(g1em325df, columns);
end

# ╔═╡ 16152d4e-9d2b-417d-ae10-4967e472250c
plnpeak = dfp.plot_xy(cs, fmax, "C (M)", "I (AU)", "I (peak) vs C (325 nm)");	

# ╔═╡ 71244d78-a021-48b1-9625-ad8fb4a403cb
plna =dfp.plot_xy(cs, Is, "C (M)", "I (AU)", "I (sum) vs C (325 nm)");

# ╔═╡ 73889cb4-49b3-4e24-8b74-441cf8edbafc
plot(plnpeak, plna, layout = (1, 2), legend=false, fmt = :png)


# ╔═╡ afafda77-71f5-447a-a1ff-d3cd5a733d15
plnpeaka = dfp.plot_xy(cs[4:end-1], fmax[4:end-1], 
	"C (M)", "I (AU)", "I (peak) vs C")	;

# ╔═╡ 3ecd0d69-08e8-491b-9eae-ff2764395857
plnpeakb = dfp.plot_xy(cs[3:end-1], fmax[3:end-1], 
	"C (M)", "I (AU)", "I (peak) vs C")	;

# ╔═╡ fd3aab21-21b9-4dab-bff6-d792f3191bcf
plnpeakc = dfp.plot_xy(cs[2:end-1], fmax[2:end-1], 
	"C (M)", "I (AU)", "I (peak) vs C")	;

# ╔═╡ 24e78bfb-97de-4f2f-8baa-2fceb0b94a97
plot(plnpeaka, plnpeakb, plnpeakc, layout = (3, 1), legend=false, fmt = :png)

# ╔═╡ 9ab91673-619d-4d29-80d7-7e2e37898f00
Is405 = afluo(g1em405df, col405);

# ╔═╡ a1924c0b-d0a9-44e9-833b-a6216e9e0627
plna405 =dfp.plot_xy(cs405, Is405, "C (M)", "I (AU)", "I (sum) vs C (405 nm)")

# ╔═╡ Cell order:
# ╠═e790f112-a73c-11eb-0997-47b38885410b
# ╠═167d8fa6-0297-49ab-83f7-49113da458a3
# ╠═e3e38f9c-668a-453a-9ce1-e4b8a98ea931
# ╠═5c75d2f1-75a1-43d1-b604-9c3b0d116c67
# ╠═3860a099-3acb-4f2c-9926-8f0221500435
# ╠═ac6e7ea2-0d99-41ce-a600-b5cdffb21b8d
# ╠═f5b87a7b-8345-4c05-a99f-48f87ae02f6a
# ╠═ac21609b-e68a-418f-a18a-857492d6533c
# ╠═69dff321-a39c-4322-ba3f-4372aef840a5
# ╠═8eb858ce-83f1-490a-b85e-6edb6a4eb193
# ╠═0540fcda-f003-428e-a2da-46764339059f
# ╠═405ad1d1-aaf0-4dae-9f2a-fac0e60fad6d
# ╠═f4a15e73-aeb1-4c90-a96c-20745d4d84ff
# ╠═0f7e6f41-d418-4388-a947-9a1110bd7148
# ╠═d09fb693-58e4-40f3-b3f2-0bda0550dffe
# ╠═88214d2e-2706-4ab4-9f11-1bba0bb421e1
# ╠═73889cb4-49b3-4e24-8b74-441cf8edbafc
# ╠═30c959a8-c36a-4088-95c6-8cadd5e28e17
# ╠═368e8ac9-ffd2-442e-a38c-a61f6e1c4356
# ╠═8586f28f-93e6-4508-8818-41a7004d4f3e
# ╠═0a5d227a-f2d7-45ae-a836-3eacfc313a98
# ╠═c8325077-fe1e-439b-b16c-6b795eb92245
# ╠═5c6c34df-c47b-474d-8762-b1b03da9451c
# ╠═16152d4e-9d2b-417d-ae10-4967e472250c
# ╠═71244d78-a021-48b1-9625-ad8fb4a403cb
# ╠═afafda77-71f5-447a-a1ff-d3cd5a733d15
# ╠═3ecd0d69-08e8-491b-9eae-ff2764395857
# ╠═fd3aab21-21b9-4dab-bff6-d792f3191bcf
# ╠═24e78bfb-97de-4f2f-8baa-2fceb0b94a97
# ╠═9ab91673-619d-4d29-80d7-7e2e37898f00
# ╠═a1924c0b-d0a9-44e9-833b-a6216e9e0627
# ╠═a180ec0f-c47d-47cb-82eb-b120a32a7ed6
# ╠═298e45c4-dc70-4088-a523-98d92f496e9b
# ╠═6b589891-16ff-42b8-95b0-a66e2f94c62b
# ╠═44469ce2-173d-48c6-bb00-fea6211b4b21
# ╠═e23bd3ad-f76d-43c2-8a91-3b36a65db389
# ╠═d21a4e66-6621-4c2d-97ff-85f8a6a09c34
# ╠═66278a87-23c2-4354-9239-2c9a577ea68f
# ╠═6bd6d4e9-0e13-420a-ac5c-cda420c392a4
# ╠═0358b5d4-4574-444e-95de-dcafdfc84bd3
# ╠═c431ca45-3c51-410f-aeff-20d8e5f3c725
# ╠═a144b6bf-9986-4246-9a24-aeaa69b82b80
# ╠═5b3d0a71-00ff-44fe-b004-9d15c92a1c2a
# ╠═c2f95a5d-3711-45af-835f-3770d1612954
# ╠═341c8e71-b73d-451d-bba9-aae64269d01d
# ╠═8ff15860-93a5-4669-9349-f6dd10f97cb7
# ╠═e616c755-4d69-422f-b694-9c76dcaf0334
# ╠═7b959d38-ae82-4355-87ee-5c2f4ff937c2
# ╠═b9a30be5-17d7-4137-af7c-4935c567e0f8
# ╠═3373fef3-08c9-4bb7-9d33-814ce4eedf24
# ╠═81a28886-603a-4618-a98c-d6562fc8f5a1
# ╠═53e4fc32-c61a-4b17-90de-45fe0bb1c233
# ╠═b345f6ed-9224-418c-8076-d69e7e48f01e
# ╠═da49a021-4bfc-43c8-91f5-a93f9ffc4270
# ╠═dd00edef-1382-4f05-a632-18bcbd6d3872
# ╠═f2d4aa10-fda1-4fdb-b134-6526b225c726
# ╠═07aafc3c-a0c5-4f50-8f29-8e686a28b828
# ╠═ae1e97ca-bcc0-4ba9-b878-1f38eafc714f
# ╠═0488a732-a7c0-4531-aeab-d7a0b41701b7
# ╠═02776897-4a63-4066-aff9-e308151f715a
# ╠═4f77b21b-b16a-41f4-bcc0-459e79db8896
# ╠═a47c2d44-b11f-4d5e-810a-dfa83957a683
# ╠═01e7991c-91ce-48b9-969e-5a71bf055f28

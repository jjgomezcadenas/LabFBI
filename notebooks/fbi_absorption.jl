### A Pluto.jl notebook ###
# v0.14.4

using Markdown
using InteractiveUtils

# ╔═╡ 3d85935c-223e-4e79-b2d5-82b907d87def
begin
	import Pkg
	Pkg.add("GLM")
end

# ╔═╡ e4d1f33c-a595-11eb-3987-c9c2518d9bdf
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

# ╔═╡ ba836a5d-f180-4776-b9ef-fc91c44dbea0
using DrWatson

# ╔═╡ f3f9df6f-5948-48a2-8d8f-e08410658aad
using GLM

# ╔═╡ 502081ac-1fd2-4c5a-9997-1ceae6536a50
md"# Absorption as a function of concentration: molar extinction coefficient"

# ╔═╡ d29f573b-38b7-47a6-ad32-9d8eaf5c025b
@quickactivate "LabFBI"

# ╔═╡ 9bd592f0-b5cc-4208-85d8-0d5951eb19c0
projectdir()

# ╔═╡ a5257183-9a75-4c7b-ab29-a382ba2a0fdd
datadir()

# ╔═╡ ed9219c1-44ce-437f-9e6f-dc23a4ff3a2d
srcdir()

# ╔═╡ 13f4d808-e369-4eb6-bd75-14ae943659fb
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

# ╔═╡ f5ecd53e-4e19-4a0d-bce1-c55f111f0c87
dff = ingredients(srcdir("fbi.jl"))

# ╔═╡ c4b2ef8a-8ed1-431a-816b-a05b0067196d
dfp = ingredients(srcdir("plotDF.jl"))

# ╔═╡ e489b62b-05dd-40ba-8d79-58506cb4175b
utils = ingredients(srcdir("utils.jl"))

# ╔═╡ ee5f3f82-630c-4d60-8e57-9f14b66c925e
math = ingredients(srcdir("math.jl"))

# ╔═╡ d0f39119-6d29-47e5-8c97-6b3cf3d254b0
cols = [:green :orange :black :purple :red  :yellow :brown :white]


# ╔═╡ d9c61a65-7372-4bcf-80e6-69879cacb4fa
md"## Notebook"

# ╔═╡ 74bedd71-3c65-4814-8985-5aec22bdb8a0
md"### Main result
- Table with the values of the molar extinction coefficient"

# ╔═╡ ca8e9e14-23d0-40f6-903b-307c61991b17
md"### Results"

# ╔═╡ a0560645-8201-4333-b777-143ca86474a1
md"""
	mec(df, rcols, lambda, title, markersize=4, legend=:topleft, digits=3)

Produce the fit object and relevant plot obects
"""

# ╔═╡ 637cbf91-c779-40f1-a13c-213cf1958467
function mec(df, rcols, lambda, title, markersize=4, legend=:topleft, digits=3)
	cg  = names(df) 
	pcg = dfp.plotdf_xy.( (df,), ("λ",), cg[rcols], 
	                 ("λ (nm)"), ("A (1/cm)",), label=cg[rcols])
	csg, eps, ft = dff.mec_at_lamda(df, lambda, rcols)
	pdt = round.(predict(ft),digits=digits)
	
	pg  = dfp.plot_data_and_model(csg, eps, pdt, "C (M)", "A(1/cm)", 
	                              title, 
		                          markersize=markersize, 
		                          legend=legend)
	return ft, pcg, pg
end

# ╔═╡ 7e510b91-f2e1-438e-aa75-052fefa3a92f
md"### Molar extinction coefficient for FbiG1 at 325 & 405 nm"

# ╔═╡ abca888d-e16d-4a82-8d37-831ce09fc005
g1adf = dff.load_df_from_csv(datadir("fluorimeter/absorption"),
	                                  "FBIG1_coeff_abs_molar.csv",  dff.spG) 


# ╔═╡ 69011717-92d4-478c-a817-d186e3355966
begin
	g1ft325, g1pcg325, g1pg325 = mec(g1adf, 2:5, 325.0, "A vs C FBI G1 (325 nm)");
	g1ft405, g1pcg405, g1pg405 = mec(g1adf, 2:5, 405.0, "A vs C FBI G1 (405 nm)");
end
	

# ╔═╡ ee7b9b20-197c-41f4-a1dc-cdc8b72e69bd
plot(g1pcg325..., layout = (2, 2), legend = false, fmt = :png)

# ╔═╡ f9fb0f66-ac7c-4c7e-af0f-1f7b1d928d72
plot(g1pg325, g1pg405, layout = (1, 2), legend=:topleft, fmt = :png)

# ╔═╡ 5d2a9080-330b-4052-83d9-40d72bd6a8db
md" #### Molar extinction coefficient for FbiG1:
- For λ=325 nm 
- ϵ  = $(utils.float_to_fstr(coef(g1ft325)[2], \"%7.3g\")) M-¹ cm-¹
- σϵ = $(utils.float_to_fstr(stderror(g1ft325)[2], \"%7.3g\")) M-¹ cm-¹

- For λ=405 nm 
- ϵ  = $(utils.float_to_fstr(coef(g1ft405)[2], \"%7.3g\")) M-¹ cm-¹
- σϵ = $(utils.float_to_fstr(stderror(g1ft405)[2], \"%7.3g\")) M-¹ cm-¹
"

# ╔═╡ d7735076-eb30-44a9-ae2d-8dde1872e7d9
md"### Molar extinction coefficient for FbiBaG1 at 325 & 405 nm"

# ╔═╡ 9d23a7f5-73e5-490a-8b40-2a001eedd260
bag1adf = dff.load_df_from_csv(datadir("fluorimeter/absorption"),
	                                  "FBIBaG1_coeff_abs_molar.csv",  dff.spG) 


# ╔═╡ 57aba06e-e9bb-4b41-816a-64e2e06cd712
begin
	bag1ft325, bag1pcg325, bag1pg325 = mec(bag1adf, 2:5, 325.0, 
	                                   "A vs C FBIBa G1 (325 nm)");
	bag1ft405, bag1pcg405, bag1pg405 = mec(bag1adf, 2:5, 405.0, 
	                                    "A vs C FBI G1 (405 nm)");
end

# ╔═╡ 49580657-7568-45b0-9132-9830ce8bb47c
plot(bag1pcg325..., layout = (2, 2), legend = false, fmt = :png)

# ╔═╡ dedf9713-8d8f-4560-8652-0ed7174b7af7
plot(bag1pg325, bag1pg405, layout = (1, 2), legend=:topleft, fmt = :png)

# ╔═╡ 9cfbf0ec-c4a8-4d73-a55b-67ac54459f59
md" #### Molar extinction coefficient for FbiG1:
- For λ=325 nm 
- ϵ  = $(utils.float_to_fstr(coef(bag1ft325)[2], \"%7.3g\")) M-¹ cm-¹
- σϵ = $(utils.float_to_fstr(stderror(bag1ft325)[2], \"%7.3g\")) M-¹ cm-¹

- For λ=405 nm 
- ϵ  = $(utils.float_to_fstr(coef(bag1ft405)[2], \"%7.3g\")) M-¹ cm-¹
- σϵ = $(utils.float_to_fstr(stderror(bag1ft405)[2], \"%7.3g\")) M-¹ cm-¹
"

# ╔═╡ 83729856-82a7-4294-9ffb-09b9bd66852d
md"### Molar extinction coefficient for FbiG2"

# ╔═╡ 6c5bf83f-c640-41b5-9095-2b3f7958823c
g2adf = dff.load_df_from_csv(datadir("fluorimeter/absorption"),
	                                  "FBIG2_coeff_abs_molar.csv",  dff.spG) 


# ╔═╡ 869ec1a4-a239-4fde-828c-88a570a14f4b
begin
	g2ft325, g2pcg325, g2pg325 = mec(g2adf, 2:7, 325.0, "A vs C FBI G2 (325 nm)");
	g2ft500, g2pcg500, g2pg500 = mec(g2adf, 2:7, 405.0, "A vs C FBI G2 (500 nm)");
end

# ╔═╡ 8a457bab-f999-49b7-8a9a-e060026fe9f9
plot(g2pcg325..., layout = (3, 2), legend = false, fmt = :png)

# ╔═╡ 0d84f24f-b222-4c67-a17d-835d5fdb1c84
plot(g2pg325, g2pg500, layout = (1, 2), legend=:topleft, fmt = :png)

# ╔═╡ b312a638-cac9-4bbf-85c5-ec1cd48db04e
md"##### Response at 5e-04 non linear"

# ╔═╡ ec0b705e-06a1-4572-83a7-284deedf179b
begin
	g2ft325x, g2pcg325x, g2pg325x = mec(g2adf, 3:6, 325.0, "A vs C FBI G2 (325 nm)");
	g2ft500x, g2pcg500x, g2pg500x = mec(g2adf, 3:6, 405.0, "A vs C FBI G2 (500 nm)");
end

# ╔═╡ 0b1c941b-0683-453c-8cdc-e37d01db7682
plot(g2pcg325x..., layout = (2, 2), legend = false, fmt = :png)

# ╔═╡ a3de88ba-5e81-41cf-b79c-a325aa54dc4b
plot(g2pg325x, g2pg500x, layout = (1, 2), legend=:topleft, fmt = :png)

# ╔═╡ 88e74cf6-41a4-4888-a327-a9f1cf701e04
md" #### Molar extinction coefficient for FbiG2:
- Linearity is lost for 5E-04 M
- Fluorimeter sensitivity near 1E-07M
- For λ=325 nm 
- ϵ  = $(utils.float_to_fstr(coef(g2ft325x)[2], \"%7.3g\")) M-¹ cm-¹
- σϵ = $(utils.float_to_fstr(stderror(g2ft325x)[2], \"%7.3g\")) M-¹ cm-¹

- For λ=500 nm 
- ϵ  = $(utils.float_to_fstr(coef(g2ft500x)[2], \"%7.3g\")) M-¹ cm-¹
- σϵ = $(utils.float_to_fstr(stderror(g2ft500x)[2], \"%7.3g\")) M-¹ cm-¹
"

# ╔═╡ ede87f21-44e6-46c2-a783-7a68bc323a58
md"### Molar extinction coefficient for FbiBaG2"

# ╔═╡ e2eac5bb-7f16-4d9e-8c2c-a4107fb38e68
bag2adf = dff.load_df_from_csv(datadir("fluorimeter/absorption"),
	                                  "FBIBaG2_coeff_abs_molar.csv",  dff.spG) 


# ╔═╡ f2807e07-0d75-4927-ac72-2bf57f542b1d
begin
	bag2ft325, bag2pcg325, bag2pg325 = mec(bag2adf, 2:5, 325.0, 
	                                   "A vs C FBIBa G1 (325 nm)");
	bag2ft500, bag2pcg500, bag2pg500 = mec(bag2adf, 2:5, 500.0, 
	                                    "A vs C FBI G1 (500 nm)");
end

# ╔═╡ d61f0e04-c246-4cf8-a543-40d0ff7a36ff
plot(bag2pcg325..., layout = (2, 2), legend = false, fmt = :png)

# ╔═╡ b72ba21f-a4b9-4e22-badf-18fc3ccd4faf
plot(bag2pg325, bag2pg500, layout = (1, 2), legend=:topleft, fmt = :png)

# ╔═╡ 886ea5cc-561f-453b-beda-be1325f929c3
md" #### Molar extinction coefficient for FbiBaG2:
- Linearity is lost for 5E-04 M
- Fluorimeter sensitivity near 1E-07M
- For λ=325 nm 
- ϵ  = $(utils.float_to_fstr(coef(bag2ft325)[2], \"%7.3g\")) M-¹ cm-¹
- σϵ = $(utils.float_to_fstr(stderror(bag2ft325)[2], \"%7.3g\")) M-¹ cm-¹

- For λ=500 nm 
- ϵ  = $(utils.float_to_fstr(coef(bag2ft500)[2], \"%7.3g\")) M-¹ cm-¹
- σϵ = $(utils.float_to_fstr(stderror(bag2ft500)[2], \"%7.3g\")) M-¹ cm-¹
"

# ╔═╡ 3f7ad328-38e0-4f52-8d7c-8816681b0c7c
begin
	ϵ325  = utils.float_to_fstr(coef(g1ft325)[2], "%7.3g")
	σϵ    = utils.float_to_fstr(stderror(g1ft325)[2], "%7.3g")
	ϵ325ba  = utils.float_to_fstr(coef(bag1ft325)[2], "%7.3g")
	σϵba    = utils.float_to_fstr(stderror(bag1ft325)[2], "%7.3g")
	ϵ325g2  = utils.float_to_fstr(coef(g2ft325x)[2], "%7.3g")
	σϵg2 = utils.float_to_fstr(stderror(g2ft325x)[2], "%7.3g")
	ϵ325bag2  = utils.float_to_fstr(coef(bag2ft325)[2], "%7.3g")
	σϵbag2    = utils.float_to_fstr(stderror(bag2ft325)[2], "%7.3g")
end

# ╔═╡ 7c00f635-5c6e-48cd-acdd-77b5079234c1
begin
	ϵ405 = utils.float_to_fstr(coef(g1ft405)[2], "%7.3g")
	σϵ405    = utils.float_to_fstr(stderror(g1ft405)[2], "%7.3g")
	ϵ405ba  = utils.float_to_fstr(coef(bag1ft405)[2], "%7.3g")
	σϵba405    = utils.float_to_fstr(stderror(bag1ft405)[2], "%7.3g")
	ϵ500g2  = utils.float_to_fstr(coef(g2ft500x)[2], "%7.3g")
	σϵg2500 = utils.float_to_fstr(stderror(g2ft500x)[2], "%7.3g")
	ϵ500bag2  = utils.float_to_fstr(coef(bag2ft500)[2], "%7.3g")
	σϵbag2500    = utils.float_to_fstr(stderror(bag2ft500)[2], "%7.3g")
end

# ╔═╡ 07cea1d6-8774-47a0-9541-eef1a9afbabc
md"""
| λ (nm) | Fbi G1 | FbiBa G1 | Fbi G2 | FbiBa G2
|:---------- | ---------- |:------------:|:------------:|:------------:|
| 325    | $ϵ325 +- $σϵ| $ϵ325ba +- $σϵba | $ϵ325g2 +- $σϵg2| $ϵ325bag2 +- $σϵbag2|
| 405 (500)    | $ϵ405 +- $σϵ405| $ϵ405ba +- $σϵba405  | $ϵ500g2 +- $σϵg2500| $ϵ500bag2 +- $σϵbag2500|
"""

# ╔═╡ 3307dc4c-69ec-46c8-97db-3d635213fe28
md"""
| λ (nm) | Fbi G1 | FbiBa G1 | Fbi G2 | FbiBa G2
|:---------- | ---------- |:------------:|:------------:|:------------:|
| 325    | $ϵ325 +- $σϵ| $ϵ325ba +- $σϵba | $ϵ325g2 +- $σϵg2| $ϵ325bag2 +- $σϵbag2|
| 405 (500)    | $ϵ405 +- $σϵ405| $ϵ405ba +- $σϵba405  | $ϵ500g2 +- $σϵg2500| $ϵ500bag2 +- $σϵbag2500|
"""

# ╔═╡ 25dbb38f-48d0-4ffb-b1a3-a6869a676d17
begin
	ϵ325f  = coef(g1ft325)[2]
	σϵf    = stderror(g1ft325)[2]
	ϵ325baf  = coef(bag1ft325)[2]
	σϵbaf    = stderror(bag1ft325)[2]
	ϵ325g2f  = coef(g2ft325x)[2]
	σϵg2f = stderror(g2ft325x)[2]
	ϵ325bag2f  = coef(bag2ft325)[2]
	σϵbag2f    = stderror(bag2ft325)[2]
end

# ╔═╡ e7a85d15-5f74-4c80-9aa9-8e899418dd60
begin
	ϵ405f = coef(g1ft405)[2]
	σϵ405f    = stderror(g1ft405)[2]
	ϵ405baf  = coef(bag1ft405)[2]
	σϵba405f    = stderror(bag1ft405)[2]
	ϵ500g2f  = coef(g2ft500x)[2]
	σϵg2500f = stderror(g2ft500x)[2]
	ϵ500bag2f  = coef(bag2ft500)[2]
	σϵbag2500f   = stderror(bag2ft500)[2]
end

# ╔═╡ 2853ce5f-270b-4e62-9d6c-317458db8461
adf = DataFrame(λ = ["325", "405(500)"], 
	            ϵFbiG1 = [ϵ325f, ϵ405f],
				ϵFbiBaG1 = [ϵ325baf, ϵ405baf],
	            ϵFbiG2 = [ϵ325g2f, ϵ500g2f],
				ϵFbiBaG2 = [ϵ325bag2f, ϵ500bag2f],
	            σFbiG1 = [σϵf, σϵ405f],
				σFbiBaG1 = [σϵbaf, σϵba405f],
	            σFbiG2 = [σϵg2f, σϵg2500f],
				σFbiBaG2 = [σϵbag2f, σϵbag2500f]
)

# ╔═╡ 65367886-c8c7-45fb-90cc-4fe65717a2af
begin
	path = datadir("fbi/molar_extinction_coefficient_G1G2.csv")
	CSV.write(path, adf)
end

# ╔═╡ Cell order:
# ╠═3d85935c-223e-4e79-b2d5-82b907d87def
# ╠═e4d1f33c-a595-11eb-3987-c9c2518d9bdf
# ╠═ba836a5d-f180-4776-b9ef-fc91c44dbea0
# ╠═f3f9df6f-5948-48a2-8d8f-e08410658aad
# ╠═502081ac-1fd2-4c5a-9997-1ceae6536a50
# ╠═d29f573b-38b7-47a6-ad32-9d8eaf5c025b
# ╠═9bd592f0-b5cc-4208-85d8-0d5951eb19c0
# ╠═a5257183-9a75-4c7b-ab29-a382ba2a0fdd
# ╠═ed9219c1-44ce-437f-9e6f-dc23a4ff3a2d
# ╟─13f4d808-e369-4eb6-bd75-14ae943659fb
# ╠═f5ecd53e-4e19-4a0d-bce1-c55f111f0c87
# ╠═c4b2ef8a-8ed1-431a-816b-a05b0067196d
# ╠═e489b62b-05dd-40ba-8d79-58506cb4175b
# ╠═ee5f3f82-630c-4d60-8e57-9f14b66c925e
# ╠═d0f39119-6d29-47e5-8c97-6b3cf3d254b0
# ╠═d9c61a65-7372-4bcf-80e6-69879cacb4fa
# ╠═74bedd71-3c65-4814-8985-5aec22bdb8a0
# ╠═07cea1d6-8774-47a0-9541-eef1a9afbabc
# ╠═ca8e9e14-23d0-40f6-903b-307c61991b17
# ╠═a0560645-8201-4333-b777-143ca86474a1
# ╠═637cbf91-c779-40f1-a13c-213cf1958467
# ╠═7e510b91-f2e1-438e-aa75-052fefa3a92f
# ╠═abca888d-e16d-4a82-8d37-831ce09fc005
# ╠═69011717-92d4-478c-a817-d186e3355966
# ╠═ee7b9b20-197c-41f4-a1dc-cdc8b72e69bd
# ╠═f9fb0f66-ac7c-4c7e-af0f-1f7b1d928d72
# ╟─5d2a9080-330b-4052-83d9-40d72bd6a8db
# ╠═d7735076-eb30-44a9-ae2d-8dde1872e7d9
# ╠═9d23a7f5-73e5-490a-8b40-2a001eedd260
# ╠═57aba06e-e9bb-4b41-816a-64e2e06cd712
# ╠═49580657-7568-45b0-9132-9830ce8bb47c
# ╠═dedf9713-8d8f-4560-8652-0ed7174b7af7
# ╟─9cfbf0ec-c4a8-4d73-a55b-67ac54459f59
# ╠═83729856-82a7-4294-9ffb-09b9bd66852d
# ╠═6c5bf83f-c640-41b5-9095-2b3f7958823c
# ╠═869ec1a4-a239-4fde-828c-88a570a14f4b
# ╠═8a457bab-f999-49b7-8a9a-e060026fe9f9
# ╠═0d84f24f-b222-4c67-a17d-835d5fdb1c84
# ╠═b312a638-cac9-4bbf-85c5-ec1cd48db04e
# ╠═ec0b705e-06a1-4572-83a7-284deedf179b
# ╠═0b1c941b-0683-453c-8cdc-e37d01db7682
# ╠═a3de88ba-5e81-41cf-b79c-a325aa54dc4b
# ╠═88e74cf6-41a4-4888-a327-a9f1cf701e04
# ╠═ede87f21-44e6-46c2-a783-7a68bc323a58
# ╠═e2eac5bb-7f16-4d9e-8c2c-a4107fb38e68
# ╠═f2807e07-0d75-4927-ac72-2bf57f542b1d
# ╠═d61f0e04-c246-4cf8-a543-40d0ff7a36ff
# ╠═b72ba21f-a4b9-4e22-badf-18fc3ccd4faf
# ╟─886ea5cc-561f-453b-beda-be1325f929c3
# ╠═3f7ad328-38e0-4f52-8d7c-8816681b0c7c
# ╠═7c00f635-5c6e-48cd-acdd-77b5079234c1
# ╠═3307dc4c-69ec-46c8-97db-3d635213fe28
# ╠═25dbb38f-48d0-4ffb-b1a3-a6869a676d17
# ╠═e7a85d15-5f74-4c80-9aa9-8e899418dd60
# ╠═2853ce5f-270b-4e62-9d6c-317458db8461
# ╠═65367886-c8c7-45fb-90cc-4fe65717a2af

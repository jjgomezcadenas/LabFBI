### A Pluto.jl notebook ###
# v0.14.4

using Markdown
using InteractiveUtils

# ╔═╡ 0af7907c-a991-11eb-04c7-7f5f1225a448
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

# ╔═╡ a1782816-94e0-4ae0-8083-ed2782b5804d
using DrWatson

# ╔═╡ 1281af5d-937a-47d7-889f-3f2297875236
md"# Test and documentation of functions in module `filters.jl`"


# ╔═╡ d3ce7f8a-cb24-4c05-bdad-d7088470ac71
@quickactivate "LabFBI"

# ╔═╡ fbf6f68d-10d6-44f2-a78a-3cb6357d9e7d
projectdir()

# ╔═╡ f949885c-85b1-4bb8-a796-4a32a1f4ba95
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

# ╔═╡ 547be1f2-93ee-402e-bc2d-c13653a61370
ft = ingredients(srcdir("filters.jl"))

# ╔═╡ df5540b0-7e37-482c-b93b-2af64a879e9c
ut = ingredients(srcdir("utils.jl"))

# ╔═╡ f7caa7c7-2e7c-425e-b721-105428e76d9b
fdf = ingredients(srcdir("dffunctions.jl"))

# ╔═╡ e8e8f7c7-090f-4101-b994-1cd12b5e12a3
pf = ingredients(srcdir("plotFilters.jl"))

# ╔═╡ fec2a2bd-44c3-481a-b188-f47296aca76e
md"## Notebook"

# ╔═╡ 66ee5b4b-1c6c-484a-841d-06485dc9c12d
function vect_to_list(vect::AbstractVector, fmt::String)
	vs = ut.to_fstr.(vect,(fmt,))
	str = ""
	for v in vs[1:end-1]
		str = string(str," - ", v,"\n")
	end
	str = string(str," - ", vs[end])
	str
end


# ╔═╡ 5bcd5d64-7d79-4d1c-a530-2e68c26affd5
md"- #### Filter names"

# ╔═╡ fdbc6713-b604-4c82-b327-b08d5efb72d3
Text(vect_to_list(ft.fnames, "%s"))

# ╔═╡ 6131a5be-d918-47f0-8521-d60b10351c04
md"- #### Filter Description"

# ╔═╡ f9f10cd1-b4f9-42ee-8628-6a0ffc569f43
Text(vect_to_list(ft.fdoc, "%s"))

# ╔═╡ c1253cc3-7fb2-4b58-9104-937fec7032a9
md"- #### Filter ranges"

# ╔═╡ 580b9151-11e5-4e15-9aa4-b11698432d05
ft.franges

# ╔═╡ e8109d9a-88e1-410a-b35e-4e1410b14b22
md"""
	path_from_name(name::String, path::String, ext::String=".csv")
Return a full file path (and extension) fron name

\# Arguments
- `name::String`     : name of the filter.
- `path::String`     : a path to the directory with filter data.
- `ext::String=csv`  : file extension

"""

# ╔═╡ f5455f20-9f3e-4262-829c-da9db557d9f6
ft.path_from_name("BandPass405", datadir("filters"))

# ╔═╡ bf0270db-127f-4f81-8e94-27179fd74a29
flt = split(ft.path_from_name("BandPass405", datadir("filters")),"/")[end]

# ╔═╡ 859a191a-f20b-4ef6-b385-b2e14726440b
function fname(flt)
	return split(ft.path_from_name(flt, datadir("filters")),"/")[end]
end

# ╔═╡ e97c90e7-a645-46b0-93b8-b63e5f7b2c86
function check_file(path, flt)
	for f in cd(readdir, path)
		if f == flt
			return true 
		end
	end
	return false
end

# ╔═╡ 7d0292e0-f67f-4e69-9b0e-5e7e9d89fc3f
check_file(datadir("filters"), flt)

# ╔═╡ c465a804-4da3-4393-b83e-9ec57e479507
@test all([check_file(datadir("filters"), flt) for flt in fname.(ft.fnames)])

# ╔═╡ b8f6ba40-02d4-4778-b736-47a1c7fec2c5
md"""
	function load_filter(name::String, path::String, ext::String=".csv")

Load a filter dataframe from a csv file

\# Arguments
- `name::String`       : name of the filter.
- `path::String`       : path to the filter folder.
- `ext::String=".csv"` : file extension.
"""

# ╔═╡ ab5f048e-7200-4d3c-8247-d0a209bfbc70
bp405, fbp405 = ft.load_filter("BandPass405", datadir("filters"));

# ╔═╡ b9c2d705-5bd9-43bb-ab1f-094e14235544
pf.plot_filter("BandPass405", bp405, fbp405, ft.filter_ranges["BandPass405"])

# ╔═╡ 473ce2cb-7bc5-44b0-9a9f-0f24d32b588f
lp425, flp425 = ft.load_filter("LongPass425", datadir("filters"));

# ╔═╡ 54bb3b5b-fab2-4dfc-9f15-8db0b01c3e39
pf.plot_filter("LongPass425", lp425, flp425, ft.filter_ranges["LongPass425"])

# ╔═╡ 57c41720-8824-4c15-ac2f-3fa05b966e06
bp430, fbp430 = ft.load_filter("BandPass430", datadir("filters"));

# ╔═╡ 1f9cc656-3934-466d-a5bb-899f68601a68
pf.plot_filter("BandPass430", bp430, fbp430, ft.filter_ranges["BandPass430"])

# ╔═╡ 677d1628-2036-4760-b894-8fb8caf4e4fc
lp450, flp450 = ft.load_filter("LongPass450", datadir("filters"));

# ╔═╡ 83ce923a-d85f-44af-8dcf-4cb594970869
pf.plot_filter("LongPass450", lp450, flp450, ft.filter_ranges["LongPass450"])

# ╔═╡ c47feabf-cc6a-4ca1-8390-48cb7538e163
dn405_522, fdn405_522 = ft.load_filter("DoubleNotch405_522", datadir("filters"));

# ╔═╡ 471a489b-b9f5-4367-9fbd-ebc0055d9324
pf.plot_filter("DoubleNotch405_522", dn405_522, fdn405_522, ft.filter_ranges["DoubleNotch405_522"])

# ╔═╡ 15408f9e-ae45-4cc8-a996-4ad35a5d3975
n405, fn405 = ft.load_filter("NF405", datadir("filters"));

# ╔═╡ a00380fe-b91f-46c9-a16b-002aef6c12f0
pf.plot_filter("NF405", n405, fn405, ft.filter_ranges["NF405"])

# ╔═╡ df146c61-9a89-4dff-a63c-ba5c09eb8f4f
md"- The length of all the DFs should be identical to the ranges defined between start and end λ but some times one λ is missing. Thus, one missing value is allowed"

# ╔═╡ 0a7ab50d-1cc6-418c-a221-5d1313c1fe1f
fdfl = [length(f[!, "λ"]) for f in [bp405,lp425, bp430, lp450, dn405_522, n405]]
	

# ╔═╡ 07953971-9afd-4256-a7b3-422588824dee
fdfl2 = [length(collect(ft.franges[i])) for (i,r) in enumerate(ft.franges)]

# ╔═╡ 9d25e91b-af7c-4289-9140-015f92699edf
dx = abs.(fdfl .-fdfl2)

# ╔═╡ d708c801-004c-4c9d-a1c9-3a418297b5d1
@test length(dx[dx.>1]) == 0

# ╔═╡ cc4893bf-93c4-4163-8488-d243356a7f78

Filters = ft.load_filter.(ft.fnames, (datadir("filters"),));

# ╔═╡ 8e417654-0057-43dd-bb21-83c1fc6d7e12
fdfs, fints = collect(zip(Filters...));

# ╔═╡ 52416da8-fca6-4e89-ade8-0f2bc53166d6
fdflt = [length(f[!, "λ"]) for f in fdfs]


# ╔═╡ a80be897-c1c6-4b54-8ad6-661324e5a89b
bp405p = pf.plot_filterset(fbp405, "BandPass405", 390:420);

# ╔═╡ 1e6b8d94-84cc-4ffa-ac04-014fe7eada3c
bp430p = pf.plot_filterset(fbp430, "BandPass430", 420:440);

# ╔═╡ 5c3a274e-d0e1-4ca4-8538-eab7dc98f345
lp425p = pf.plot_filterset(flp425, "LongPass425", 420:520);

# ╔═╡ 3bdab1fd-9248-477f-b8d2-7d543fb2d99f
lp450p = pf.plot_filterset(flp450, "LongPass450", 440:540);

# ╔═╡ 2c252a08-b6ef-42dd-a16b-b38e08e53dcc
n405p = pf.plot_filterset(fn405, "NF405", 390:420);

# ╔═╡ 0d0d0f28-2c23-4f0e-a481-519b7ad90d9c
dn405_522p = pf.plot_filterset(fdn405_522, "DoubleNotch405_522", 390:450);

# ╔═╡ 9ac02495-860d-4953-8cb7-35a17095957f
plot(bp405p, bp430p, lp425p, lp450p, n405p, dn405_522p, layout = (3, 2), legend=false, fmt = :png)


# ╔═╡ Cell order:
# ╠═0af7907c-a991-11eb-04c7-7f5f1225a448
# ╠═a1782816-94e0-4ae0-8083-ed2782b5804d
# ╠═1281af5d-937a-47d7-889f-3f2297875236
# ╠═d3ce7f8a-cb24-4c05-bdad-d7088470ac71
# ╠═fbf6f68d-10d6-44f2-a78a-3cb6357d9e7d
# ╠═f949885c-85b1-4bb8-a796-4a32a1f4ba95
# ╠═547be1f2-93ee-402e-bc2d-c13653a61370
# ╠═df5540b0-7e37-482c-b93b-2af64a879e9c
# ╠═f7caa7c7-2e7c-425e-b721-105428e76d9b
# ╠═e8e8f7c7-090f-4101-b994-1cd12b5e12a3
# ╠═fec2a2bd-44c3-481a-b188-f47296aca76e
# ╠═66ee5b4b-1c6c-484a-841d-06485dc9c12d
# ╠═5bcd5d64-7d79-4d1c-a530-2e68c26affd5
# ╠═fdbc6713-b604-4c82-b327-b08d5efb72d3
# ╠═6131a5be-d918-47f0-8521-d60b10351c04
# ╠═f9f10cd1-b4f9-42ee-8628-6a0ffc569f43
# ╠═c1253cc3-7fb2-4b58-9104-937fec7032a9
# ╠═580b9151-11e5-4e15-9aa4-b11698432d05
# ╟─e8109d9a-88e1-410a-b35e-4e1410b14b22
# ╠═f5455f20-9f3e-4262-829c-da9db557d9f6
# ╠═bf0270db-127f-4f81-8e94-27179fd74a29
# ╠═859a191a-f20b-4ef6-b385-b2e14726440b
# ╠═e97c90e7-a645-46b0-93b8-b63e5f7b2c86
# ╠═7d0292e0-f67f-4e69-9b0e-5e7e9d89fc3f
# ╠═c465a804-4da3-4393-b83e-9ec57e479507
# ╟─b8f6ba40-02d4-4778-b736-47a1c7fec2c5
# ╠═ab5f048e-7200-4d3c-8247-d0a209bfbc70
# ╠═b9c2d705-5bd9-43bb-ab1f-094e14235544
# ╠═473ce2cb-7bc5-44b0-9a9f-0f24d32b588f
# ╠═54bb3b5b-fab2-4dfc-9f15-8db0b01c3e39
# ╠═57c41720-8824-4c15-ac2f-3fa05b966e06
# ╠═1f9cc656-3934-466d-a5bb-899f68601a68
# ╠═677d1628-2036-4760-b894-8fb8caf4e4fc
# ╠═83ce923a-d85f-44af-8dcf-4cb594970869
# ╠═c47feabf-cc6a-4ca1-8390-48cb7538e163
# ╠═471a489b-b9f5-4367-9fbd-ebc0055d9324
# ╠═15408f9e-ae45-4cc8-a996-4ad35a5d3975
# ╠═a00380fe-b91f-46c9-a16b-002aef6c12f0
# ╠═df146c61-9a89-4dff-a63c-ba5c09eb8f4f
# ╠═0a7ab50d-1cc6-418c-a221-5d1313c1fe1f
# ╠═07953971-9afd-4256-a7b3-422588824dee
# ╠═9d25e91b-af7c-4289-9140-015f92699edf
# ╠═d708c801-004c-4c9d-a1c9-3a418297b5d1
# ╠═cc4893bf-93c4-4163-8488-d243356a7f78
# ╠═8e417654-0057-43dd-bb21-83c1fc6d7e12
# ╠═52416da8-fca6-4e89-ade8-0f2bc53166d6
# ╠═a80be897-c1c6-4b54-8ad6-661324e5a89b
# ╠═1e6b8d94-84cc-4ffa-ac04-014fe7eada3c
# ╠═5c3a274e-d0e1-4ca4-8538-eab7dc98f345
# ╠═3bdab1fd-9248-477f-b8d2-7d543fb2d99f
# ╠═2c252a08-b6ef-42dd-a16b-b38e08e53dcc
# ╠═0d0d0f28-2c23-4f0e-a481-519b7ad90d9c
# ╠═9ac02495-860d-4953-8cb7-35a17095957f

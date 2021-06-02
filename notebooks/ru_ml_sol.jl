### A Pluto.jl notebook ###
# v0.14.2

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
md"# Ru++: Response in ITO and in Solution"


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
Analyzes the data corresponding to Ru++ in solution and compares it with Ru++ in ITO with different techniques.

In the first data set, three samples are studied.

1. ITO\_Ru++ : Covalent deposition of Rupp directly on ITO
2. ITO\_APTES\_RU++: Covalent deposition of Rupp on APTES
3. Ru++\_1E-5\_sol : Rupp on solution at a concentratoin of 1E-5 M

The conditions of the fluorimeter (opening slits) are the same for the three setups, and tunes up for the solution. Thus, the ITO samples are noisy, due to scarce irradiation.

In the second data set, we consider also three data sets:

1. ITO\_APTES\_Ru++	
2. ITO\_Ru++	
3. ITO\_APTES\_Ru++\_unwashed

The first two samples are described above. The sample _unwashed is similar to (1.) but excess fluorophores where not eliminated by washing. 



"""

# ╔═╡ bcd854b0-9918-478e-8db6-656b87652eed
md"#### Read data
"

# ╔═╡ bc2b045a-e211-4294-b79e-0cb37bc5cd41
rudf = dff.load_df_from_csv(datadir("RuIrPy"),
	                          "Rupp_ITO_solution.csv", 
	                          dff.spG)

# ╔═╡ b81f7da7-0489-4468-a0d1-fa825472ff94
begin
    wf = 450.
    we = 900.
    ws = 2.
    wl = wf:ws:we
	LL = collect(wl)
	println("")
end

# ╔═╡ cc56a6a6-f215-464e-87db-f5bafbd8edc0
begin
	gitoru = lfi.LabFbi.dftogf(wl, rudf, "ITO_Ru++")
	gitoapru = lfi.LabFbi.dftogf(wl, rudf, "ITO_APTES_Ru++")
	gsolru = lfi.LabFbi.dftogf(wl, rudf, "Ru++_1E-5_sol")
end

# ╔═╡ 064d1e19-df72-4b74-8782-4c9fdd888a5a


# ╔═╡ 000d0296-6a59-47c0-9b4f-78a8fc173424


# ╔═╡ 53b0ca40-6111-4463-afc6-47035d8cbcf7
pru =lpi.LabPlots.plotdf_gfs([gitoru,gitoapru,gsolru], 500.0:850.0,
	                          ["ITO_Ru++","ITO_APTES_Ru++","Ru++_1E-5_sol"], 
	                          markercolors, 
	                    "λ (nm)", "I (au)", "Ru++", pdf=false,
		                legend=:topright)

# ╔═╡ 196065ec-caf6-4d09-9276-8dc10c086ad1
prusol =lpi.LabPlots.plotdf_gfs([gsolru], 500.0:850.0,
	                          ["Ru++_1E-5_sol"], 
	                          [:orange], 
	                    "λ (nm)", "I (au)", "Ru++", pdf=false,
		                legend=:topright)

# ╔═╡ 8b67f469-280c-4f6d-acf2-1bb5ebb3f6c8
pitoru =lpi.LabPlots.plotdf_gfs([gitoru], 500.0:800.0,
	                          ["ITO_Ru++"], 
	                          markercolors, 
	                    "λ (nm)", "I (au)", "ITO_Ru++", pdf=false,
		                legend=:topright)

# ╔═╡ 4ff066a1-6766-43da-a0f7-f19c766e7ed3
pitoapru =lpi.LabPlots.plotdf_gfs([gitoapru], 500.0:800.0,
	                          ["ITO_APTES_Ru++"], 
	                          markercolors, 
	                    "λ (nm)", "I (au)", "ITO_APTES_Ru++", pdf=false,
		                legend=:topright)

# ╔═╡ 33b28d8a-4b25-4f49-b2c1-3d2febfe574b
rudf2 = dff.load_df_from_csv(datadir("RuIrPy"),
	                          "Rupp_ITO.csv", 
	                          dff.spG)

# ╔═╡ 0e2fd52f-8edf-4553-a1d8-c4b7634d991b
begin
    wf2 = 420.
    we2 = 900.
    ws2 = 2.
    wl2 = wf2:ws2:we2
	LL2 = collect(wl2)
	println("")
end

# ╔═╡ ea1a3251-187b-45f6-b466-7a2455de69bd
begin
	gitoru2 = lfi.LabFbi.dftogf(wl2, rudf2, "ITO_Ru++")
	gitoapru2 = lfi.LabFbi.dftogf(wl2, rudf2, "ITO_APTES_Ru++")
	gitoapruw = lfi.LabFbi.dftogf(wl2, rudf2, "ITO_APTES_Ru++_unwashed")
end

# ╔═╡ 2b58534a-3a9d-4f36-8e85-204acd2c1d27
pru2 =lpi.LabPlots.plotdf_gfs([gitoapru2,gitoru2,gitoapruw], 500.0:800.0,
	                          ["ITO_APTES_Ru++","ITO_Ru++","ITO_APTES_Ru++_unwashed"], 
	                          markercolors, 
	                    "λ (nm)", "I (au)", "Ru++", pdf=false,
		                legend=:topleft)

# ╔═╡ 6ececd05-01ed-4865-a0ce-9d6ba6378b92
pitoapru2 =lpi.LabPlots.plotdf_gfs([gitoapru2], 500.0:850.0,
	                          ["ITO_APTES_Ru++"], 
	                          [:red], 
	                    "λ (nm)", "I (au)", "ITO_APTES_Ru++", pdf=false,
		                legend=:topright)

# ╔═╡ 2e1d8d9a-d09b-405b-80ea-e5a887cf6b20
pitoru2 =lpi.LabPlots.plotdf_gfs([gitoru2], 500.0:850.0,
	                          ["ITO_Ru++"], 
	                          [:blue], 
	                    "λ (nm)", "I (au)", "ITO_APTES_Ru++", pdf=false,
		                legend=:topright)

# ╔═╡ ee832a2f-a5e0-4067-a57f-811c51ad24e9


# ╔═╡ b670323e-9a0a-46cd-85c4-cf3d72b43dc3
plot(prusol, pitoapru2, layout = (1, 2), legend = :topright, fmt = :png)

# ╔═╡ 416070ec-2c7f-42b0-a4fd-3b8e0def84a2
plot(pitoru2, pitoapru2, layout = (1, 2), legend = :topright, fmt = :png)

# ╔═╡ eae25daf-cdda-4d3d-8eb3-63e934bb326e
lpi.LabPlots.merge_plots!(prusol, pitoapru2)

# ╔═╡ 872954ee-9edc-4b40-b99b-635c12ac2279
lpi.LabPlots.merge_plots!(pitoru2, pitoapru2)

# ╔═╡ Cell order:
# ╠═26badf86-379f-4e66-b3c9-9dc122294a91
# ╠═c5119b0c-66e3-4fcc-b8f4-253b6e2d0c89
# ╠═c60716c4-108c-4c1e-85c1-024520065952
# ╠═f3ec8d6e-a8bf-11eb-1f98-636c338d90e7
# ╠═955d45a1-9d25-4a23-86c2-cc279d14e710
# ╠═0447a463-0ceb-49e9-a565-a47d9e881455
# ╟─173989ad-8f8e-4f7e-a88a-9647bf63c844
# ╠═312e7f8b-9df5-4752-8de8-afff2484dad7
# ╟─fe595465-577c-4352-ae1a-66ce5b5edbf8
# ╟─21712b78-ccf6-47fa-b450-13186d35cd5f
# ╠═14f34300-f974-432d-8a53-91da2a0dc298
# ╠═68c4adce-58c0-4807-b87b-10c4e8d4b9b4
# ╠═4c9c1952-9a1b-4bfa-b665-a05ad399e112
# ╠═bcd854b0-9918-478e-8db6-656b87652eed
# ╠═bc2b045a-e211-4294-b79e-0cb37bc5cd41
# ╠═b81f7da7-0489-4468-a0d1-fa825472ff94
# ╠═cc56a6a6-f215-464e-87db-f5bafbd8edc0
# ╠═064d1e19-df72-4b74-8782-4c9fdd888a5a
# ╠═000d0296-6a59-47c0-9b4f-78a8fc173424
# ╠═53b0ca40-6111-4463-afc6-47035d8cbcf7
# ╠═196065ec-caf6-4d09-9276-8dc10c086ad1
# ╠═8b67f469-280c-4f6d-acf2-1bb5ebb3f6c8
# ╠═4ff066a1-6766-43da-a0f7-f19c766e7ed3
# ╠═33b28d8a-4b25-4f49-b2c1-3d2febfe574b
# ╠═0e2fd52f-8edf-4553-a1d8-c4b7634d991b
# ╠═ea1a3251-187b-45f6-b466-7a2455de69bd
# ╠═2b58534a-3a9d-4f36-8e85-204acd2c1d27
# ╠═6ececd05-01ed-4865-a0ce-9d6ba6378b92
# ╠═2e1d8d9a-d09b-405b-80ea-e5a887cf6b20
# ╠═ee832a2f-a5e0-4067-a57f-811c51ad24e9
# ╠═b670323e-9a0a-46cd-85c4-cf3d72b43dc3
# ╠═416070ec-2c7f-42b0-a4fd-3b8e0def84a2
# ╠═eae25daf-cdda-4d3d-8eb3-63e934bb326e
# ╠═872954ee-9edc-4b40-b99b-635c12ac2279

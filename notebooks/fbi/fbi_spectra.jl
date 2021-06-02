### A Pluto.jl notebook ###
# v0.14.4

using Markdown
using InteractiveUtils

# ╔═╡ 5d5e5072-7ad8-11eb-1559-9763c2eff6a0
begin
	using Pkg
	Pkg.add.(["CSV", "DataFrames", "PlutoUI", "Shapefile", "ZipFile", "LsqFit", 	"Plots","Interpolations","QuadGK"])

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
end

# ╔═╡ f58fd9ea-7b52-11eb-2ad5-c5e14f503d9d
md"# FBI and FBI + Ba spectra on Silica Powder"

# ╔═╡ bbfe5d26-7ce2-11eb-3ddc-0dae0b16d55a
md" ### Globals"

# ╔═╡ cc5be582-7ce2-11eb-0982-73d56ffafc9e
begin
	markershapes = [:circle]
	markercolors = [:green :orange :black :purple :red   :yellow :brown :white]
end

# ╔═╡ baa564ce-7cee-11eb-3436-e7055db2ac84
f_430 = (420, 440)

# ╔═╡ bf4739d8-7ce4-11eb-1f1f-7343c35d34a5
md"### Functions"

# ╔═╡ 2395bfe8-7ce4-11eb-0589-5de38fc62bcc
function read_csv(filename::String)
	csv_data = CSV.File(filename)
	data = DataFrame(csv_data)
	return data
end

# ╔═╡ 90c8d7e6-7ce7-11eb-0b10-b5b481623f24
"""
In a DataFrame one has often "XY" variables, that is, a pair of columns "X" and "Y" which represent correlated variables (e.g, intensity and wavelength). In such cases one often wants the XY maximum, that is, finding the maximum of "Y" (call it ymax) and the corresponding value in "X" (x for which y is ymax). This corresponds, for example, to the wavelength at which the intensity is maximal.
"""
function find_max_xy(df::DataFrame, xc::String, yc::String)
	#find ymax
	ymax, imax = findmax(df[!, yc]) # function findmax returns value and index of max
	x_ymax = df[imax, xc]
	return ymax, x_ymax
end

# ╔═╡ 7b645b80-7b7c-11eb-02ea-5106775dd55b
md" ### Data

Fluorimeter data taken March 3 for FBI and FBI + Ba+"

# ╔═╡ 87cf9b72-7ad8-11eb-3cac-f706589b8d5d
begin
	path="/Users/jj/JuliaProjects/LabFBI/data/exp_raw/EDI_043_FBIpowder_silica"
	ffbi = string(path,"/", "EDI_43_FBI.csv")
	ffbiba = string(path,"/", "EDI_43_FBI_Ba.csv")
end

# ╔═╡ 951d3140-7ad8-11eb-1161-3de7e719f309
begin
csv_fbi = CSV.File(ffbi);
fbidf = DataFrame(csv_fbi);
csv_fbiba = CSV.File(ffbiba);
fbibadf = DataFrame(csv_fbiba);
end

# ╔═╡ 3494f86e-7ced-11eb-1748-6d9c053db821
md"### Interpolate the data"

# ╔═╡ e4451f02-7ceb-11eb-2d11-292d679d1e80
begin
	wf = 300.
	we = 700.
	ws = 2.
	wl = wf:ws:we
end

# ╔═╡ ef77d974-7cea-11eb-1692-018323485ae6
fbi = CubicSplineInterpolation(wl, fbidf.I)

# ╔═╡ 43874a0c-7ced-11eb-1d0c-676d88af3b5d
fbiba = CubicSplineInterpolation(wl, fbibadf.I)

# ╔═╡ a27857be-7cec-11eb-29d8-793f631b2a1b
ifbi = [fbi(x) for x in wl]

# ╔═╡ 996a520c-7ced-11eb-05dd-9ffbb67c055a
ifbiba = [fbiba(x) for x in wl]

# ╔═╡ a02dfb78-7d03-11eb-3cdf-73861d043c13
md" ### Plot data and interpolations"

# ╔═╡ 728bb154-7ce2-11eb-3d4b-f784e172ded1
begin
	plot(fbidf.L, fbidf.I,
		label = "FBI data",
    	shape = :circle,
    	color = :black,
    	markersize = 3, leg=:topleft)
	plot!(wl, ifbi,
		label = "FBI-BA interpolation",
    	color = :green,
    	markersize = 3, leg=:topleft)
	plot!(fbibadf.L, fbibadf.I,
		label = "FBI-BA data",
    	shape = :circle,
    	color = :black,
    	markersize = 3, leg=:topleft)
	plot!(wl, ifbiba,
		label = "FBI-BA interpolation",
    	color = :blue,
    	markersize = 3, leg=:topleft)

	xlabel!("Wavelength (nm)")
	ylabel!("Intensity (au) ")
	title!("FBI fluoremeter response")
end


# ╔═╡ c179a3d6-7d03-11eb-3a4f-87cff0bea714
md" ### Find the maxima in L and I"

# ╔═╡ f867c546-7ce8-11eb-26bd-1bbfdbc973eb
imxFbi, wimxFbi = find_max_xy(fbidf, "L", "I")

# ╔═╡ ccdc0ce4-7ce7-11eb-2b52-472696956add
imxFbiBa, wimxFbiBa = find_max_xy(fbibadf, "L", "I")

# ╔═╡ cf176de8-7d03-11eb-28f2-871837439779
md"### Compute the ratios of intensity"

# ╔═╡ f106b17a-7d03-11eb-208b-492841b3ff1e
md" #### Integrated intensity in blue filter"

# ╔═╡ 93389690-7cec-11eb-036d-fd864e3b76be
I430_fbiba, eps = quadgk(fbiba, f_430[1], f_430[2])

# ╔═╡ 0aeef56e-7cef-11eb-0fb2-c5ecd504abbd
I430_fbi, eps2 = quadgk(fbi, f_430[1], f_430[2])

# ╔═╡ 04e492c8-7d04-11eb-0a83-f349205b785c
md" #### Integrated total intensity"

# ╔═╡ 7c5e69e4-7cef-11eb-07e0-b571995e5b10
It_fbiba, eps3 = quadgk(fbiba, wf, we)

# ╔═╡ 5cdf3bc4-7cf0-11eb-31bf-1fc6de0c29a0
It_fbi, eps4 = quadgk(fbi, wf, we)

# ╔═╡ 20b331c8-7d04-11eb-3a0a-cd384a3ee871
md"#### Ratios

- The fraction of FBI signal w.r.t. FBI+BA signal is 12 %
- The fraction of FBI signal w.r.t. total FBI signal is 1 %
- The fraction of FBO+Ba signal w.r.t. total FBI+Ba is 18 %
"

# ╔═╡ 6a0db9fa-7cec-11eb-3d2b-097823ed716a
R430 = I430_fbi / I430_fbiba

# ╔═╡ 69f45766-7cef-11eb-2735-5755a45f80d6
F430_fbiba = I430_fbiba / It_fbiba

# ╔═╡ 0edf9516-7ceb-11eb-09c9-2f4036de2371
F430_fbi = I430_fbi / It_fbi

# ╔═╡ bec56674-7cf0-11eb-2c9a-f9a4ab1b7722
R2_430 = F430_fbi / F430_fbiba

# ╔═╡ d1547dd4-7cf0-11eb-0168-7f2b100a5a22


# ╔═╡ Cell order:
# ╠═f58fd9ea-7b52-11eb-2ad5-c5e14f503d9d
# ╠═5d5e5072-7ad8-11eb-1559-9763c2eff6a0
# ╠═bbfe5d26-7ce2-11eb-3ddc-0dae0b16d55a
# ╠═cc5be582-7ce2-11eb-0982-73d56ffafc9e
# ╠═baa564ce-7cee-11eb-3436-e7055db2ac84
# ╠═bf4739d8-7ce4-11eb-1f1f-7343c35d34a5
# ╠═2395bfe8-7ce4-11eb-0589-5de38fc62bcc
# ╠═90c8d7e6-7ce7-11eb-0b10-b5b481623f24
# ╠═7b645b80-7b7c-11eb-02ea-5106775dd55b
# ╠═87cf9b72-7ad8-11eb-3cac-f706589b8d5d
# ╠═951d3140-7ad8-11eb-1161-3de7e719f309
# ╠═3494f86e-7ced-11eb-1748-6d9c053db821
# ╠═e4451f02-7ceb-11eb-2d11-292d679d1e80
# ╠═ef77d974-7cea-11eb-1692-018323485ae6
# ╠═43874a0c-7ced-11eb-1d0c-676d88af3b5d
# ╠═a27857be-7cec-11eb-29d8-793f631b2a1b
# ╠═996a520c-7ced-11eb-05dd-9ffbb67c055a
# ╠═a02dfb78-7d03-11eb-3cdf-73861d043c13
# ╠═728bb154-7ce2-11eb-3d4b-f784e172ded1
# ╠═c179a3d6-7d03-11eb-3a4f-87cff0bea714
# ╠═f867c546-7ce8-11eb-26bd-1bbfdbc973eb
# ╠═ccdc0ce4-7ce7-11eb-2b52-472696956add
# ╠═cf176de8-7d03-11eb-28f2-871837439779
# ╠═f106b17a-7d03-11eb-208b-492841b3ff1e
# ╠═93389690-7cec-11eb-036d-fd864e3b76be
# ╠═0aeef56e-7cef-11eb-0fb2-c5ecd504abbd
# ╠═04e492c8-7d04-11eb-0a83-f349205b785c
# ╠═7c5e69e4-7cef-11eb-07e0-b571995e5b10
# ╠═5cdf3bc4-7cf0-11eb-31bf-1fc6de0c29a0
# ╠═20b331c8-7d04-11eb-3a0a-cd384a3ee871
# ╠═6a0db9fa-7cec-11eb-3d2b-097823ed716a
# ╠═69f45766-7cef-11eb-2735-5755a45f80d6
# ╠═0edf9516-7ceb-11eb-09c9-2f4036de2371
# ╠═bec56674-7cf0-11eb-2c9a-f9a4ab1b7722
# ╠═d1547dd4-7cf0-11eb-0168-7f2b100a5a22

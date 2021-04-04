### A Pluto.jl notebook ###
# v0.14.0

using Markdown
using InteractiveUtils

# ╔═╡ 4f35dcfc-8644-11eb-1c86-75608b38bf5f
begin
	using Pkg
	Pkg.add.(["Unitful", "UnitfulEquivalences", "Plots","Interpolations"])
end

# ╔═╡ 2e3c7fda-904f-11eb-2988-25604b1caad0
begin
	Pkg.add("PhysicalConstants")
end

# ╔═╡ 4ea80a97-44a6-4205-ae30-eaf239309313
Pkg.add.(["Images", "ImageIO", "ImageView"])

# ╔═╡ 4a239876-4667-48e5-8f9e-619f2b30765a
Pkg.add("Peaks")

# ╔═╡ 3207f446-8643-11eb-37ba-c9aec47fcb8f
begin
	using Markdown
	using InteractiveUtils
	using PlutoUI
	using Test
	using Plots
	using LsqFit
	using Interpolations
	using Images
	using ImageIO
	using ImageView
	using CSV
	using DataFrames
end

# ╔═╡ 5115917a-8644-11eb-19fc-0528741ca75d
begin
	using Unitful
	using UnitfulEquivalences
end

# ╔═╡ 8ba78a98-904f-11eb-1977-2750643d2f9f
using PhysicalConstants.CODATA2018

# ╔═╡ f2faa594-9316-11eb-2451-33e4b4d3a4f5
using StatsPlots

# ╔═╡ 5356699e-93ca-11eb-0afd-e7d0edc4231e
using QuadGK

# ╔═╡ 887cfa56-2fe0-4173-83fd-0ade21c8ff0d
using Peaks

# ╔═╡ fc8b79a2-8728-11eb-2da7-e3ffa3ceef08
using DrWatson

# ╔═╡ 79cfd2fc-9046-11eb-2b13-1b877d57645d
md"# FBIG2

- First studies of new FBIG2 molecule"

# ╔═╡ 0f2f4c78-8729-11eb-2bab-27812ce8c47e
@quickactivate "LabFBI"

# ╔═╡ 5ee27d52-86fd-11eb-365e-9f2b1a095575
;pwd()

# ╔═╡ 621ec96c-86fd-11eb-1c41-379cc17180dc
datadir()

# ╔═╡ 0b1cedc7-8ada-45a5-adef-fbae794dee3e
markercolors = [:green :orange :black :purple :red  :yellow :brown :white]

# ╔═╡ 3c01d7c8-8645-11eb-372d-4d66ae46ae72
md"## Units"

# ╔═╡ 46b5465a-8645-11eb-0291-612455795518
import Unitful:
    nm, μm, mm, cm, m, km, inch, ft, mi,
    ac,
    mg, g, kg,
    Ra, °F, °C, K,
    rad, °,
    ns, μs, ms, s, minute, hr, d, yr, Hz,
    eV,
    μJ, mJ, J,
	mW, W,
    A, N, mol, mmol, V, L, M

# ╔═╡ ee7f63de-904f-11eb-0e8f-1f8ec9d50811
md"## Physical constants"

# ╔═╡ 631833b8-0649-4f45-b1dd-e4a53091d321
md"#### import Avogadro Number"

# ╔═╡ f8d3f48a-904f-11eb-2659-51f7508b434d
import PhysicalConstants.CODATA2018: N_A

# ╔═╡ 9fdae9e2-6d2e-4901-8f48-14764b6800c2
N_A

# ╔═╡ 45bdd766-9131-11eb-31dc-f9ee1ced0db9
md"# Notebook"

# ╔═╡ cadd0620-9315-11eb-228e-135bd066038d
md"## Molecule spectra"

# ╔═╡ 3581a751-8c1b-4ac7-a2d0-d3f993cf3c8e
md"### FBI"

# ╔═╡ d65c78aa-9315-11eb-094c-01765ce6cc2f
begin
	path=string(datadir(),"/fluorimeter")
	selFbiBa =["EDI044_FBI_Ba_round4_em325_375_405.csv",
             "EDI044_FBI_Ba_round5_em325_375_405.csv"]
	selFbi =["EDI044_FBI_round4_em325_375_405.csv",
               "EDI044_FBI_round5_em325_375_405.csv"]

	fsfbi =   [string(path,"/", fbi) for fbi in selFbi]
    fsfbiba = [string(path,"/", fbi) for fbi in selFbiBa]

	csvSelFbi = CSV.File(fsfbi[1]; delim=';', decimal=',')
	fbidf = DataFrame(csvSelFbi);
	csvSelFbiBa = CSV.File(fsfbiba[1]; delim=';', decimal=',')
	fbibadf = DataFrame(csvSelFbiBa);
end

# ╔═╡ 0eb3c064-9316-11eb-3c97-afba4053eb94
function plot_fbi(λex, fbidf, fbibadf, fbiint, fbibaint, wl,
		          labelFbi, labelFbiBa, colFbi, colFbiBa, title, legend)
	L =collect(wl)
	FBI = fbiint.(L)
	FBIBA = fbibaint.(L);
	p1 = plot(fbidf.λ, fbidf[!,λex],
              	   colour = colFbi,
                   shape  = :circle,
               	   label  = labelFbi,
               	   legend = legend)
	
	lfbi = string(labelFbi," Interpol")
	lfbiba = string(labelFbiBa," Interpol")
	
	p2 = plot!(p1, L, FBI,
    	  colour = colFbi,
    	  label  = lfbi,
    	  legend = legend)

	p3 = plot!(p2, fbibadf.λ, fbibadf[!,λex],
               		  colour = colFbiBa,
               		  shape  = :circle,
               		  label  =labelFbiBa,
               		  legend =legend)

	p4 = plot!(p3, L, FBIBA,
    	  colour = colFbiBa,
    	  label  = lfbiba,
    	  legend = legend)


	xlabel!("λ (nm)")
	ylabel!("Intensity(AU)")
	title!(title)
	return p4
end

# ╔═╡ 322216ee-93cb-11eb-272c-b7f5e5d87ee5
function plot_fbin(fbin, fbiban, wl, 
		           labelFbi, labelFbiBa, colFbi, colFbiBa, title, legend)
	L =collect(wl)
	FBI = fbin.(L)
	FBIBA = fbiban.(L);
	p1 = plot(L, FBI,
    		  colour = colFbi,
    		  label  =labelFbi,
    		  legend=legend)

	p2 = plot!(p1, L, FBIBA,
    		   colour = colFbiBa,
    		   label  =labelFbiBa,
    		   legend=legend)

	xlabel!("λ (nm)")
	ylabel!("Intensity(normalised)")
	title!(title)
	return p2
end

# ╔═╡ 8ad489e2-9318-11eb-1531-7515b778e6b8
begin
    wf = 335.
    we = 701.
    ws = 2.
    wl = wf:ws:we
end

# ╔═╡ 9989f79a-9318-11eb-07a3-a9986579342d
begin
	fbi325 = LinearInterpolation(wl, fbidf[!,"E325"]);
	fbiba325 = LinearInterpolation(wl, fbibadf[!,"E325"]);
end

# ╔═╡ 5361bf8e-9318-11eb-1acd-e9a10c2b30bf
pfbiE325G1 = plot_fbi("E325",fbidf, fbibadf,fbi325,fbiba325,wl,"FbiG1", "FbiBaG1", :green, :blue, "FBI/BA G1 325 nm", :topright)

# ╔═╡ 2ecbd182-93ca-11eb-08f6-3d0c1710f9f5
md"#### Normalize to the area of fbi and fbiba"

# ╔═╡ eba136e2-93c9-11eb-21f8-8b36a5158624
fbiE325G1N = quadgk(fbi325, 335., 700.)[1]

# ╔═╡ 6deec574-93ca-11eb-0b5f-1b0d7ca04b61
fbibaE325G1N = quadgk(fbiba325, 335., 700.)[1]

# ╔═╡ 46411c0e-93cd-11eb-1039-8fe43a7679d5
md"`fpdf(fi, xmin, xmax, N)` takes an interpolated function fi, lower (xmin) and upper (xmax) values of interpolation and a normalisation constant and returns a function equal to fi(x) / N in the range xmin-xmax and zero otherwise."

# ╔═╡ 353b7b37-c037-49f9-b149-56541db6cc28
function gfpdf(fi, xmin, xmax, N)
	function fn(x)
		if x < xmin || x > xmax
			return 0
		else
			return fi(x) / N
		end
	end

	return fn
end

# ╔═╡ 566a984a-f0fc-482f-9b04-029cba608aef
fpdf(fi, xmin::Float64, xmax::Float64, N::Float64) = gfpdf(fi, xmin, xmax, N)

# ╔═╡ fee9eace-93cd-11eb-114c-c174cc6637ce
fbipdf = fpdf(fbi325, 335.0, 700.0, 1.0)

# ╔═╡ 4d17046e-ce2e-4863-8231-2ded5722ff44
fbibapdf = fpdf(fbiba325, 335.0, 700.0, 1.0)

# ╔═╡ 2d0e3c95-fcd6-4d96-9ca0-bd8c8be01edb
@test quadgk(fbipdf, 335., 700.)[1] / fbiE325G1N ≈ 1

# ╔═╡ 7738f929-c288-4c90-92e4-17078184c2de
@test quadgk(fbibapdf, 335., 700.)[1] / fbibaE325G1N ≈ 1

# ╔═╡ 76a569c4-93cb-11eb-3bb4-55e0c9221e53
pfbi = plot_fbin(fbipdf, fbibapdf, wl, "FBI", "FBIBa2+", :green, :blue, "FBI G1 (325 nm)", :topleft)

# ╔═╡ 4d5bdc05-0f5d-403f-b303-93e95ec408ca
md"### FBIG2"

# ╔═╡ 1b1f52a8-ddf1-4f16-ab4a-196af7938450
begin
	fsfbig2    =   string(path,"/", "FBIG2.csv") 
	csvFbiG2 = CSV.File(fsfbig2; delim=';', decimal=',')
	fbig2df = DataFrame(csvFbiG2)
end

# ╔═╡ b75e15f4-c70b-4e8e-8786-b96307e69f16
begin
	fsfbibag2    =   string(path,"/", "FBIBaG2.csv") 
	csvFbiBaG2 = CSV.File(fsfbibag2; delim=';', decimal=',')
	fbibag2df = DataFrame(csvFbiBaG2)
end

# ╔═╡ ecb20e55-e9ae-46e4-928e-35a9382a320c
begin
    w2f = 300.
    w2e = 800.
    w2l = w2f:ws:w2e
end

# ╔═╡ 7cd1e229-96d0-4786-bbe4-d63b524e6619
begin
	fbi325g2 = LinearInterpolation(w2l, fbig2df[!,"E325"]);
	fbiba325g2 = LinearInterpolation(w2l, fbibag2df[!,"E325"]);
end

# ╔═╡ bd811db1-28b8-4869-9b8d-d10e5fa0324e
pfbiE325G2 = plot_fbi("E325", fbig2df, fbibag2df,fbi325g2,fbiba325g2,w2l,"FbiG2", "FbiBaG2", :red, :orange, "FBI/BA G2 325 nm", :topright)

# ╔═╡ 93b3e901-84be-462a-96c6-c6cec4604b50
fbiE325G2N = quadgk(fbi325g2, 480., 800.)[1]

# ╔═╡ d3fd15f2-4e32-4f3f-a5fa-dfc6df4721d1
fbibaE325G2N = quadgk(fbiba325g2, 480., 800.)[1]

# ╔═╡ 1db9de6b-4572-48ba-ab3a-84b4ddc56696
fbig2pdf = fpdf(fbi325g2, 480.0, 800.0, 1.0)

# ╔═╡ 7b9c32f2-a314-4523-a27a-5fee782616d0
fbibag2pdf = fpdf(fbiba325g2, 480.0, 800.0, 1.0)

# ╔═╡ dba59cde-7305-4aff-80d2-0d83fc20d3b1
@test quadgk(fbig2pdf, 480., 800.)[1] / fbiE325G2N ≈ 1

# ╔═╡ 845af2e4-791b-44dd-af4b-2f8fcfa0056d
@test quadgk(fbibag2pdf, 480., 800.)[1] / fbibaE325G2N ≈ 1

# ╔═╡ 7dc9afbf-074a-4149-bbae-7b182916e7f2
pfbig2 = plot_fbin(fbig2pdf, fbibag2pdf, 400.0:800.0, "FBI", "FBIBa2+", :red, :orange, "FBI G2 (325 nm)", :topleft)

# ╔═╡ 6a0ba3b4-ed64-4c21-bbcd-ad88685b6453
plot(pfbi, pfbig2, layout = (1, 2), legend = false, fmt = :png)

# ╔═╡ a2cc0f38-f774-489b-98a0-9c713e99b11a
pfbig2v1 = plot_fbin(fbig2pdf, fbibag2pdf, 480.0:570.0, "FBI", "FBIBa2+", :red, :orange, "FBI G2 (325 nm)", :topleft)

# ╔═╡ ebc85de9-e9f1-4be0-8d63-1982006c047e
fbibag2pdf(540.) / fbig2pdf(540.)

# ╔═╡ 2ee3952b-99d4-4e7f-aa92-e32c16aa008a
pfbig1v1 = plot_fbin(fbipdf, fbibapdf, 400.0:450.0, "FBI", "FBIBa2+", :green, :blue, "FBI G1 (325 nm)", :topleft)

# ╔═╡ 95ce45d1-d000-46c4-9996-2aa6db2876cf
md"### Ratios and double ratio for G1 and G2"

# ╔═╡ b6a12dc8-a173-4c57-a30e-5d8cbd7e74ec
md"##### Maximum discrimination region:
- for g1 is between 420 and 430
- for g2 is between 520 and 540"

# ╔═╡ dfda96c5-9e1a-4ff4-bc41-127ab9332dcf
begin
	λming1 = 420
	λmaxg1 = 430
	λming2 = 520
	λmaxg2 = 540
end

# ╔═╡ 895210fa-ab74-4021-8e5f-e5bae0e8fc5b
md"#### Ratios for G2"

# ╔═╡ 51af996b-3b7a-4156-b4bd-6d4ce60ef653
fbig2bp500 = quadgk(fbig2pdf, 520.0,  540.0)[1]

# ╔═╡ 6da39665-c055-4b15-808a-ba154c148967
fbibag2bp500 = quadgk(fbibag2pdf, 520.0,  540.0)[1]

# ╔═╡ 47909f5c-6808-4f37-8893-26c3003ee47f
md"##### Ratio of total recorded signal: FBIBa2+/FBI"

# ╔═╡ 91b3b0e8-833a-453d-a88b-fdb9ba19650b
raG2 = fbibaE325G2N / fbiE325G2N

# ╔═╡ 3b1b2e70-0bb8-4591-a4ef-3aa6fa9d683f
md"##### Fraction of FBI - FBIBa2+ in selected region"

# ╔═╡ fde798e2-9e72-4aa9-a091-fa04b2b920af
xfbibaG2 = fbibag2bp500 / fbibaE325G2N

# ╔═╡ f459b50f-9639-456e-9939-5e13d8e042d2
xfbiG2 = fbig2bp500 / fbiE325G2N

# ╔═╡ a5a186ba-d2e7-4c4d-840c-1798e8ac59bb
md"##### Signal over background in selected region (double ratio)"

# ╔═╡ 1a044409-d13a-4a0e-bc38-5d0d8d4366b8
r2g2 = fbibag2bp500 / fbig2bp500

# ╔═╡ 1a1b4a30-3966-40ab-8469-1509539defa2
md"### Ratios for G1"

# ╔═╡ e0b1c1d0-1d65-48c0-b09f-844bd27187fb
fbig1bp430 = quadgk(fbipdf, 400.0,  430.0)[1]

# ╔═╡ dea8bc12-e165-471f-8689-2d329a92b364
fbibag1bp430 = quadgk(fbibapdf, 400.0,  430.0)[1]

# ╔═╡ d5326ccc-3e0c-4d5a-bada-060236ea521a
md"##### Ratio of total recorded signal: FBIBa2+/FBI"

# ╔═╡ a5069c59-3069-4be5-b154-7a82f8d27e96
raG1 = fbibaE325G1N / fbiE325G1N

# ╔═╡ 68769f50-d271-4495-9747-b09dee800db2
md"##### Fraction of FBI - FBIBa2+ in selected region"

# ╔═╡ 42ab1511-6c45-46f5-b88c-3d7c3b897a0a
xfbibaG1 = fbibag1bp430 / fbibaE325G1N

# ╔═╡ 0d25107d-45c8-4a5c-bfb1-b94c5e7e122e
xfbiG1 = fbig1bp430 / fbiE325G1N

# ╔═╡ 416a6cc5-06d7-4749-ab07-c074f16fa41a
md"##### Signal over background in selected region (double ratio)"

# ╔═╡ e1451597-f05e-4cce-bb2f-b41969638027
r2g1 = fbibag1bp430 / fbig1bp430

# ╔═╡ 09c95810-3983-4007-82a1-275aeb8e07cc
md"### Summary of results:
- Double ratio for G1 r2g1 =$r2g1
- Double ratio for G2 r2g1 =$r2g2
"

# ╔═╡ 4c8107d9-327f-4ee1-bfb9-fb1384c7cc8b


# ╔═╡ Cell order:
# ╟─79cfd2fc-9046-11eb-2b13-1b877d57645d
# ╠═4f35dcfc-8644-11eb-1c86-75608b38bf5f
# ╠═2e3c7fda-904f-11eb-2988-25604b1caad0
# ╠═4ea80a97-44a6-4205-ae30-eaf239309313
# ╠═4a239876-4667-48e5-8f9e-619f2b30765a
# ╠═3207f446-8643-11eb-37ba-c9aec47fcb8f
# ╠═5115917a-8644-11eb-19fc-0528741ca75d
# ╠═8ba78a98-904f-11eb-1977-2750643d2f9f
# ╠═f2faa594-9316-11eb-2451-33e4b4d3a4f5
# ╠═5356699e-93ca-11eb-0afd-e7d0edc4231e
# ╠═887cfa56-2fe0-4173-83fd-0ade21c8ff0d
# ╠═fc8b79a2-8728-11eb-2da7-e3ffa3ceef08
# ╠═0f2f4c78-8729-11eb-2bab-27812ce8c47e
# ╠═5ee27d52-86fd-11eb-365e-9f2b1a095575
# ╠═621ec96c-86fd-11eb-1c41-379cc17180dc
# ╠═0b1cedc7-8ada-45a5-adef-fbae794dee3e
# ╟─3c01d7c8-8645-11eb-372d-4d66ae46ae72
# ╠═46b5465a-8645-11eb-0291-612455795518
# ╠═ee7f63de-904f-11eb-0e8f-1f8ec9d50811
# ╠═631833b8-0649-4f45-b1dd-e4a53091d321
# ╠═f8d3f48a-904f-11eb-2659-51f7508b434d
# ╠═9fdae9e2-6d2e-4901-8f48-14764b6800c2
# ╟─45bdd766-9131-11eb-31dc-f9ee1ced0db9
# ╠═cadd0620-9315-11eb-228e-135bd066038d
# ╠═3581a751-8c1b-4ac7-a2d0-d3f993cf3c8e
# ╠═d65c78aa-9315-11eb-094c-01765ce6cc2f
# ╠═0eb3c064-9316-11eb-3c97-afba4053eb94
# ╠═322216ee-93cb-11eb-272c-b7f5e5d87ee5
# ╠═8ad489e2-9318-11eb-1531-7515b778e6b8
# ╠═9989f79a-9318-11eb-07a3-a9986579342d
# ╠═5361bf8e-9318-11eb-1acd-e9a10c2b30bf
# ╠═2ecbd182-93ca-11eb-08f6-3d0c1710f9f5
# ╠═eba136e2-93c9-11eb-21f8-8b36a5158624
# ╠═6deec574-93ca-11eb-0b5f-1b0d7ca04b61
# ╠═46411c0e-93cd-11eb-1039-8fe43a7679d5
# ╠═353b7b37-c037-49f9-b149-56541db6cc28
# ╠═566a984a-f0fc-482f-9b04-029cba608aef
# ╠═fee9eace-93cd-11eb-114c-c174cc6637ce
# ╠═4d17046e-ce2e-4863-8231-2ded5722ff44
# ╠═2d0e3c95-fcd6-4d96-9ca0-bd8c8be01edb
# ╠═7738f929-c288-4c90-92e4-17078184c2de
# ╠═76a569c4-93cb-11eb-3bb4-55e0c9221e53
# ╠═4d5bdc05-0f5d-403f-b303-93e95ec408ca
# ╠═1b1f52a8-ddf1-4f16-ab4a-196af7938450
# ╠═b75e15f4-c70b-4e8e-8786-b96307e69f16
# ╠═ecb20e55-e9ae-46e4-928e-35a9382a320c
# ╠═7cd1e229-96d0-4786-bbe4-d63b524e6619
# ╠═bd811db1-28b8-4869-9b8d-d10e5fa0324e
# ╠═93b3e901-84be-462a-96c6-c6cec4604b50
# ╠═d3fd15f2-4e32-4f3f-a5fa-dfc6df4721d1
# ╠═1db9de6b-4572-48ba-ab3a-84b4ddc56696
# ╠═7b9c32f2-a314-4523-a27a-5fee782616d0
# ╠═dba59cde-7305-4aff-80d2-0d83fc20d3b1
# ╠═845af2e4-791b-44dd-af4b-2f8fcfa0056d
# ╠═7dc9afbf-074a-4149-bbae-7b182916e7f2
# ╠═6a0ba3b4-ed64-4c21-bbcd-ad88685b6453
# ╠═a2cc0f38-f774-489b-98a0-9c713e99b11a
# ╠═ebc85de9-e9f1-4be0-8d63-1982006c047e
# ╠═2ee3952b-99d4-4e7f-aa92-e32c16aa008a
# ╠═95ce45d1-d000-46c4-9996-2aa6db2876cf
# ╠═b6a12dc8-a173-4c57-a30e-5d8cbd7e74ec
# ╠═dfda96c5-9e1a-4ff4-bc41-127ab9332dcf
# ╠═895210fa-ab74-4021-8e5f-e5bae0e8fc5b
# ╠═51af996b-3b7a-4156-b4bd-6d4ce60ef653
# ╠═6da39665-c055-4b15-808a-ba154c148967
# ╠═47909f5c-6808-4f37-8893-26c3003ee47f
# ╠═91b3b0e8-833a-453d-a88b-fdb9ba19650b
# ╠═3b1b2e70-0bb8-4591-a4ef-3aa6fa9d683f
# ╠═fde798e2-9e72-4aa9-a091-fa04b2b920af
# ╠═f459b50f-9639-456e-9939-5e13d8e042d2
# ╠═a5a186ba-d2e7-4c4d-840c-1798e8ac59bb
# ╠═1a044409-d13a-4a0e-bc38-5d0d8d4366b8
# ╠═1a1b4a30-3966-40ab-8469-1509539defa2
# ╠═e0b1c1d0-1d65-48c0-b09f-844bd27187fb
# ╠═dea8bc12-e165-471f-8689-2d329a92b364
# ╠═d5326ccc-3e0c-4d5a-bada-060236ea521a
# ╠═a5069c59-3069-4be5-b154-7a82f8d27e96
# ╠═68769f50-d271-4495-9747-b09dee800db2
# ╠═42ab1511-6c45-46f5-b88c-3d7c3b897a0a
# ╠═0d25107d-45c8-4a5c-bfb1-b94c5e7e122e
# ╠═416a6cc5-06d7-4749-ab07-c074f16fa41a
# ╠═e1451597-f05e-4cce-bb2f-b41969638027
# ╠═09c95810-3983-4007-82a1-275aeb8e07cc
# ╠═4c8107d9-327f-4ee1-bfb9-fb1384c7cc8b

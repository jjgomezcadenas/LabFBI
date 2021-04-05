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

# ╔═╡ 20a461d5-16ac-4707-adec-ebfce6fc89c6
include(srcdir("fbi.jl"))

# ╔═╡ 79cfd2fc-9046-11eb-2b13-1b877d57645d
md"# FBIG2

- First studies of new FBIG2 molecule"

# ╔═╡ 0f2f4c78-8729-11eb-2bab-27812ce8c47e
@quickactivate "LabFBI"

# ╔═╡ 5ee27d52-86fd-11eb-365e-9f2b1a095575
;pwd()

# ╔═╡ 621ec96c-86fd-11eb-1c41-379cc17180dc
datadir()

# ╔═╡ 9913df39-7aa4-42a3-8cb9-f66707185e22
srcdir("fbi.jl")

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

# ╔═╡ a1af2f43-a6f4-4b5d-a882-abf91f109a44
md"## Notebook functions"

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
               	   legend = legend,
				   fmt = :png)
	
	lfbi = string(labelFbi," Interpol")
	lfbiba = string(labelFbiBa," Interpol")
	
	p2 = plot!(p1, L, FBI,
    	  colour = colFbi,
    	  label  = lfbi,
    	  legend = legend,
	   	  fmt = :png)

	p3 = plot!(p2, fbibadf.λ, fbibadf[!,λex],
               		  colour = colFbiBa,
               		  shape  = :circle,
               		  label  =labelFbiBa,
               		  legend =legend,
				      fmt = :png)

	p4 = plot!(p3, L, FBIBA,
    	  colour = colFbiBa,
    	  label  = lfbiba,
    	  legend = legend,
		  fmt = :png)


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
    		  legend=legend,
		      fmt = :png)

	p2 = plot!(p1, L, FBIBA,
    		   colour = colFbiBa,
    		   label  =labelFbiBa,
    		   legend=legend,
		       fmt = :png)

	xlabel!("λ (nm)")
	ylabel!("Intensity(normalised)")
	title!(title)
	return p2
end

# ╔═╡ 8e63260a-0d0b-4551-96ee-88322e6e389f


# ╔═╡ cadd0620-9315-11eb-228e-135bd066038d
md"## Molecule spectra"

# ╔═╡ 660c8c2d-6b3c-4a12-8e44-ea9dd1c461ff
begin
fbiname      = "EDI044_FBI_round4_em325_375_405.csv"
fbibaname    = "EDI044_FBI_Ba_round4_em325_375_405.csv"
fbiabsname   = "EDI_029_absorption.csv"
fbig2name    = "FBIG2.csv"
fbibag2name  = "FBIBaG2.csv"
	
fbiabsG2name = 	"FBIG2_absorption.csv"
fbidir       = "fluorimeter"
end

# ╔═╡ 9cc0ece2-f2b9-4b54-9603-75e78d564b76
fbidf   = load_fbidf(fbidir, fbiname)

# ╔═╡ d4e70ec0-b6b9-47a0-a575-6d15fc2d25e1
fbig2df = load_fbidf(fbidir, fbig2name)

# ╔═╡ 871a129f-d016-4088-b1e6-c0a7d74e0393
fbibadf = load_fbidf(fbidir, fbibaname)

# ╔═╡ e7d591ac-fd39-45e6-8c61-e5c74cdfd597
fbibag2df = load_fbidf(fbidir, fbibag2name)

# ╔═╡ 5a401dd6-05a4-47c1-af10-d67a4d9b40dc
fbiabsdf = load_fbidf(fbidir, fbiabsname)

# ╔═╡ d2bc643b-c247-4f55-ad67-1992bf25a30f
fbiabsg2df = load_fbidf(fbidir, fbiabsG2name)

# ╔═╡ 82b9a9fc-feaa-4b67-98b8-577c8bdf3378
md" ## G1"

# ╔═╡ 8ad489e2-9318-11eb-1531-7515b778e6b8
begin
    wf = 335.
    we = 701.
    ws = 2.
    wl = wf:ws:we
end

# ╔═╡ e47a9702-710b-44f4-9adf-136e04b3b04e
fbi325 = dftof(wl, fbidf, "E325")

# ╔═╡ 2cb6d5aa-3bdd-4ce6-be0a-dcaac538b86f
fbiba325 = dftof(wl, fbibadf, "E325")

# ╔═╡ 772e3bc5-a169-4c1d-aa22-2ea0720ef707
fbiE325G1N, fbipdf = ftopdf(wl, fbi325)

# ╔═╡ a1bfbae2-b284-417d-953f-bc5a2493d56c
@test quadgk(fbipdf, wl[1], wl[end])[1] ≈ 1

# ╔═╡ 8cf9ece5-c345-4e65-8426-2b71f9f7fea8
fbibaE325G1N, fbibapdf = ftopdf(wl, fbiba325)

# ╔═╡ ff5cdffb-0ce9-4f1b-a0ad-ce74e9d8f5e9
@test quadgk(fbibapdf, wl[1], wl[end])[1] ≈ 1

# ╔═╡ 627d372c-3b7a-45c0-a416-b8ffe4e5a02b
md"## G2"

# ╔═╡ e6565f94-4419-4359-941d-bbddd56381fc
begin
    w2f = 300.
    w2e = 800.
    w2l = w2f:ws:w2e
end

# ╔═╡ 0b07b62b-62d1-43da-949e-77a02d74856c
fbi325g2 = dftof(w2l, fbig2df, "E325")

# ╔═╡ ea30ed3d-a19c-4d3f-8806-737df8954a15
fbiba325g2 = dftof(w2l, fbibag2df, "E325")

# ╔═╡ 4a57223f-d47b-4d3c-b07c-5121585d20d5
fbiE325G2N, fbig2pdf   = ftopdf(w2l, fbi325g2)

# ╔═╡ 4eee17f1-117a-4fb9-88b8-5abd5bc6c69a
fbibaE325G2N, fbibag2pdf = ftopdf(w2l, fbiba325g2)

# ╔═╡ bd9421ed-8fe5-45dd-9f9e-b9cf99077694
@test quadgk(fbig2pdf, w2l[1], w2l[end])[1]  ≈ 1

# ╔═╡ 1f0ec58c-3184-4ceb-b019-28d9bffbc16f
@test quadgk(fbibag2pdf, w2l[1], w2l[end])[1]  ≈ 1

# ╔═╡ c9671526-27fe-41fd-aaf4-a63a5562a726
md"### Absorption"

# ╔═╡ 1481eac4-8e0a-42db-822b-399af881abb4
md"#### G1"

# ╔═╡ c33899c4-ab79-4bc2-9511-c99900a0d2f8
begin
	@df fbiabsdf plot(:λ, [:FBI :FBI_Ba], 
				   colour = [:green :blue],
				   legend=:topleft, fmt = :png)
	xlabel!("λ (nm)")
	ylabel!("Absorption probability (AU)")
	title!("Absorption probability FBI and FBI+BA G1 ")
end

# ╔═╡ dbc8668e-26f4-43e6-9474-29d04108a7ef
md"#### G2"

# ╔═╡ c4d19122-f3cc-4129-ac00-95fe0365d097
begin
	@df fbiabsg2df plot(:λ, [:FBIG2 :BIG2Ba], 
				   colour = [:green :blue],
				   legend=:topleft, fmt = :png)
	xlabel!("λ (nm)")
	ylabel!("Absorption probability (AU)")
	title!("Absorption probability FBI and FBI+BA G2")
end

# ╔═╡ cb03c3d4-c7f8-4247-9a66-99d3deda9fbb
begin
	@df fbiabsg2df plot(:λ[400:600], [:FBIG2[400:600] :BIG2Ba[400:600]], 
				   colour = [:green :blue],
				   legend=:topleft, fmt = :png)
	xlabel!("λ (nm)")
	ylabel!("Absorption probability (AU)")
	title!("Absorption probability FBI and FBI+BA G2")
end

# ╔═╡ 5361bf8e-9318-11eb-1acd-e9a10c2b30bf
pfbiE325G1 = plot_fbi("E325",fbidf, fbibadf,fbi325,fbiba325,wl,"FbiG1", "FbiBaG1", :green, :blue, "FBI/BA G1 325 nm", :topright)

# ╔═╡ 971f16da-62af-49b0-bb38-5f4c39b2b886
pfbi = plot_fbin(fbipdf, fbibapdf, wl, 
		  "Fbi G1 pdf", "FbiBa G1 pdf", :green, :blue, "PDF FBI G1 325 nm", :topright)

# ╔═╡ 4d5bdc05-0f5d-403f-b303-93e95ec408ca
md"### FBIG2"

# ╔═╡ 7cd1e229-96d0-4786-bbe4-d63b524e6619


# ╔═╡ bd811db1-28b8-4869-9b8d-d10e5fa0324e
pfbiE325G2 = plot_fbi("E325", fbig2df, fbibag2df,fbi325g2,fbiba325g2,w2l,"FbiG2", "FbiBaG2", :red, :orange, "FBI/BA G2 325 nm", :topright)

# ╔═╡ 7dc9afbf-074a-4149-bbae-7b182916e7f2
pfbig2 = plot_fbin(fbig2pdf, fbibag2pdf, w2l, "FBI", "FBIBa2+", :red, :orange, "FBI G2 (325 nm)", :topleft)

# ╔═╡ 6a0ba3b4-ed64-4c21-bbcd-ad88685b6453
plot(pfbi, pfbig2, layout = (1, 2), legend = false, fmt = :png)

# ╔═╡ 9de1cf4b-a383-4bfe-a45b-351b94c75385
plot(pfbiE325G1, pfbiE325G2, layout = (1, 2), legend = false, fmt = :png)

# ╔═╡ b6a12dc8-a173-4c57-a30e-5d8cbd7e74ec
md"### Maximum discrimination region:
- for g1 is between 420 and 430
- for g2 is between 520 and 540"

# ╔═╡ a2cc0f38-f774-489b-98a0-9c713e99b11a
pfbig2v1 = plot_fbin(fbig2pdf, fbibag2pdf, 480.0:570.0, "FBI", "FBIBa2+", :red, :orange, "FBI G2 (325 nm)", :topleft)

# ╔═╡ ebc85de9-e9f1-4be0-8d63-1982006c047e
fbibag2pdf(540.) / fbig2pdf(540.)

# ╔═╡ 2ee3952b-99d4-4e7f-aa92-e32c16aa008a
pfbig1v1 = plot_fbin(fbipdf, fbibapdf, 400.0:450.0, "FBI", "FBIBa2+", :green, :blue, "FBI G1 (325 nm)", :topleft)

# ╔═╡ 95ce45d1-d000-46c4-9996-2aa6db2876cf
md"### Ratios and double ratio for G1 and G2"

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
fbig2bp500 = quadgk(fbi325g2, λming2,  λmaxg2)[1]

# ╔═╡ 6da39665-c055-4b15-808a-ba154c148967
fbibag2bp500 = quadgk(fbiba325g2, λming2,  λmaxg2)[1]

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
fbig1bp430 = quadgk(fbi325, 400.0,  430.0)[1]

# ╔═╡ dea8bc12-e165-471f-8689-2d329a92b364
fbibag1bp430 = quadgk(fbiba325, 400.0,  430.0)[1]

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
- ##### Double ratio for G1 
- r2g1 =$r2g1
- ##### Double ratio for G2 
- r2g1 =$r2g2
- ##### Fraction of FBI - FBIBa2+ in selected region for G1:
- signal : $xfbibaG1
- bkgnd  : $xfbiG1
- ##### Fraction of FBI - FBIBa2+ in selected region for G2:
- signal : $xfbibaG2
- bkgnd  : $xfbiG2
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
# ╠═9913df39-7aa4-42a3-8cb9-f66707185e22
# ╠═20a461d5-16ac-4707-adec-ebfce6fc89c6
# ╠═0b1cedc7-8ada-45a5-adef-fbae794dee3e
# ╟─3c01d7c8-8645-11eb-372d-4d66ae46ae72
# ╠═46b5465a-8645-11eb-0291-612455795518
# ╠═ee7f63de-904f-11eb-0e8f-1f8ec9d50811
# ╠═631833b8-0649-4f45-b1dd-e4a53091d321
# ╠═f8d3f48a-904f-11eb-2659-51f7508b434d
# ╠═9fdae9e2-6d2e-4901-8f48-14764b6800c2
# ╟─45bdd766-9131-11eb-31dc-f9ee1ced0db9
# ╠═a1af2f43-a6f4-4b5d-a882-abf91f109a44
# ╠═0eb3c064-9316-11eb-3c97-afba4053eb94
# ╠═322216ee-93cb-11eb-272c-b7f5e5d87ee5
# ╠═8e63260a-0d0b-4551-96ee-88322e6e389f
# ╠═cadd0620-9315-11eb-228e-135bd066038d
# ╠═660c8c2d-6b3c-4a12-8e44-ea9dd1c461ff
# ╠═9cc0ece2-f2b9-4b54-9603-75e78d564b76
# ╠═d4e70ec0-b6b9-47a0-a575-6d15fc2d25e1
# ╠═871a129f-d016-4088-b1e6-c0a7d74e0393
# ╠═e7d591ac-fd39-45e6-8c61-e5c74cdfd597
# ╠═5a401dd6-05a4-47c1-af10-d67a4d9b40dc
# ╠═d2bc643b-c247-4f55-ad67-1992bf25a30f
# ╠═82b9a9fc-feaa-4b67-98b8-577c8bdf3378
# ╠═8ad489e2-9318-11eb-1531-7515b778e6b8
# ╠═e47a9702-710b-44f4-9adf-136e04b3b04e
# ╠═2cb6d5aa-3bdd-4ce6-be0a-dcaac538b86f
# ╠═772e3bc5-a169-4c1d-aa22-2ea0720ef707
# ╠═a1bfbae2-b284-417d-953f-bc5a2493d56c
# ╠═8cf9ece5-c345-4e65-8426-2b71f9f7fea8
# ╠═ff5cdffb-0ce9-4f1b-a0ad-ce74e9d8f5e9
# ╠═627d372c-3b7a-45c0-a416-b8ffe4e5a02b
# ╠═e6565f94-4419-4359-941d-bbddd56381fc
# ╠═0b07b62b-62d1-43da-949e-77a02d74856c
# ╠═ea30ed3d-a19c-4d3f-8806-737df8954a15
# ╠═4a57223f-d47b-4d3c-b07c-5121585d20d5
# ╠═4eee17f1-117a-4fb9-88b8-5abd5bc6c69a
# ╠═bd9421ed-8fe5-45dd-9f9e-b9cf99077694
# ╠═1f0ec58c-3184-4ceb-b019-28d9bffbc16f
# ╠═c9671526-27fe-41fd-aaf4-a63a5562a726
# ╠═1481eac4-8e0a-42db-822b-399af881abb4
# ╠═c33899c4-ab79-4bc2-9511-c99900a0d2f8
# ╠═dbc8668e-26f4-43e6-9474-29d04108a7ef
# ╠═c4d19122-f3cc-4129-ac00-95fe0365d097
# ╠═cb03c3d4-c7f8-4247-9a66-99d3deda9fbb
# ╠═5361bf8e-9318-11eb-1acd-e9a10c2b30bf
# ╠═971f16da-62af-49b0-bb38-5f4c39b2b886
# ╠═4d5bdc05-0f5d-403f-b303-93e95ec408ca
# ╠═7cd1e229-96d0-4786-bbe4-d63b524e6619
# ╠═bd811db1-28b8-4869-9b8d-d10e5fa0324e
# ╠═7dc9afbf-074a-4149-bbae-7b182916e7f2
# ╠═6a0ba3b4-ed64-4c21-bbcd-ad88685b6453
# ╠═9de1cf4b-a383-4bfe-a45b-351b94c75385
# ╠═b6a12dc8-a173-4c57-a30e-5d8cbd7e74ec
# ╠═a2cc0f38-f774-489b-98a0-9c713e99b11a
# ╠═ebc85de9-e9f1-4be0-8d63-1982006c047e
# ╠═2ee3952b-99d4-4e7f-aa92-e32c16aa008a
# ╠═95ce45d1-d000-46c4-9996-2aa6db2876cf
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

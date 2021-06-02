### A Pluto.jl notebook ###
# v0.14.4

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
Pkg.add.(["Images", "ImageIO", "ImageView", "Peaks","Formatting", "StrLiterals", "StrFormat", "Format"])

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
	using PhysicalConstants.CODATA2018
	using StatsPlots
	using QuadGK
	using Peaks
	using Formatting
	using Printf
	using StrLiterals
	using StrFormat
	using Format
end

# ╔═╡ fc8b79a2-8728-11eb-2da7-e3ffa3ceef08
using DrWatson

# ╔═╡ 79cfd2fc-9046-11eb-2b13-1b877d57645d
md"# FBI\_325G12\_linTest

- Tests of linearity G1 and G2 at 325 nm"

# ╔═╡ 7f6b6c13-60b2-462d-9db5-d1421212e94e
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

# ╔═╡ 0f2f4c78-8729-11eb-2bab-27812ce8c47e
@quickactivate "LabFBI"

# ╔═╡ 5ee27d52-86fd-11eb-365e-9f2b1a095575
;pwd()

# ╔═╡ 621ec96c-86fd-11eb-1c41-379cc17180dc
datadir()

# ╔═╡ 9913df39-7aa4-42a3-8cb9-f66707185e22
srcdir("fbi.jl")

# ╔═╡ 20a461d5-16ac-4707-adec-ebfce6fc89c6
fbi = ingredients(srcdir("fbi.jl"))

# ╔═╡ cf8a41e9-ce11-4543-b67c-175ac1add17f
lfbi = ingredients(srcdir("labFbi.jl"))

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
    A, N, mol, mmol, V, L, mL, μL, M

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

# ╔═╡ f28f308a-4f5c-4dfe-8831-94de433c8fb9
md"## Discussion"

# ╔═╡ a1af2f43-a6f4-4b5d-a882-abf91f109a44
md"## Notebook functions"

# ╔═╡ 4128badd-1aa4-4997-939f-df2320f04b4e
function plot_fbi_cols(fbidf, cols, fbis, wl, colors, title, legend)
	L =collect(wl)
	FBI = fbis[1].(L)
	p = plot(fbidf.λ, fbidf[!,cols[1]],
              	   colour = colors[1],
                   shape  = :circle,
               	   label  = cols[1],
               	   legend = legend,
				   fmt = :png)
	
	lfbi = string(cols[1],"-f")
	
	p = plot!(p, L, FBI,
    	      colour = colors[1],
    	      label  = lfbi,
    	      legend = legend,
	   	      fmt = :png)
	
	for (i, col) in enumerate(cols[2:end])
		FBI = fbis[i+1].(L)
		p = plot!(p, fbidf.λ, fbidf[!,col],
               		  colour = colors[i+1],
               		  shape  = :circle,
               		  label  = cols[i+1],
               		  legend = legend,
				      fmt = :png)

		lfbi = string(cols[i+1],"-f")
		
		p = plot!(p, L, FBI,
    	       	  colour = colors[i+1],
    	          label  = lfbi,
    	          legend = legend,
	   	          fmt = :png)
	end
	
	xlabel!("λ (nm)")
	ylabel!("Intensity(AU)")
	title!(title)
	return p

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

# ╔═╡ 4af1e07a-13e4-450f-85b6-cedc5cdeda52
function plot_fbi325(fbidf, fbig, fbiint, fbibaint, wl,
		          labelFbi, labelFbiBa, colFbi, colFbiBa, title, legend)
	L =collect(wl)
	FBI = fbiint.(L)
	FBIBA = fbibaint.(L);
	fbilbl = string("FBI", fbig)
	fbibalbl = string("FBIBa", fbig)
	p1 = plot(fbidf.λ, fbidf[!,fbilbl],
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

	p3 = plot!(p2, fbidf.λ, fbidf[!,fbibalbl],
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
function plot_eps_band(fbipdf, fbibapdf, wl, title)
	L =collect(wl)
	
	EPS = [eps_band(fbipdf, fbibapdf, wl[1], λ) for λ in wl[2]:wl[end]]
	             
	p1 = plot(L[2:end], EPS,
    		  colour = :blue, legend=false,
		      fmt = :png)

	xlabel!("λ (nm)")
	ylabel!("ϵ")
	title!(title)
	return p1
end

# ╔═╡ cadd0620-9315-11eb-228e-135bd066038d
md"## Molecule spectra"

# ╔═╡ 660c8c2d-6b3c-4a12-8e44-ea9dd1c461ff
begin
fbidir         = "fluorimeter"
fbi325g2siName = "FBIG2_325nm_si_ltest.csv"
fbi325g1siName = "FBIG1_325nm_si_ltest.csv"
fbiba325g1siName = "FBIGBa1_325nm_si_ltest.csv"
end

# ╔═╡ 14202b4e-3ea0-4d36-b331-5e993bf1ec4a
fbi325g1sidf   = load_fbidf(fbidir, fbi325g1siName);

# ╔═╡ 8a4abf1a-543a-4e3c-92ce-de37d6646931
fbiba325g1sidf   = load_fbidf(fbidir, fbiba325g1siName);

# ╔═╡ 993eaa5a-91c5-4d7d-97d4-67993cb3d990
fbi325g2sidf   = load_fbidf(fbidir, fbi325g2siName);

# ╔═╡ 82b9a9fc-feaa-4b67-98b8-577c8bdf3378
md" ## Parameterized functions"

# ╔═╡ 8ad489e2-9318-11eb-1531-7515b778e6b8
begin
    wf = 375.
    we = 801.
    ws = 2.
    wl = wf:ws:we
	Ll = collect(wl)
end

# ╔═╡ 17b14af4-bf8f-4ead-a151-ca79decf40e0
colg1df = ["FBIG1",	"FBIG1_rep", "FBIG1x10"];

# ╔═╡ 6a02a440-4166-46aa-a1fc-6cc6d6c3c511
colbag1df = ["FBIG1_Ba",	"FBIG1_Ba_rep", "FBIG1x10_Ba"];

# ╔═╡ b49daeb8-c7d7-44b9-a6e8-aa18d30ffdb8
colg2df = ["FBIG2",	"FBIG2_rep", "FBIG2dil"];

# ╔═╡ 5be376d2-4464-45a1-b749-4c1b3c9ac202
md"- ##### fbi, silica, λ = 325 nm, G1.  
- FBIG1 is the reference with a concentration of 2.27E-6
- FBIG1_rep repeats after 24 h with a concentration of 2.27E-6
- FBIG1x10 uses a concentration of 2.27E-5
"

# ╔═╡ b5d306e9-78ae-4749-87a7-940b090cc47c
fbig1, fbig1rep, fbig1x10 = fbi.fbigen(wl, fbi325g1sidf,colg1df);

# ╔═╡ 53fb0f06-1389-4296-a59a-aa3e2d0f3bbe
@test fbi.test_fbis(wl, [fbig1, fbig1rep, fbig1x10])

# ╔═╡ 53b61f8d-1521-445b-96ed-fe778d94a505
@test fbi.test_fbi_pdfs(wl, [fbig1, fbig1rep, fbig1x10])

# ╔═╡ 7a47ec2b-0bcc-430b-a103-661220dc9df0
md"- ##### FbiBa, silica, λ = 325 nm, G1.  
- FBIG1 is the reference with a concentration of 7.4E-9
- FBIG1_rep repeats after 24 h with a concentration of 7.4E-9
- FBIG1x10 uses a concentration of 7.4E-8
"

# ╔═╡ 9fdeef5b-6d76-4c6d-983b-e1bcd9c51a5a
fbibag1, fbibag1rep, fbibag1x10 = fbi.fbigen(wl, fbiba325g1sidf,colbag1df);

# ╔═╡ 94113c9d-1af2-482f-a7b1-a0a659f186a0
@test fbi.test_fbis(wl, [fbibag1, fbibag1rep, fbibag1x10])

# ╔═╡ 2f8f134e-4dab-4a8c-b235-2bf22b090e81
@test fbi.test_fbi_pdfs(wl, [fbibag1, fbibag1rep, fbibag1x10])

# ╔═╡ 7c6a3a57-017e-4c9a-8bdc-57ad237225c1
md"- ##### fbi, silica, λ = 325 nm, G2.  
- FBIG2 is the reference with a concentration of 2.27E-6
- FBIG2_rep repeats after 24 h with a concentration of 2.27E-6
- FBIG1_dil uses a concentration of 2.27E-7
"

# ╔═╡ d9a66e9d-c1ec-4e6a-b1c8-4e423f0f1efe
fbig2, fbig2rep, fbig2dil = fbi.fbigen(wl, fbi325g2sidf,colg2df);

# ╔═╡ 93801b5e-23bc-4f0f-95f6-ba6c27776f03
@test fbi.test_fbis(wl, [fbig2, fbig2rep, fbig2dil])

# ╔═╡ 95e5907b-45ea-4b5c-8580-f4a0dbb6596b
@test fbi.test_fbi_pdfs(wl, [fbig2, fbig2rep, fbig2dil])

# ╔═╡ aab20e3b-2a5e-46a4-9cf8-0c8ec5576639
md"### FBIG1"

# ╔═╡ f31fa1f0-2ae4-405e-b4d9-829c7b0ac80c
pfbig1 = plot_fbi_cols(fbi325g1sidf, colg1df, [fbig1.f, fbig1rep.f, fbig1x10.f], wl, 
	[:blue, :green, :red], "FBI G1 in Silica", :topright)

# ╔═╡ ded6eb6d-22b3-4b71-bd30-cae9e80d3a3a
pfbig2 = plot_fbi_cols(fbi325g2sidf, colg2df, [fbig2.f, fbig2rep.f, fbig2dil.f], wl, 
	[:blue, :green, :red], "FBI G2 in Silica", :topright)

# ╔═╡ 6a0ba3b4-ed64-4c21-bbcd-ad88685b6453
plot(pfbig1, pfbig2, layout = (1, 2), legend = true, fmt = :png)

# ╔═╡ 996dd414-4f3a-4acb-9f85-404fba1836c5
md"### FbiBa G1 in Silica"

# ╔═╡ 7eb00152-088a-452d-9e94-33bdc6a3e6c0
pfbibag1 = plot_fbi_cols(fbiba325g1sidf,colbag1df, [fbibag1.f, fbibag1rep.f, fbibag1x10.f], wl, 
	[:blue, :green, :red], "FBIBa G1 in Silica", :topright)

# ╔═╡ 724283ba-97d3-4e56-b607-333533aa2e93
fbibag1.N /fbibag1rep.N

# ╔═╡ 37d6a78a-754a-46fc-af88-0a495d928781
fbibag1x10.N / fbibag1.N

# ╔═╡ 250a16d9-4c28-45fb-9a7d-3b3348d95bea
fbibag1x10.N / fbibag1rep.N

# ╔═╡ 890cebbe-025a-4c2b-a223-1982ac040f7c
plot_fbin(fbibag1.pdf, fbibag1x10.pdf, wl, 
		           "FbiBa ref", "FbiBa x 10", :blue, :red, "FBI-BA G1 lin", :topright)

# ╔═╡ 9f06ce2d-f0b1-4ab6-931d-700d0176742a
md"### A surprising result ?

The observed effect is that the intensity increases and the free molecules shifts to blue (or yellow) **when the concentration decreases**. Apart of that, G1 and G2 behave different with time. G1 decreases a bit its luminosity while G2 increases it. 

Can we explain this? It appears that increasing the concentration (more molecules per unit surface) may result in collective effects (molecules interfering with each other but also preventing the other molecules to be adsorbed in the surface) both of which go in the same direction. The free molecule shines less and its spectrum shifts to green (or red). **But in a monolayer the number of molecules per unit surface will be high** 

For the chelated molecule we observe that the shape changes a bit (the double pleak is lost, but this feature may be recovered after letting the molecule rest) but there is no shift towards the blue and the ratio of concentration is of the order of the ratio of recorded intensity (within a factor 2!). **But the concentration of the chelated molecule is much smaller (3 orders of magnitude) than that of the free molecule.** Thus, we could be in a 'free molecule' regime
"

# ╔═╡ 2970a155-e267-41a8-a1c8-b6b6c869d7fb
md" ### How many molecules for the different concentrations?"

# ╔═╡ 1c583a52-a335-4fe7-aa9a-fecc250db99f
md"#### FBI and FBIBa2+ in solution: Reference concentration"

# ╔═╡ 4bf66c91-0541-422b-ac7a-a6461b84d334
solfbi = lfbi.Solution("FBI-Standard-Solution", 5e-5M);

# ╔═╡ f1b0beb9-2502-4dc8-a1a6-6dfa7a900f34
solfbiba = lfbi.Solution("FBIBa-Standard-Solution", 5e-5mol/L);

# ╔═╡ 0cb2c5f3-f128-4e05-a534-c9f4d6d1be86
@test solfbi.c == solfbiba.c

# ╔═╡ 9736476c-c317-43dc-b7c9-71d62155b10d
@test lfbi.nof_molecules_volume(solfbi, 1μL) ≈ lfbi.nof_molecules_volume(5e-5M, 1μL)

# ╔═╡ ee0f3796-1474-4f10-be97-f64c8237cdf8
begin
	C=exp10.(range(-7,stop=-3,length=100)) * M 
	NOF = lfbi.nof_molecules_volume.(C, (1μL,))
	plot(C/M, NOF,
		 xaxis=:log, yaxis=:log,
    	 colour = :blue,
		 lw = 2,
		 legend=false,
	   	 fmt = :png)
	xlabel!(" Concentration (M)")
	ylabel!("nof molecules in a μL")
	title!("number of molecules vs concentration")

end

# ╔═╡ 975542df-8d55-405c-adfb-3a310611e374
md"#### FBI and FBIBa2+ in silica: Reference concentration"

# ╔═╡ 410fc4ad-ea57-4372-86ed-deabb2ffb826
sifbi = lfbi.Powder("FBI-Standard", 2.27e-5mmol/mg, 50mg/cm^2)

# ╔═╡ fbde913d-a3c4-46fe-b834-3d44c5737e37
sifbiba = lfbi.Powder("FBIBa-Standard", 7.38e-8mmol/mg, 50mg/cm^2)

# ╔═╡ 10b2ad4e-c0e6-4655-9f7c-cbdf84744e85
nfbi = f"\%7.2g(lfbi.nof_molecules_area(sifbi, 1μm^2))"

# ╔═╡ 69ff3a77-40c1-4eec-bfa8-4333d50d64a7
nfbiba = f"\%7.2g(lfbi.nof_molecules_area(sifbiba, 1μm^2))"

# ╔═╡ 8fee2554-8c20-41ea-8dac-1656ccdb24f7
md" Notice that the reference concentration, assuming that silica is packed with a density of 50 mg/cm2 is high. In a monolayer with a density of 1 mol/nm^2 one expects 10^6 mol/μ^2, while the reference concentration for FBI is $nfbi mol/μ^2 and the reference concentration for the chelated FBI is $nfbiba mol/μ^2 "

# ╔═╡ Cell order:
# ╟─79cfd2fc-9046-11eb-2b13-1b877d57645d
# ╟─7f6b6c13-60b2-462d-9db5-d1421212e94e
# ╠═4f35dcfc-8644-11eb-1c86-75608b38bf5f
# ╠═2e3c7fda-904f-11eb-2988-25604b1caad0
# ╠═4ea80a97-44a6-4205-ae30-eaf239309313
# ╠═3207f446-8643-11eb-37ba-c9aec47fcb8f
# ╠═5115917a-8644-11eb-19fc-0528741ca75d
# ╠═fc8b79a2-8728-11eb-2da7-e3ffa3ceef08
# ╠═0f2f4c78-8729-11eb-2bab-27812ce8c47e
# ╠═5ee27d52-86fd-11eb-365e-9f2b1a095575
# ╠═621ec96c-86fd-11eb-1c41-379cc17180dc
# ╠═9913df39-7aa4-42a3-8cb9-f66707185e22
# ╠═20a461d5-16ac-4707-adec-ebfce6fc89c6
# ╠═cf8a41e9-ce11-4543-b67c-175ac1add17f
# ╠═0b1cedc7-8ada-45a5-adef-fbae794dee3e
# ╟─3c01d7c8-8645-11eb-372d-4d66ae46ae72
# ╠═46b5465a-8645-11eb-0291-612455795518
# ╟─ee7f63de-904f-11eb-0e8f-1f8ec9d50811
# ╟─631833b8-0649-4f45-b1dd-e4a53091d321
# ╠═f8d3f48a-904f-11eb-2659-51f7508b434d
# ╠═9fdae9e2-6d2e-4901-8f48-14764b6800c2
# ╟─45bdd766-9131-11eb-31dc-f9ee1ced0db9
# ╟─f28f308a-4f5c-4dfe-8831-94de433c8fb9
# ╟─a1af2f43-a6f4-4b5d-a882-abf91f109a44
# ╠═4128badd-1aa4-4997-939f-df2320f04b4e
# ╟─0eb3c064-9316-11eb-3c97-afba4053eb94
# ╟─4af1e07a-13e4-450f-85b6-cedc5cdeda52
# ╟─322216ee-93cb-11eb-272c-b7f5e5d87ee5
# ╟─8e63260a-0d0b-4551-96ee-88322e6e389f
# ╟─cadd0620-9315-11eb-228e-135bd066038d
# ╠═660c8c2d-6b3c-4a12-8e44-ea9dd1c461ff
# ╠═14202b4e-3ea0-4d36-b331-5e993bf1ec4a
# ╠═8a4abf1a-543a-4e3c-92ce-de37d6646931
# ╠═993eaa5a-91c5-4d7d-97d4-67993cb3d990
# ╟─82b9a9fc-feaa-4b67-98b8-577c8bdf3378
# ╠═8ad489e2-9318-11eb-1531-7515b778e6b8
# ╠═17b14af4-bf8f-4ead-a151-ca79decf40e0
# ╠═6a02a440-4166-46aa-a1fc-6cc6d6c3c511
# ╠═b49daeb8-c7d7-44b9-a6e8-aa18d30ffdb8
# ╟─5be376d2-4464-45a1-b749-4c1b3c9ac202
# ╠═b5d306e9-78ae-4749-87a7-940b090cc47c
# ╠═53fb0f06-1389-4296-a59a-aa3e2d0f3bbe
# ╠═53b61f8d-1521-445b-96ed-fe778d94a505
# ╟─7a47ec2b-0bcc-430b-a103-661220dc9df0
# ╠═9fdeef5b-6d76-4c6d-983b-e1bcd9c51a5a
# ╠═94113c9d-1af2-482f-a7b1-a0a659f186a0
# ╠═2f8f134e-4dab-4a8c-b235-2bf22b090e81
# ╟─7c6a3a57-017e-4c9a-8bdc-57ad237225c1
# ╠═d9a66e9d-c1ec-4e6a-b1c8-4e423f0f1efe
# ╠═93801b5e-23bc-4f0f-95f6-ba6c27776f03
# ╠═95e5907b-45ea-4b5c-8580-f4a0dbb6596b
# ╟─aab20e3b-2a5e-46a4-9cf8-0c8ec5576639
# ╠═f31fa1f0-2ae4-405e-b4d9-829c7b0ac80c
# ╠═ded6eb6d-22b3-4b71-bd30-cae9e80d3a3a
# ╠═6a0ba3b4-ed64-4c21-bbcd-ad88685b6453
# ╟─996dd414-4f3a-4acb-9f85-404fba1836c5
# ╠═7eb00152-088a-452d-9e94-33bdc6a3e6c0
# ╠═724283ba-97d3-4e56-b607-333533aa2e93
# ╠═37d6a78a-754a-46fc-af88-0a495d928781
# ╠═250a16d9-4c28-45fb-9a7d-3b3348d95bea
# ╠═890cebbe-025a-4c2b-a223-1982ac040f7c
# ╟─9f06ce2d-f0b1-4ab6-931d-700d0176742a
# ╟─2970a155-e267-41a8-a1c8-b6b6c869d7fb
# ╟─1c583a52-a335-4fe7-aa9a-fecc250db99f
# ╠═4bf66c91-0541-422b-ac7a-a6461b84d334
# ╠═f1b0beb9-2502-4dc8-a1a6-6dfa7a900f34
# ╠═0cb2c5f3-f128-4e05-a534-c9f4d6d1be86
# ╠═9736476c-c317-43dc-b7c9-71d62155b10d
# ╠═ee0f3796-1474-4f10-be97-f64c8237cdf8
# ╟─975542df-8d55-405c-adfb-3a310611e374
# ╠═410fc4ad-ea57-4372-86ed-deabb2ffb826
# ╠═fbde913d-a3c4-46fe-b834-3d44c5737e37
# ╠═10b2ad4e-c0e6-4655-9f7c-cbdf84744e85
# ╠═69ff3a77-40c1-4eec-bfa8-4333d50d64a7
# ╟─8fee2554-8c20-41ea-8dac-1656ccdb24f7

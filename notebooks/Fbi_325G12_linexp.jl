### A Pluto.jl notebook ###
# v0.14.1

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
md"# FBI\_325G12\_linexp

- Linearity as a function of concentration for G1 and G2 FBI(Ba) molecules in silica and solution. 

- All studies at 325 nm
"

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

# ╔═╡ 63c6d4b2-6643-4901-a9f3-9cd1b8ef3c38
gr(size=(700,700), xtickfontsize=12, ytickfontsize=12, xguidefontsize=14, yguidefontsize=14, legendfontsize=10, dpi=100, grid=(:y, :gray, :solid, 2, 0.4));

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

# ╔═╡ f4f2ce5a-5ab7-4a9f-8c64-5db8f37a672f
md"- Plots the columns defined by cols from dataframe df"

# ╔═╡ 708d8cc6-40e2-4d5d-8039-7502ccb32d45
function plot_df(df, cols, norm, colors, xlabel, ylabel, title, legend)
	function csum(dfc)
		if norm == true
			xsum = sum(dfc)
		else
			xsum = 1.0
		end
		return xsum
	end
	
	p = plot(df.λ, df[!,cols[1]]/csum(df[!,cols[1]]),
              	   colour = colors[1],
                   #shape  = :circle,
               	   label  = cols[1],
               	   legend = legend,
				   fmt = :png)


	for (i, col) in enumerate(cols[2:end])
		p = plot!(p, df.λ, df[!,col]/csum(df[!,col]),
               		  colour = colors[i+1],
               		  #shape  = :circle,
               		  label  = cols[i+1],
               		  legend = legend,
				      fmt = :png)

	end

	xlabel!(xlabel)
	ylabel!(ylabel)
	title!(title)
	return p

end

# ╔═╡ d8c4f687-b853-4780-ba57-b09b0cacd46f
md"- Plot intensity vs concentration"

# ╔═╡ a1532ead-5044-4234-8cde-45408aa2a95a
function plot_df_ivsc(df, cols, title)
	function csum(dfc)
		if norm == true
			xsum = sum(dfc)
		else
			xsum = 1.0
		end
		return xsum
	end
	
	function ccfromcol(cols)
		CC = []
		for col in cols
			xcol = split(col,"_")
			push!(CC, parse(Float64, xcol[2]))
		end
		return CC
	end
	
	areas = [sum(df[!,col]) for col in cols]
	
	areas = areas./maximum(areas)
	
	ccs = ccfromcol(cols)
	
	p = plot(ccs, areas,
             shape  = :circle,
		     legend = :false,
			 xaxis=:log,
			 fmt = :png)


	xlabel!("Concentration (M)")
	ylabel!("Intensity (au)")
	title!(title)
	return p

end

# ╔═╡ 2462feb3-664a-4270-8674-17675e04e20f
md"- Takes a list of dfs and a list of cols and plots them. One column per df"

# ╔═╡ 9b52cdd3-8e2c-4f5d-9511-a559f64eb91e
function plot_dfs(dfs, dflbls, cols, norm, colors, xlabel, ylabel, title, legend)
	function csum(dfc)
		if norm == true
			xsum = sum(dfc)
		else
			xsum = 1.0
		end
		return xsum
	end
	lbl = string(dflbls[1],"-",cols[1])
	
	p = plot(dfs[1].λ, dfs[1][!,cols[1]]/csum(dfs[1][!,cols[1]]),
              	   colour = colors[1],
			       lw     = 2, 
                   #shape  = :circle,
               	   label  = lbl,
               	   legend = legend,
				   fmt = :png)


	for (i, col) in enumerate(cols[2:end])
		lbl = string(dflbls[i+1],"-",cols[i+1])
		p = plot!(p, dfs[i+1].λ, dfs[i+1][!,col]/csum(dfs[i+1][!,col]),
               		  colour = colors[i+1],
               		  #shape  = :circle,
				      lw = 2, 
               		  label  = lbl,
               		  legend = legend,
				      fmt = :png)

	end

	xlabel!(xlabel)
	ylabel!(ylabel)
	title!(title)
	return p

end

# ╔═╡ 4789c2b9-c23e-4833-bfd7-d3b5df7545d5
md" - Plots selected cols from two dataframes"

# ╔═╡ ec9eee46-ec8b-4e97-ad0e-c4a4d46ec8e7
function plot_df12(df1, df2, cols1, cols2, norm, colors1, colors2, style1, style2,
		           xlabel, ylabel, title, legend)

	function csum(dfc)
		if norm == true
			xsum = sum(dfc)
		else
			xsum = 1.0
		end
		return xsum
	end


	
	p = plot(df1.λ, df1[!,cols1[1]]/csum(df1[!,cols1[1]]),
             colour = colors1[1],
			 lw     = 2, 
		     style  = style1[1],
             label  = cols1[1],
             legend = legend,
		     fmt = :png)

	for (i, col) in enumerate(cols1[2:end])

			p = plot!(p, df1.λ, df1[!,col]/csum(df1[!,col]),
					  colour = colors1[i+1],
					  lw = 2, 
					  style  = style1[i+1],
					  label  = cols1[i+1],
					  legend = legend,
					  fmt = :png)

	end

	for (i, col) in enumerate(cols2[1:end])
		
		p = plot!(p, df2.λ, df2[!,col]/csum(df2[!,col]),
               	  colour = colors2[i],
				  lw = 2, 
				  style  = style2[i],
               	  label  = cols2[i],
               	  legend = legend,
				  fmt = :png)

	end

	xlabel!(xlabel)
	ylabel!(ylabel)
	title!(title)
	return p

end

# ╔═╡ bedccecd-ca4d-4597-95d8-e7f4877803c8
md"- Plot cols of a DF"

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

# ╔═╡ cadd0620-9315-11eb-228e-135bd066038d
md"## Molecule spectra"

# ╔═╡ 660c8c2d-6b3c-4a12-8e44-ea9dd1c461ff
begin
	fbidir         = "fluorimeter/linexp"
	fbig1solName   = "FBIG1_disol_linexp.csv"
	fbig2solName   = "FBIG2_disol_linexp.csv"
	fbig1siName    = "FBIG1_silica_linexp.csv"
	fbig2siName    = "FBIG2_silica_linexp.csv"
	fbibag1siName  = "FBIBaG1_silica_linexp.csv"
	fbibag2siName  = "FBIBaG2_silica_linexp.csv"
	
end

# ╔═╡ 14202b4e-3ea0-4d36-b331-5e993bf1ec4a
fbibag1sidf   = fbi.load_fbidf(fbidir, fbibag1siName);

# ╔═╡ a0a07dcd-dda3-4b87-bc12-178d8f3c47ab
fbibag2sidf   = fbi.load_fbidf(fbidir, fbibag2siName);

# ╔═╡ 715465c2-1415-4f42-abb3-20e3edf821f3
fbig1soldf    = fbi.load_fbidf(fbidir, fbig1solName);

# ╔═╡ 233199b7-8988-442e-9960-6b94270f548a
fbig2soldf    = fbi.load_fbidf(fbidir, fbig2solName);

# ╔═╡ 15f2fc46-acd9-46d8-8fc3-97342cf1f552
fbig1sidf     = fbi.load_fbidf(fbidir, fbig1siName);

# ╔═╡ 66176c3c-c127-4132-a2cc-ea6f864c365a
fbig2sidf     = fbi.load_fbidf(fbidir, fbig2siName);

# ╔═╡ 8ce1b48d-8ef6-4cf0-8460-bb366e657fb4
md"- Data Frames columns"

# ╔═╡ 17b14af4-bf8f-4ead-a151-ca79decf40e0
colg1soldf = ["FBIG1_5E-6","FBIG1_5E-5", "FBIG1_1E-4",	"FBIG1_5E-4"];

# ╔═╡ 9afb1a36-a5d2-4718-b926-55c77bb3b26d
colg1sidf = ["FBIG1_2.3E-5", "FBIG1_2.3E-6", "FBIG1_2.3E-7", "FBIG1_7.4E-8",	"FBIG1_2.3E-8"];

# ╔═╡ c485625b-c653-4b37-be27-598ba3c8010c
colbag1sidf = ["FBIBaG1_7.4E-7", "FBIBaG1_7.4E-8", "FBIBaG1_7.4E-09",	"FBIBaG1_7.4E-10"];

# ╔═╡ a4bd50fb-0268-4564-b9a1-9721ba592e60
colg2soldf = ["FBIG2_5E-6",	"FBIG2_5E-5", "FBIG2_1E-4",	"FBIG2_5E-4"];

# ╔═╡ ff0d10c5-89e2-4cbc-81d9-7276986168dd
colg2sidf = ["FBIG2_2.3E-6", "FBIG2_2.3E-7", "FBIG2_2.3E-8"];

# ╔═╡ ab7210ff-8365-488a-b127-d4ad9e0e2c3c
colbag2sidf = ["FBIBaG2_7.4E-8",	"FBIBaG2_7.4E-9",	"FBIBaG2_7.4E-10"];

# ╔═╡ 2bc8d78f-1412-45e3-9c96-04e3683e3c39
md"#### Dependence of FBI G1 with concentration 
- Observe saturation (with huge decrease of light emission) at 5E-4 
- Saturation also occurs in Silica at around 2.3E-5
- The chelated molecule behaves linearly, except at very low concentration (but the concentrations are much smaller than those for FBI
- The free molecule moves towards the blue, developing 'shoulders' clearly visible in solution and even more in silica, **as the concentration descreases**. Conversely, the molecule spectrum becomes smooth and shifts towards green as the concentation increases.
- The chelated molecule also shifts towards blue as the concentration decreases. 
"

# ╔═╡ 9820feab-3598-470a-a6dd-5ba581e84716
pfbig1soldf = plot_df(fbig1soldf, colg1soldf, false, markercolors, "λ(nm)", "I (au)", "G1 Solution", :topright);

# ╔═╡ 5aa22150-0f63-44c7-8a9c-7cbdf305a11a
pafbig1sol =  plot_df_ivsc(fbig1soldf, colg1soldf, "I vs C FBI-G1 SOL");

# ╔═╡ 70f940ca-5532-4432-83c8-eabd0b3af8a6
pfbig1sidf = plot_df(fbig1sidf, colg1sidf, false, markercolors, "λ(nm)", "I (au)", "G1 Silica", :topright);

# ╔═╡ ada90fd8-a31e-4750-ae41-fa5f47b0377e
pafbig1si =  plot_df_ivsc(fbig1sidf, colg1sidf, "I vs C FBI-G1 Si");

# ╔═╡ 98185a08-d81c-48c0-9fbb-ea736b4106a6
pfbibag1sidf = plot_df(fbibag1sidf, colbag1sidf, false, markercolors, "λ(nm)", "I (au)", "G1 FbiBa Silica", :topright);

# ╔═╡ 095db23e-fe28-455f-ba0b-82ed4d8158ad
pafbibag1si=  plot_df_ivsc(fbibag1sidf, colbag1sidf, "I vs C FBIBa-G1 Si");

# ╔═╡ 7cbc12e3-3973-4b50-a18f-57bb83517d7c
pnfbig1soldf = plot_df(fbig1soldf, colg1soldf, true, markercolors, "λ(nm)", "I (au)", "G1 Solution PDF", :topright);

# ╔═╡ 4e05e0aa-823a-4365-a348-fda7f7d2e4fd
pnfbig1sidf = plot_df(fbig1sidf, colg1sidf, true, markercolors, "λ(nm)", "I (au)", "G1 Silica PDF", :topright);

# ╔═╡ 3bc0b72a-fe70-4bb9-8266-7ef45cf8e7a9
pnfbibag1sidf = plot_df(fbibag1sidf, colbag1sidf, true, markercolors, "λ(nm)", "I (au)", "G1 FbiBa Silica PDF", :topright);

# ╔═╡ a1cbf361-13c3-46bb-958d-6b999c399d17
md" #### G1: solutions vs Silica"

# ╔═╡ b36e7daf-5d3a-4e68-a6e7-a70fac9125bd
plot(pfbig1soldf, pnfbig1soldf, pfbig1sidf, pnfbig1sidf, layout = (2, 2), legend = true, fmt = :png)

# ╔═╡ 47d9a9cd-bf9a-40c2-94b7-62ff2722bb40
md" #### G1: in Silica: FBi vs FbiBa"

# ╔═╡ 64772623-2d88-4b43-b1ce-f20d8bb5288c
plot(pfbig1sidf, pnfbig1sidf, pfbibag1sidf, pnfbibag1sidf, layout = (2, 2), legend = true, fmt = :png)

# ╔═╡ ca2af6eb-f830-4ab7-b378-a6fab79e751d
md"#### Dependence of FBI G2  with concentration 
- Observe saturation (with huge decrease of light emission) already at 1E-4 
- Saturation already occurs in Silica at around 1E-6
- The chelated molecule behaves quite linearly
- The free molecule moves towards the blue, in solution and even more in silica, **as the concentration descreases**, but shoulders are less pronounced. Conversely, the molecule spectrum becomes smooth and shifts towards green as the concentation increases.
- The chelated molecule does not shift too much with concentration. 
"

# ╔═╡ f114b802-5f6d-4a52-a4e1-54a4173779fc
pfbig2soldf = plot_df(fbig2soldf, colg2soldf, false, markercolors, "λ(nm)", "I (au)", "G2 Solution", :topright);

# ╔═╡ 0121ab11-145f-4985-92dc-e8f1d67cbe0c
pafbig2sol=  plot_df_ivsc(fbig2soldf, colg2soldf, "I vs C FBI-G2 SOL");

# ╔═╡ a733d819-deb9-4f2b-96cf-7af90e0cb4f5
pfbig2sidf = plot_df(fbig2sidf, colg2sidf, false, markercolors, "λ(nm)", "I (au)", 
	"G2 Silica", :topright);

# ╔═╡ 60666285-6fd5-4cff-884d-2e32a139b748
pafbig2si=  plot_df_ivsc(fbig2sidf, colg2sidf, "I vs C FBI-G2 Si");

# ╔═╡ 5de64351-5877-4c39-ad96-63b97f7b2f33
pfbibag2sidf = plot_df(fbibag2sidf, colbag2sidf, false, markercolors, "λ(nm)", "I (au)", "G2 FbiBa Silica", :topright);

# ╔═╡ a314fde3-d779-475f-8e11-6bdd87ad44e7
pafbibag2si=  plot_df_ivsc(fbibag2sidf, colbag2sidf, "I vs C FBIBa-G2 Si");

# ╔═╡ c4f72191-2c19-4423-95d2-06bcc5b8ff7f
pnfbig2soldf = plot_df(fbig2soldf, colg2soldf, true, markercolors, "λ(nm)", "I (au)", "G2 Solution PDF", :topright);

# ╔═╡ 97e7a84b-90b7-48e4-9c9d-95f4224f5c7c
pnfbig2sidf = plot_df(fbig2sidf, colg2sidf, true, markercolors, "λ(nm)", "I (au)", 
	"G2 Silica PDF", :topright);

# ╔═╡ 8eb4af33-e288-45bf-9907-05db62b28a74
pnfbibag2sidf = plot_df(fbibag2sidf, colbag2sidf, true, markercolors, "λ(nm)", "I (au)", "G2 FbiBa Silica PDF", :topright);

# ╔═╡ 8682ab56-92d3-47ec-aadd-9d3f15810b8a
pnfbi_free_ba_g1_si = plot_df12(fbig1sidf, fbibag1sidf, ["FBIG1_2.3E-5"], ["FBIBaG1_7.4E-7", "FBIBaG1_7.4E-8", "FBIBaG1_7.4E-09",	"FBIBaG1_7.4E-10"], true, [:green], [:blue, :orange, :yellow, :red], [:solid], [:dash,:dash,:dash,:dash],
		           "λ(nm)", "I (au)", "G1 Fbi and FbiBa Silica PDF", :topright);

# ╔═╡ 353692fa-1eff-4610-be20-4866659e4c91
pnfbi_free_ba_g2_si = plot_df12(fbig2sidf, fbibag2sidf, ["FBIG2_2.3E-6"], ["FBIBaG2_7.4E-8",	"FBIBaG2_7.4E-9",	"FBIBaG2_7.4E-10"], true, [:green], [:blue, :orange, :red], [:solid], [:dash,:dash,:dash],
		           "λ(nm)", "I (au)", "G2 Fbi and FbiBa Silica PDF", :topright);

# ╔═╡ ebf12f7f-3cae-4e3c-9f80-a24e714c9505
pnfbi_free_ba_g2_sol_si = plot_df12(fbig2soldf, fbibag2sidf, ["FBIG2_5E-4"], ["FBIBaG2_7.4E-8",	"FBIBaG2_7.4E-9",	"FBIBaG2_7.4E-10"], true, [:green], [:blue, :orange, :red], [:solid], [:dash,:dash,:dash],
		           "λ(nm)", "I (au)", "G2 Fbi (Sol) and FbiBa Silica PDF", :topright);

# ╔═╡ 1b284a60-bfa1-4d40-8a6c-b8e117498440
md"#### G2 in Silica"

# ╔═╡ 688ff111-4e73-428c-a95b-59e6dbb53b6b
plot(pfbig2sidf, pnfbig2sidf, pfbibag2sidf, pnfbibag2sidf, layout = (2, 2), legend = :topright,  fmt = :png)

# ╔═╡ eed3f804-3712-4c5b-aa2d-ab2f6527766a
md" #### FBI and FbiBa, G1 and G2

- There is more apparent separation for G1, but the reason is that we are comparing with FBI concentrated.

- If one uses the FBI G2 from solution (which is more concentrated) separation is much better.
"

# ╔═╡ c0133195-98e1-410f-8fe6-d001f89ad32f
plot(pnfbi_free_ba_g1_si, pnfbi_free_ba_g2_si,  layout = (1, 2), legend = true, fmt = :png)

# ╔═╡ 6ae99831-0640-400c-90d8-8131a0359251
plot(pnfbi_free_ba_g2_si, pnfbi_free_ba_g2_sol_si,  layout = (1, 2), legend = true, fmt = :png)

# ╔═╡ 9df9ee8e-ede3-4c48-8dca-3f84e3731386
md"#### Linearity (and lack of) for FBI"

# ╔═╡ eca19f3f-efbc-473e-a9ed-f43c0f53a697
plot(pafbig1sol, pafbig1si, pafbibag1si, pafbig2sol, pafbig2si, pafbibag2si,  layout = (3, 2), legend = false, fmt = :png)

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
sifbi = lfbi.Powder("FBI-Standard", 2.27e-5mmol/mg, 30mg/cm^2)

# ╔═╡ fbde913d-a3c4-46fe-b834-3d44c5737e37
sifbiba = lfbi.Powder("FBIBa-Standard", 7.38e-8mmol/mg, 30mg/cm^2)

# ╔═╡ 10b2ad4e-c0e6-4655-9f7c-cbdf84744e85
nfbi = f"\%7.2g(lfbi.nof_molecules_area(sifbi, 1μm^2))"

# ╔═╡ 69ff3a77-40c1-4eec-bfa8-4333d50d64a7
nfbiba = f"\%7.2g(lfbi.nof_molecules_area(sifbiba, 1μm^2))"

# ╔═╡ 8fee2554-8c20-41ea-8dac-1656ccdb24f7
md" Notice that the reference concentration, assuming that silica is packed with a density of 30 mg/cm2 is high. In a monolayer with a density of 1 mol/nm^2 one expects 10^6 mol/μ^2, while the reference concentration for FBI is $nfbi mol/μ^2 and the reference concentration for the chelated FBI is $nfbiba mol/μ^2 "

# ╔═╡ 67d8955d-55fe-4de2-baae-95aa9078aed5
md"## Discussion

The non-linearity of response will be likely due to self-absorption at large concentrations. At the same time, large concentration appears to shift the spectrum of FBI towards green (red), which in principle is positive.

On the other hand, the behaviour of a single FBI chelated molecule surrounded by unchelated FBI molecules is hard to predict. It could be that the molecule does not suffer self-absorption, but is difficult to quantify. 
"

# ╔═╡ Cell order:
# ╠═79cfd2fc-9046-11eb-2b13-1b877d57645d
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
# ╠═63c6d4b2-6643-4901-a9f3-9cd1b8ef3c38
# ╟─3c01d7c8-8645-11eb-372d-4d66ae46ae72
# ╠═46b5465a-8645-11eb-0291-612455795518
# ╟─ee7f63de-904f-11eb-0e8f-1f8ec9d50811
# ╟─631833b8-0649-4f45-b1dd-e4a53091d321
# ╠═f8d3f48a-904f-11eb-2659-51f7508b434d
# ╠═9fdae9e2-6d2e-4901-8f48-14764b6800c2
# ╟─45bdd766-9131-11eb-31dc-f9ee1ced0db9
# ╟─f28f308a-4f5c-4dfe-8831-94de433c8fb9
# ╟─a1af2f43-a6f4-4b5d-a882-abf91f109a44
# ╟─f4f2ce5a-5ab7-4a9f-8c64-5db8f37a672f
# ╟─708d8cc6-40e2-4d5d-8039-7502ccb32d45
# ╟─d8c4f687-b853-4780-ba57-b09b0cacd46f
# ╟─a1532ead-5044-4234-8cde-45408aa2a95a
# ╟─2462feb3-664a-4270-8674-17675e04e20f
# ╟─9b52cdd3-8e2c-4f5d-9511-a559f64eb91e
# ╟─4789c2b9-c23e-4833-bfd7-d3b5df7545d5
# ╟─ec9eee46-ec8b-4e97-ad0e-c4a4d46ec8e7
# ╟─bedccecd-ca4d-4597-95d8-e7f4877803c8
# ╟─4128badd-1aa4-4997-939f-df2320f04b4e
# ╟─cadd0620-9315-11eb-228e-135bd066038d
# ╠═660c8c2d-6b3c-4a12-8e44-ea9dd1c461ff
# ╠═14202b4e-3ea0-4d36-b331-5e993bf1ec4a
# ╠═a0a07dcd-dda3-4b87-bc12-178d8f3c47ab
# ╠═715465c2-1415-4f42-abb3-20e3edf821f3
# ╠═233199b7-8988-442e-9960-6b94270f548a
# ╠═15f2fc46-acd9-46d8-8fc3-97342cf1f552
# ╠═66176c3c-c127-4132-a2cc-ea6f864c365a
# ╟─8ce1b48d-8ef6-4cf0-8460-bb366e657fb4
# ╠═17b14af4-bf8f-4ead-a151-ca79decf40e0
# ╠═9afb1a36-a5d2-4718-b926-55c77bb3b26d
# ╠═c485625b-c653-4b37-be27-598ba3c8010c
# ╠═a4bd50fb-0268-4564-b9a1-9721ba592e60
# ╠═ff0d10c5-89e2-4cbc-81d9-7276986168dd
# ╠═ab7210ff-8365-488a-b127-d4ad9e0e2c3c
# ╠═2bc8d78f-1412-45e3-9c96-04e3683e3c39
# ╠═9820feab-3598-470a-a6dd-5ba581e84716
# ╠═5aa22150-0f63-44c7-8a9c-7cbdf305a11a
# ╠═70f940ca-5532-4432-83c8-eabd0b3af8a6
# ╠═ada90fd8-a31e-4750-ae41-fa5f47b0377e
# ╠═98185a08-d81c-48c0-9fbb-ea736b4106a6
# ╠═095db23e-fe28-455f-ba0b-82ed4d8158ad
# ╠═7cbc12e3-3973-4b50-a18f-57bb83517d7c
# ╠═4e05e0aa-823a-4365-a348-fda7f7d2e4fd
# ╠═3bc0b72a-fe70-4bb9-8266-7ef45cf8e7a9
# ╟─a1cbf361-13c3-46bb-958d-6b999c399d17
# ╠═b36e7daf-5d3a-4e68-a6e7-a70fac9125bd
# ╟─47d9a9cd-bf9a-40c2-94b7-62ff2722bb40
# ╠═64772623-2d88-4b43-b1ce-f20d8bb5288c
# ╟─ca2af6eb-f830-4ab7-b378-a6fab79e751d
# ╠═f114b802-5f6d-4a52-a4e1-54a4173779fc
# ╠═0121ab11-145f-4985-92dc-e8f1d67cbe0c
# ╠═a733d819-deb9-4f2b-96cf-7af90e0cb4f5
# ╠═60666285-6fd5-4cff-884d-2e32a139b748
# ╠═5de64351-5877-4c39-ad96-63b97f7b2f33
# ╠═a314fde3-d779-475f-8e11-6bdd87ad44e7
# ╠═c4f72191-2c19-4423-95d2-06bcc5b8ff7f
# ╠═97e7a84b-90b7-48e4-9c9d-95f4224f5c7c
# ╠═8eb4af33-e288-45bf-9907-05db62b28a74
# ╠═8682ab56-92d3-47ec-aadd-9d3f15810b8a
# ╠═353692fa-1eff-4610-be20-4866659e4c91
# ╠═ebf12f7f-3cae-4e3c-9f80-a24e714c9505
# ╟─1b284a60-bfa1-4d40-8a6c-b8e117498440
# ╠═688ff111-4e73-428c-a95b-59e6dbb53b6b
# ╟─eed3f804-3712-4c5b-aa2d-ab2f6527766a
# ╠═c0133195-98e1-410f-8fe6-d001f89ad32f
# ╠═6ae99831-0640-400c-90d8-8131a0359251
# ╟─9df9ee8e-ede3-4c48-8dca-3f84e3731386
# ╠═eca19f3f-efbc-473e-a9ed-f43c0f53a697
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
# ╠═67d8955d-55fe-4de2-baae-95aa9078aed5

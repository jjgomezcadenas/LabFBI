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
Pkg.add.(["Images", "ImageIO", "ImageView"])

# ╔═╡ 4a239876-4667-48e5-8f9e-619f2b30765a
Pkg.add("Peaks")

# ╔═╡ 0787a9f6-d81f-4862-9f3a-5974129977bd
Pkg.add("Formatting")

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

# ╔═╡ aef7b341-bf25-4024-8ec6-7f9311736c72
using Formatting

# ╔═╡ ce6cc961-96cf-450a-b029-bf31f2e0c019
using Printf

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

# ╔═╡ f28f308a-4f5c-4dfe-8831-94de433c8fb9
md"## Discussion"

# ╔═╡ 861b21ea-6fdd-4988-af50-b38f801de3a0
md"
Given an irreducible area (e.g, a pixel) with $n_c$ chelated molecules
and $n_u$ unchelated molecules iluminated by a laser beam of photon density I during a time t, we find that:

$\Gamma_c  = n_c \sigma_c Q_c \int_0^t I dt = n_c \gamma_c$

$\Gamma_u  = n_u \sigma_u Q_u \int_0^t I dt = n_u \gamma_u$

Where $\Gamma_c$ is the number of photons emitted by the chelated molecules in the time t and $\Gamma_u$ is the number of photons emitted by the unchelated molecules in the time t. By construction

$\gamma_c = \sigma_c Q_c \int_0^t I dt$
$\gamma_u = \sigma_u Q_u \int_0^t I dt$

are the number of photons per molecule emitted by the chelated and unchelated species during time t.

If the optical system has an overall geometrical acceptance of $\epsilon_g$, and a detection efficiency of $\epsilon_d$ the number of recorded photons in the corresponding detector (e.g, a pixel in the CCD) is:

$N_c =  n_c \gamma_c \epsilon$

$N_u =  n_u \gamma_u \epsilon$

Where $\epsilon = \epsilon_g \epsilon_d$.

The emission of light by the molecules follows a spectrum $\phi(\lambda)$, where
$\lambda$ is the wavelength of the emission. Then:

$\gamma_{u,c} = \int \phi_{u,c}(\lambda) d \lambda$

Since the emission of chelated and unchelated molecules are chromatically shifted, one can define a selection band where the response of the chelated signal is maximised with respect to the response of the unchelated background. The band is defined between wavelengths $\lambda_{min}$ and $\lambda_{max}$. The band efficiency is defined as:

$\epsilon^b_{u,c} = \frac{\int_{\lambda_{min}}^{\lambda_{max}} \phi_{u,c}(\lambda) d \lambda}{\int \phi_{u,c}(\lambda) d \lambda} = \frac{\Phi_{u,c}}{\gamma_{u,c}}$

where:

$\Phi_{u,c}= \int_{\lambda_{min}}^{\lambda_{max}} \phi_{u,c}(\lambda) d \lambda$


The observed signal in the selection color band is:

$S =  N_c \epsilon_c^b = n_c \gamma_c \epsilon_c^b \epsilon = n_c \Phi_{c} \epsilon$
$B =  N_u \epsilon_u^b = n_u \gamma_u \epsilon_u^b \epsilon= n_u \Phi_{u} \epsilon$


The SNR is:

$SNR = \frac{S}{\sqrt{B}} =\frac{n_s}{\sqrt{n_b}} \frac{\Phi_c}{\sqrt{\Phi_u}}\sqrt{\epsilon}$

In a fluorimeter (or laser setup) measurement in which the concentration of signal and background is known, we have:

$N_c =  C_c V \gamma_c \epsilon$

$N_u =  C_u V \gamma_u \epsilon$

Where $V$ is the volume of the sample and $C_c, C_u$ are the concentrations of the chelated  and unchelated species. We define the ratio of the total areas measured in the fluorimeter or laser as:

$r_{cu} = \frac{N_c}{N_u} = \frac{C_c}{C_u} \frac{\gamma_c}{\gamma_u} =C_{cu} \gamma_{cu}$

where $C_{cu} = \frac{C_c}{C_u}$  $\gamma_{cu} = \frac{\gamma_c}{\gamma_u}$. Thus:

$\gamma_{cu}  = \frac{r_{cu}}{C_{cu}}$

We can also define the Signal-to-Background ratio (SBR) as:

$SBR  =\frac{S}{B} = \frac{N_c}{N_u} \frac{\epsilon_c^b}{\epsilon_u^b} = \frac{\Phi_{cu}}{C_{cu}}$

Notice that SNR and SBR have to be as large a possible. In order to quantify SNR one needs to know the density of the sensor layer and the setup efficiency, as well as the fluorescence cross sections, laser intensity, etc. SBR can be quantified directly from fluorimeter of laser measurements. 
"

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
fbiname      = "EDI044_FBI_round4_em325_375_405.csv"
fbibaname    = "EDI044_FBI_Ba_round4_em325_375_405.csv"
fbiabsname   = "EDI_029_absorption.csv"
fbig2name    = "FBIG2.csv"
fbibag2name  = "FBIBaG2.csv"
	
fbiabsG2name = 	"FBIG2_absorption.csv"
fbidir       = "fluorimeter"
fbi325name   =	"FBI_G1_vs_FBI_G2_em325nm.csv"
end

# ╔═╡ 9cc0ece2-f2b9-4b54-9603-75e78d564b76
fbidf   = load_fbidf(fbidir, fbiname);

# ╔═╡ d4e70ec0-b6b9-47a0-a575-6d15fc2d25e1
fbig2df = load_fbidf(fbidir, fbig2name);

# ╔═╡ 871a129f-d016-4088-b1e6-c0a7d74e0393
fbibadf = load_fbidf(fbidir, fbibaname);

# ╔═╡ e7d591ac-fd39-45e6-8c61-e5c74cdfd597
fbibag2df = load_fbidf(fbidir, fbibag2name);

# ╔═╡ 5a401dd6-05a4-47c1-af10-d67a4d9b40dc
fbiabsdf = load_fbidf(fbidir, fbiabsname);

# ╔═╡ d2bc643b-c247-4f55-ad67-1992bf25a30f
fbiabsg2df = load_fbidf(fbidir, fbiabsG2name);

# ╔═╡ 993eaa5a-91c5-4d7d-97d4-67993cb3d990
fbi325df   = load_fbidf(fbidir, fbi325name);

# ╔═╡ 82b9a9fc-feaa-4b67-98b8-577c8bdf3378
md" ## G1"

# ╔═╡ 8ad489e2-9318-11eb-1531-7515b778e6b8
begin
    wf = 350.
    we = 800.
    ws = 2.
    wl = wf:ws:we
end

# ╔═╡ e47a9702-710b-44f4-9adf-136e04b3b04e
fbi325 = dftof(wl, fbi325df, "FBIG1");

# ╔═╡ 2cb6d5aa-3bdd-4ce6-be0a-dcaac538b86f
fbiba325 = dftof(wl, fbi325df, "FBIBaG1");

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

# ╔═╡ 0b07b62b-62d1-43da-949e-77a02d74856c
fbi325g2 = dftof(wl, fbi325df, "FBIG2");

# ╔═╡ ea30ed3d-a19c-4d3f-8806-737df8954a15
fbiba325g2 = dftof(wl, fbi325df, "FBIBaG2");

# ╔═╡ 4a57223f-d47b-4d3c-b07c-5121585d20d5
fbiE325G2N, fbig2pdf   = ftopdf(wl, fbi325g2)

# ╔═╡ 4eee17f1-117a-4fb9-88b8-5abd5bc6c69a
fbibaE325G2N, fbibag2pdf = ftopdf(wl, fbiba325g2)

# ╔═╡ bd9421ed-8fe5-45dd-9f9e-b9cf99077694
@test quadgk(fbig2pdf, wl[1], wl[end])[1]  ≈ 1

# ╔═╡ 1f0ec58c-3184-4ceb-b019-28d9bffbc16f
@test quadgk(fbibag2pdf, wl[1], wl[end])[1]  ≈ 1

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

# ╔═╡ a2c85efc-6aca-4794-844e-d9be1919c16e
pfbiE325G1 = plot_fbi325(fbi325df, "G1", fbi325, fbiba325, wl,
		          "FbiG1", "FbiBaG1", :green, :blue, "FBI/BA G1 325 nm", :topright)

# ╔═╡ 971f16da-62af-49b0-bb38-5f4c39b2b886
pfbi = plot_fbin(fbipdf, fbibapdf, wl, 
		  "Fbi G1 pdf", "FbiBa G1 pdf", :green, :blue, "PDF FBI G1 325 nm", :topright)

# ╔═╡ 4d5bdc05-0f5d-403f-b303-93e95ec408ca
md"### FBIG2"

# ╔═╡ 1ffd3dcd-87e7-4bee-8743-91d6bb35f077
pfbiE325G2 = plot_fbi325(fbi325df, "G2", fbi325g2, fbiba325g2, wl,
		          "FbiG2", "FbiBaG2", :green, :blue, "FBI/BA G2 325 nm", :topright)

# ╔═╡ 7dc9afbf-074a-4149-bbae-7b182916e7f2
pfbig2 = plot_fbin(fbig2pdf, fbibag2pdf, wl, "FBI", "FBIBa2+", :red, :orange, "FBI G2 (325 nm)", :topleft)

# ╔═╡ 6a0ba3b4-ed64-4c21-bbcd-ad88685b6453
plot(pfbi, pfbig2, layout = (1, 2), legend = false, fmt = :png)

# ╔═╡ 9de1cf4b-a383-4bfe-a45b-351b94c75385
plot(pfbiE325G1, pfbiE325G2, layout = (1, 2), legend = false, fmt = :png)

# ╔═╡ 9f06ce2d-f0b1-4ab6-931d-700d0176742a
md"### Figure of merit

The SNR (see discussion above) depends on $\frac{\Phi_c}{\sqrt{\Phi_u}}$, while the SBR is proportional to $\frac{\Phi_c}{\Phi_u}$, and

$\Phi_{u,c}= \int_{\lambda_{min}}^{\lambda_{max}} \phi_{u,c}(\lambda) d \lambda$

In a real experiment, one chooses the selection band to maximise the figure of merit $\frac{\Phi_c}{\sqrt{\Phi_u}}$ (rahter than to maximise $\Phi_{u,c}$. Below we find the maxima of the fom for G1 and G2 and then compare the resulting $\Phi_{u,c}$.
"

# ╔═╡ d99b60c2-e651-4ed5-9aec-44b3e410863b
md"#### SBR for G1"

# ╔═╡ a6f9aef3-01ba-46aa-a960-41242e14a2c8
pfbig1v1 = plot_fbin(fbi325, fbiba325, 400.0:450.0, "FBI", "FBIBa2+", :green, :blue, "FBI G1 (325 nm)", :topleft)

# ╔═╡ c995629f-4728-460b-ab64-9fa47773ce58
plot_eps_band(fbi325, fbiba325, 400.0:480.0, "Band efficiency ratio  FBI G1")

# ╔═╡ 0c88b0b8-94b1-4ce6-88f2-0a813ed6fce4
Φg1u = qpdf(fbi325, 400.0, 450.)

# ╔═╡ 34e6fa38-9a1c-4c64-9287-60f7dbc7a17c
ϵg1u = Φg1u / fbiE325G1N

# ╔═╡ e1c26455-1ba4-438b-9022-cb419dde4371
Φg1c = qpdf(fbiba325, 400.0, 450.)

# ╔═╡ 52b773e3-a5d7-4f15-81d8-3607078b72eb
Φg1uc = Φg1c / Φg1u

# ╔═╡ 097a9891-3c86-4162-a6b8-3c2ac94dc12b
fomg1 = Φg1c / sqrt(Φg1u)

# ╔═╡ 0e2a6d5f-ec49-4460-ab87-8ba293032d40
ϵg1c = Φg1c / fbibaE325G1N

# ╔═╡ bca3578d-3e45-4a31-894e-94aac3d20839
md" For G1:

i) $\Phi_{u,c}=$ $Φg1uc 

ii) $\frac{\Phi_c}{\sqrt{\Phi_u}}=$ $fomg1

iii) $\epsilon^b_c=$ $ϵg1c"

# ╔═╡ 6651e306-f2fb-49ae-8632-9481f5f4b246
md"#### SBR for G2"

# ╔═╡ b0a7893f-856f-42ec-bf8b-22a4ca2907d4
pfbig2v1 = plot_fbin(fbi325g2, fbiba325g2, 480.0:570.0, "FBI", "FBIBa2+", :red, :orange, "FBI G2 (325 nm)", :topleft)

# ╔═╡ dbd925fd-37b3-46fe-ac4b-883681dcb915
plot_eps_band(fbi325g2, fbiba325g2, 480.0:570.0, "Band efficiency ratio  FBI G2")

# ╔═╡ 87fdfcaa-ce6d-4029-b94e-8e9441399a54
Φg2u = qpdf(fbi325g2, 480.0, 550.)

# ╔═╡ 04d1d9f5-87dc-40f0-a609-c3a7e2f98b00
Φg2c = qpdf(fbiba325g2, 480.0, 550.)

# ╔═╡ fb2f1b5e-2af2-403c-8fa1-0de77669887b
Φg2uc = Φg2c / Φg2u

# ╔═╡ 88a52fc2-2cac-4bed-bb02-a77d6a678c6f
ϵg2c = Φg2c / fbibaE325G2N

# ╔═╡ 7b6bf81f-a919-4d9f-a16d-e97f2dc3af95
fomg2 = Φg2c / sqrt(Φg2u)

# ╔═╡ 956c71a3-bf0c-49cc-849a-71d5bd701ebc
md" For G2:

i) $\Phi_{u,c}=$ $Φg2uc 

ii) $\frac{\Phi_c}{\sqrt{\Phi_u}}=$ $fomg2

iii) $\epsilon^b_c=$ $ϵg2c"

# ╔═╡ 4495993d-115f-4e20-9806-8793cf40d777
ΦucR2 = Φg2uc / Φg1uc

# ╔═╡ bab49c82-ca0c-4f25-9d4a-4ca38b25054d
fomR2 = fomg2 / fomg1

# ╔═╡ ff9e4c96-1353-423f-af2b-8fc92b04b5da
md"### G1 vs G2

i) $\frac{\Phi_{u,c}^{g2}}{\Phi_{u,c}^{g1}}=$ $ΦucR2

ii) fom (g2/g1) = $fomR2
"

# ╔═╡ Cell order:
# ╟─79cfd2fc-9046-11eb-2b13-1b877d57645d
# ╠═4f35dcfc-8644-11eb-1c86-75608b38bf5f
# ╠═2e3c7fda-904f-11eb-2988-25604b1caad0
# ╠═4ea80a97-44a6-4205-ae30-eaf239309313
# ╠═4a239876-4667-48e5-8f9e-619f2b30765a
# ╠═0787a9f6-d81f-4862-9f3a-5974129977bd
# ╠═3207f446-8643-11eb-37ba-c9aec47fcb8f
# ╠═5115917a-8644-11eb-19fc-0528741ca75d
# ╠═8ba78a98-904f-11eb-1977-2750643d2f9f
# ╠═f2faa594-9316-11eb-2451-33e4b4d3a4f5
# ╠═5356699e-93ca-11eb-0afd-e7d0edc4231e
# ╠═887cfa56-2fe0-4173-83fd-0ade21c8ff0d
# ╠═aef7b341-bf25-4024-8ec6-7f9311736c72
# ╠═ce6cc961-96cf-450a-b029-bf31f2e0c019
# ╠═fc8b79a2-8728-11eb-2da7-e3ffa3ceef08
# ╠═0f2f4c78-8729-11eb-2bab-27812ce8c47e
# ╠═5ee27d52-86fd-11eb-365e-9f2b1a095575
# ╠═621ec96c-86fd-11eb-1c41-379cc17180dc
# ╠═9913df39-7aa4-42a3-8cb9-f66707185e22
# ╟─20a461d5-16ac-4707-adec-ebfce6fc89c6
# ╠═0b1cedc7-8ada-45a5-adef-fbae794dee3e
# ╟─3c01d7c8-8645-11eb-372d-4d66ae46ae72
# ╠═46b5465a-8645-11eb-0291-612455795518
# ╟─ee7f63de-904f-11eb-0e8f-1f8ec9d50811
# ╟─631833b8-0649-4f45-b1dd-e4a53091d321
# ╠═f8d3f48a-904f-11eb-2659-51f7508b434d
# ╠═9fdae9e2-6d2e-4901-8f48-14764b6800c2
# ╟─45bdd766-9131-11eb-31dc-f9ee1ced0db9
# ╟─f28f308a-4f5c-4dfe-8831-94de433c8fb9
# ╟─861b21ea-6fdd-4988-af50-b38f801de3a0
# ╟─a1af2f43-a6f4-4b5d-a882-abf91f109a44
# ╟─0eb3c064-9316-11eb-3c97-afba4053eb94
# ╟─4af1e07a-13e4-450f-85b6-cedc5cdeda52
# ╟─322216ee-93cb-11eb-272c-b7f5e5d87ee5
# ╟─8e63260a-0d0b-4551-96ee-88322e6e389f
# ╟─cadd0620-9315-11eb-228e-135bd066038d
# ╟─660c8c2d-6b3c-4a12-8e44-ea9dd1c461ff
# ╠═9cc0ece2-f2b9-4b54-9603-75e78d564b76
# ╠═d4e70ec0-b6b9-47a0-a575-6d15fc2d25e1
# ╠═871a129f-d016-4088-b1e6-c0a7d74e0393
# ╠═e7d591ac-fd39-45e6-8c61-e5c74cdfd597
# ╠═5a401dd6-05a4-47c1-af10-d67a4d9b40dc
# ╠═d2bc643b-c247-4f55-ad67-1992bf25a30f
# ╠═993eaa5a-91c5-4d7d-97d4-67993cb3d990
# ╠═82b9a9fc-feaa-4b67-98b8-577c8bdf3378
# ╠═8ad489e2-9318-11eb-1531-7515b778e6b8
# ╠═e47a9702-710b-44f4-9adf-136e04b3b04e
# ╠═2cb6d5aa-3bdd-4ce6-be0a-dcaac538b86f
# ╠═772e3bc5-a169-4c1d-aa22-2ea0720ef707
# ╠═a1bfbae2-b284-417d-953f-bc5a2493d56c
# ╠═8cf9ece5-c345-4e65-8426-2b71f9f7fea8
# ╠═ff5cdffb-0ce9-4f1b-a0ad-ce74e9d8f5e9
# ╠═627d372c-3b7a-45c0-a416-b8ffe4e5a02b
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
# ╠═a2c85efc-6aca-4794-844e-d9be1919c16e
# ╠═971f16da-62af-49b0-bb38-5f4c39b2b886
# ╠═4d5bdc05-0f5d-403f-b303-93e95ec408ca
# ╠═1ffd3dcd-87e7-4bee-8743-91d6bb35f077
# ╠═7dc9afbf-074a-4149-bbae-7b182916e7f2
# ╠═6a0ba3b4-ed64-4c21-bbcd-ad88685b6453
# ╠═9de1cf4b-a383-4bfe-a45b-351b94c75385
# ╟─9f06ce2d-f0b1-4ab6-931d-700d0176742a
# ╟─d99b60c2-e651-4ed5-9aec-44b3e410863b
# ╠═a6f9aef3-01ba-46aa-a960-41242e14a2c8
# ╠═c995629f-4728-460b-ab64-9fa47773ce58
# ╠═0c88b0b8-94b1-4ce6-88f2-0a813ed6fce4
# ╠═34e6fa38-9a1c-4c64-9287-60f7dbc7a17c
# ╠═e1c26455-1ba4-438b-9022-cb419dde4371
# ╠═52b773e3-a5d7-4f15-81d8-3607078b72eb
# ╠═097a9891-3c86-4162-a6b8-3c2ac94dc12b
# ╠═0e2a6d5f-ec49-4460-ab87-8ba293032d40
# ╟─bca3578d-3e45-4a31-894e-94aac3d20839
# ╟─6651e306-f2fb-49ae-8632-9481f5f4b246
# ╠═b0a7893f-856f-42ec-bf8b-22a4ca2907d4
# ╠═dbd925fd-37b3-46fe-ac4b-883681dcb915
# ╠═87fdfcaa-ce6d-4029-b94e-8e9441399a54
# ╠═04d1d9f5-87dc-40f0-a609-c3a7e2f98b00
# ╠═fb2f1b5e-2af2-403c-8fa1-0de77669887b
# ╠═88a52fc2-2cac-4bed-bb02-a77d6a678c6f
# ╠═7b6bf81f-a919-4d9f-a16d-e97f2dc3af95
# ╟─956c71a3-bf0c-49cc-849a-71d5bd701ebc
# ╠═4495993d-115f-4e20-9806-8793cf40d777
# ╠═bab49c82-ca0c-4f25-9d4a-4ca38b25054d
# ╟─ff9e4c96-1353-423f-af2b-8fc92b04b5da

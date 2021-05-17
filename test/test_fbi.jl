
function test1_df_functions(w, fbis)
	L = collect(wl)
	test1 = [(F.f).(L) / (F.N) ≈ (F.pdf).(L) for F in fbis]
	return all(test1)
end

function test2_df_functions(w, fbis)
	L = collect(wl)
	LBL =["FBIG1", "FBIBaG1", "FBIG2", "FBIBaG2"]
	test2 =[(F.f).(L) ≈ fbidf[!,LBL[i]] for (i,F) in enumerate(fbis)]
	return all(test2)
end


function test_fbi_pdfs(wr, fbis)
	Ll = collect(wr)
	test = [ qpdf(f.pdf, wr[1], wr[end]) ≈ 1 for f in fbis]
	return all(test)
end


function test_maximum_fbi(wl, fbis)
	L = collect(wl)
	LBL =["FBIG1", "FBIBaG1", "FBIG2", "FBIBaG2"]
	test =[maximum(F.(L)) ≈  find_max_xy(fbidf, "λ",
	                                    LBL[i])[1] for (i,F) in enumerate(fbis)]
	return all(test)
end


function fname(flt)
	return split(path_from_name(flt, datadir("filters")),"/")[end]
end

function check_file(path, flt)
	for f in cd(readdir, path)
		if f == flt
			return true
		end
	end
	return false
end

fbidf = load_df_from_csv(datadir("fluorimeter/325nm"),
	                     "FBI_G1_vs_FBI_G2_em325nm.csv",  spG)

wf = 350.
we = 800.
ws = 2.
wl = wf:ws:we

λminbp = 425.0  # in nm
λmaxbp = 435.0  # in nm
λminlp = 450.0  # in nm
λmaxlp = 800.0  # in nm


fbig1   =  dftof(wl, fbidf, "FBIG1")
fbibag1 =  dftof(wl, fbidf, "FBIBaG1")
fbig2   =  dftof(wl, fbidf, "FBIG2")
fbibag2 =  dftof(wl, fbidf, "FBIBaG2")

fbigfg1   =  dftogf(wl, fbidf, "FBIG1")
fbibagfg1 =  dftogf(wl, fbidf, "FBIBaG1")
fbigfg2   =  dftogf(wl, fbidf, "FBIG2")
fbibagfg2 =  dftogf(wl, fbidf, "FBIBaG2")

Phi = fom(fbigfg1, fbibagfg1, λminbp, λmaxbp)

adf = load_df_from_csv(datadir("fbi"),
					   "molar_extinction_coefficient_G1G2.csv", enG)
gs = ["g1", "g2"]
ls = [325,405]
qs = [0.67,0.67]

fbiFluo = fbi_fluorophores(adf, gs, ls, qs)

lp    = 0.1
cs    = 5E-05
l400  = Laser(405.0nm, lp*mW)
l325  = Laser(325.0nm, lp*mW)
obj   = Objective("Topatu", 0.6, 75.0)
gl400 = GaussianLaser(l400, obj)
gl325 = GaussianLaser(l325, obj)
dl400 = diffraction_limit(gl400)
dl325 = diffraction_limit(gl325)
f400  = Fov(dl400, 2*gl400.zr)
f325  = Fov(dl325, 2*gl325.zr)
I400  = photon_density(405nm, lp*mW, f400.a)
I325  = photon_density(325nm, lp*mW, f325.a)

lT   = LaserSetup(Dict("l405" => l400, "l325" => l325), lp*mW, obj)
solfbi = Solution("FBI solution", cs*M)

include("test_dffunctions.jl")
include("test_fluorophores.jl")
include("test_setup.jl")
include("test_filters.jl")
include("test_tools.jl")

# @testset "dffunctions_test" begin
# @test test_maximum_fbi(wl, [fbig1,fbibag1,fbig2,fbibag2])
# @test test1_df_functions(wl, [fbigfg1,fbibagfg1,fbigfg2,fbibagfg2])
# @test test2_df_functions(wl, [fbigfg1,fbibagfg1,fbigfg2,fbibagfg2])
# @test test_fbi_pdfs(wl, [fbigfg1,fbibagfg1,fbigfg2,fbibagfg2])
# end
#
# @testset "fbi_test" begin
# 	@test stn_band(fbibagfg1.f, fbigfg1.f, λminbp, λmaxbp) ≈ qpdf(fbibag1,
# 	λminbp, λmaxbp) / sqrt(qpdf(fbig1, λminbp, λmaxbp))
# 	@test stn_band(fbibagfg2.f, fbigfg2.f, λminbp, λmaxbp) ≈ qpdf(fbibag2,
# 	λminbp, λmaxbp) / sqrt(qpdf(fbig2, λminbp, λmaxbp))
# 	@test double_band_ratio(fbibagfg1.f, fbigfg1.f,
# 	                        λminbp, λmaxbp,
# 				            λminbp, λmaxbp) ≈ qpdf(fbibag1,
# 							λminbp, λmaxbp) / qpdf(fbig1, λminbp, λmaxbp)
# 	@test Phi.Φu ≈ qpdf(fbig1, λminbp, λmaxbp)
# 	@test Phi.ϵu ≈ qpdf(fbig1, λminbp, λmaxbp) / qpdf(fbig1, 0.0, 800.0)
# 	@test Phi.Φc ≈ qpdf(fbibag1, λminbp, λmaxbp)
# 	@test Phi.ϵc ≈ qpdf(fbibag1, λminbp, λmaxbp) / qpdf(fbibag1, 0.0, 800.0)
# 	@test Phi.Φuc ≈ double_band_ratio(fbibagfg1.f, fbigfg1.f,
# 	                                  λminbp, λmaxbp,
# 				                      λminbp, λmaxbp)
# 	@test Phi.fuc ≈ stn_band(fbibagfg1.f, fbigfg1.f, λminbp, λmaxbp)
# end

# @testset "fluorophores_test" begin
# 	fbi = Fluorophore("FBI", 325nm, 525nm, 13e+4*(M^-1*cm^-1), 0.6)
# 	fbit = Solution("FBI Standard", 1M)
# 	Pt = Powder("Powder test", 1mol/mg, 1mg/mm^2)
# 	ps = Powder("Standard", 1mmol/mg, 1mg/cm^2)
# 	pt = Pellet(1mm^2, 1mg)
# 	tml = Mlayer("test", N_A*mol/cm^2)
# 	pe = Pellet(1μm^2, 1mg)
# 	ml = Mlayer("test", 1e+6/μm^2)
#
# 	@test unit(fbi.σ) == cm^2
# 	@test uconvert(mol, fbit.c*1L) == 1mol
# 	@test nofv(1M, 1L) ≈ N_A*1.0*mol
# 	@test nofv(fbit, 1L) ≈ N_A*1.0*mol
# 	@test nofa(Pt, pt, 1mm^2) ≈ N_A*1.0*mol
# 	@test nofa(Pt,1mm^2) ≈ N_A*1.0*mol
# 	@test pe.area * ml.σm ≈ 1e+6
# 	@test ps.cs * pe.mass ≈ 1mmol
# 	@test nofa(tml, 1cm^2) ≈ N_A*1.0*mol
# end
#
# @testset "FbiFluorophores_test" begin
#
# 	@test keys(fbiFluo.fbi) == keys(fbiFluo.fbiba)
# 	#
# 	for gn in ["g1", "g2"]
# 		efbi    = string("ϵFbi",titlecase(gn))
# 		efbiba  = string("ϵFbiBa",titlecase(gn))
# 		for l in ls
# 			lfbi = string("l", l, gn)
# 			@test fbiFluo.fbi[lfbi].ϵ*(M*cm) ≈ select_element(adf,"λ",l,efbi)
# 			@test fbiFluo.fbiba[lfbi].ϵ*(M*cm) ≈ select_element(adf,"λ",l,efbiba)
# 	   	end
# 	end
#
#
# end

# @testset "setup_test" begin
# 	fov = Fov(1mm, 1mm)
# 	ohna = Objective("High NA", 0.9, 100.0)
# 	l405 = Laser(405nm, 1mW)
# 	gl =   GaussianLaser(l405,ohna)
#
# 	@test fov.a ≈ π * (fov.d/2.0)^2
# 	@test fov.v ≈ fov.a * fov.z
#
# 	@test l405.λ == 405.0nm
# 	@test l405.P == 1.0mW
#
# 	@test gl.w0/gl.zr ≈ gl.obj.NA
# 	gz0(r)  = gf(gl, 0*nm, r*nm)
# 	gr0(z)  = gf(gl, z*nm, 0*nm)
# 	gir0(z) = gI(gl, z*nm, 0*nm)
# 	giz0(r) = gI(gl, 0*nm, r*nm)
# 	rg0(x)  = gir0(x) / gr0(x)
#
# 	@test diffraction_limit(l405, ohna) ≈ diffraction_limit(gl)
#
# 	@test qpdf(giz0, 0.0, 100.0) ≈ qpdf(gz0, 0.0, 100.0)
# 	@test qpdf(rg0, 0.0, 100.0) ≈ qpdf(gr0, 0.0, 100.0)
#
# 	@test uconvert(mJ, delivered_energy(l405, 1s)) ≈ 1mJ
# 	@test n_photons(405nm, 1mW) ≈ n_photons(l405)
# 	@test photon_density(l405, fov) ≈ photon_density(405nm,
# 	                                                 1mW, fov.a)
#
# 	@test transmission(0.9) ≈ transmission(ohna)
# end
#
# @testset "lasersetup_test" begin
#
# 	@test lT.Is["l325"] ≈ I325
# 	@test lT.Is["l405"] ≈ I400
#
# end

# @testset "fluorescence_test" begin
# 	efluo = emitted_fluorescence(fbiFluo, lT)
# 	@test keys(fbiFluo.fbi) == keys(efluo.fbi)
# 	@test keys(fbiFluo.fbiba) == keys(efluo.fbiba)
# end

# @testset "filters_test" begin
# 	flt = split(path_from_name("BandPass405",
# 	      datadir("filters")),"/")[end]
# 	@test all([check_file(datadir("filters"),
#                               flt) for flt in fname.(fnames)])
# 	Filters = load_filter.(fnames, (datadir("filters"),))
# 	fdfs, fints = collect(zip(Filters...))
# 	fdfl = [length(f[!, "λ"]) for f in fdfs]
# 	fdfl2 = [length(collect(franges[i])) for (i,r) in enumerate(franges)]
# 	dx = abs.(fdfl .-fdfl2)
# 	@test length(dx[dx.>1]) == 0
#
# 	fbp405, flp425, fbp430, flp450, fdn405_522, fn405 = fints
# 	filterSetBand430(λ) = map(x->flp425(x)^2 * fdn405_522(x) * fbp430(x), λ)
# 	filterSetLP450(λ) = map(x->flp425(x)^2 * fdn405_522(x) * flp450(x), λ)
# 	f430 = flp425(430.0)^2 * fdn405_522(430.0) * fbp430(430.0)
# 	f470 = flp425(470.0)^2 * fdn405_522(470.0) * flp450(470.0)
#
# 	@test f430 ≈ filterSetBand430(430.0)
# 	@test f470 ≈ filterSetLP450(470.0)
# end
#

println("Tests run successfully")


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
	test =[maximum(F.(L)) ≈  find_max_xy(fbidf, "λ", LBL[i])[1] for (i,F) in enumerate(fbis)]
	return all(test)
end

fbidf = load_df_from_csv(datadir("fluorimeter/325nm"),
	                     "FBI_G1_vs_FBI_G2_em325nm.csv",  spG)

wf = 350.
we = 800.
ws = 2.
wl = wf:ws:we

fbig1   =  dftof(wl, fbidf, "FBIG1")
fbibag1 =  dftof(wl, fbidf, "FBIBaG1")
fbig2   =  dftof(wl, fbidf, "FBIG2")
fbibag2 =  dftof(wl, fbidf, "FBIBaG2")

fbigfg1   =  dftogf(wl, fbidf, "FBIG1")
fbibagfg1 =  dftogf(wl, fbidf, "FBIBaG1")
fbigfg2   =  dftogf(wl, fbidf, "FBIG2")
fbibagfg2 =  dftogf(wl, fbidf, "FBIBaG2")

@testset "dffunctions_test" begin
@test test_maximum_fbi(wl, [fbig1,fbibag1,fbig2,fbibag2])
@test test1_df_functions(wl, [fbigfg1,fbibagfg1,fbigfg2,fbibagfg2])
@test test2_df_functions(wl, [fbigfg1,fbibagfg1,fbigfg2,fbibagfg2])
@test test_fbi_pdfs(wl, [fbigfg1,fbibagfg1,fbigfg2,fbibagfg2])
end

println("Tests run successfully")

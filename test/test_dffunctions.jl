@testset "dffunctions_test" begin
@test test_maximum_fbi(wl, [fbig1,fbibag1,fbig2,fbibag2])
@test test1_df_functions(wl, [fbigfg1,fbibagfg1,fbigfg2,fbibagfg2])
@test test2_df_functions(wl, [fbigfg1,fbibagfg1,fbigfg2,fbibagfg2])
@test test_fbi_pdfs(wl, [fbigfg1,fbibagfg1,fbigfg2,fbibagfg2])
end

@testset "fbi_test" begin
	@test stn_band(fbibagfg1.f, fbigfg1.f, λminbp, λmaxbp) ≈ qpdf(fbibag1,
	λminbp, λmaxbp) / sqrt(qpdf(fbig1, λminbp, λmaxbp))
	@test stn_band(fbibagfg2.f, fbigfg2.f, λminbp, λmaxbp) ≈ qpdf(fbibag2,
	λminbp, λmaxbp) / sqrt(qpdf(fbig2, λminbp, λmaxbp))
	@test double_band_ratio(fbibagfg1.f, fbigfg1.f,
	                        λminbp, λmaxbp,
				            λminbp, λmaxbp) ≈ qpdf(fbibag1,
							λminbp, λmaxbp) / qpdf(fbig1, λminbp, λmaxbp)
	@test Phi.Φu ≈ qpdf(fbig1, λminbp, λmaxbp)
	@test Phi.ϵu ≈ qpdf(fbig1, λminbp, λmaxbp) / qpdf(fbig1, 0.0, 800.0)
	@test Phi.Φc ≈ qpdf(fbibag1, λminbp, λmaxbp)
	@test Phi.ϵc ≈ qpdf(fbibag1, λminbp, λmaxbp) / qpdf(fbibag1, 0.0, 800.0)
	@test Phi.Φuc ≈ double_band_ratio(fbibagfg1.f, fbigfg1.f,
	                                  λminbp, λmaxbp,
				                      λminbp, λmaxbp)
	@test Phi.fuc ≈ stn_band(fbibagfg1.f, fbigfg1.f, λminbp, λmaxbp)
end

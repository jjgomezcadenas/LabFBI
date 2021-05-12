@testset "filters_test" begin
	flt = split(path_from_name("BandPass405",
	      datadir("filters")),"/")[end]
	@test all([check_file(datadir("filters"),
                              flt) for flt in fname.(fnames)])
	Filters = load_filter.(fnames, (datadir("filters"),))
	fdfs, fints = collect(zip(Filters...))
	fdfl = [length(f[!, "λ"]) for f in fdfs]
	fdfl2 = [length(collect(franges[i])) for (i,r) in enumerate(franges)]
	dx = abs.(fdfl .-fdfl2)
	@test length(dx[dx.>1]) == 0

	fbp405, flp425, fbp430, flp450, fdn405_522, fn405 = fints
	filterSetBand430(λ) = map(x->flp425(x)^2 * fdn405_522(x) * fbp430(x), λ)
	filterSetLP450(λ) = map(x->flp425(x)^2 * fdn405_522(x) * flp450(x), λ)
	f430 = flp425(430.0)^2 * fdn405_522(430.0) * fbp430(430.0)
	f470 = flp425(470.0)^2 * fdn405_522(470.0) * flp450(470.0)

	@test f430 ≈ filterSetBand430(430.0)
	@test f470 ≈ filterSetLP450(470.0)
end

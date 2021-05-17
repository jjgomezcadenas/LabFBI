

@testset "uitls_test" begin
	@test to_fstr(π, "%7.2f")== @sprintf "%7.2f" π
	@test to_fstr(10^π, "%7.3g") =="1.39e+03"
	@test to_fstr(1, "%d") =="1"
	@test to_fstr("cat", "%s") =="cat"
	@test vect_to_fstr([1,2,3,4], "%d") == "1, 2, 3, 4"
	@test vect_to_fstr([π,log(π),π^2,sqrt(π)], "%7.4f") ==" 3.1416,  1.1447,  9.8696,  1.7725"
	@test vect_to_fstr(["the", "sun", "shines"], "%s") == "the, sun, shines"
	c1 = collect(logrange(1,10,10))
	c2 = collect((10^y for y in range(log10(1), log10(10), length=10)))
	@test c1 ≈ c2
end

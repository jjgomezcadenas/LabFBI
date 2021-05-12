@testset "setup_test" begin
	fov = Fov(1mm, 1mm)
	ohna = Objective("High NA", 0.9, 100.0)
	l405 = Laser(405nm, 1mW)
	gl =   GaussianLaser(l405,ohna)

	@test fov.a ≈ π * (fov.d/2.0)^2
	@test fov.v ≈ fov.a * fov.z

	@test l405.λ == 405.0nm
	@test l405.P == 1.0mW

	@test gl.w0/gl.zr ≈ gl.obj.NA
	gz0(r)  = gf(gl, 0*nm, r*nm)
	gr0(z)  = gf(gl, z*nm, 0*nm)
	gir0(z) = gI(gl, z*nm, 0*nm)
	giz0(r) = gI(gl, 0*nm, r*nm)
	rg0(x)  = gir0(x) / gr0(x)

	@test diffraction_limit(l405, ohna) ≈ diffraction_limit(gl)

	@test qpdf(giz0, 0.0, 100.0) ≈ qpdf(gz0, 0.0, 100.0)
	@test qpdf(rg0, 0.0, 100.0) ≈ qpdf(gr0, 0.0, 100.0)

	@test uconvert(mJ, delivered_energy(l405, 1s)) ≈ 1mJ
	@test n_photons(405nm, 1mW) ≈ n_photons(l405)
	@test photon_density(l405, fov) ≈ photon_density(405nm,
	                                                 1mW, fov.a)

	@test transmission(0.9) ≈ transmission(ohna)
end

@testset "lasersetup_test" begin

	@test lT.Is["l325"] ≈ I325
	@test lT.Is["l405"] ≈ I400

end

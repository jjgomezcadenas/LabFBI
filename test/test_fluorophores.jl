


@testset "fluorophores_test" begin
@test fl.nofv(1M, 1L) ≈ N_A*1.0*mol
end



println("Tests run successfully")

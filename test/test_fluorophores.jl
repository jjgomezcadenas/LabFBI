@testset "fluorophores_test" begin
	fbi = Fluorophore("FBI", 325nm, 525nm, 13e+4*(M^-1*cm^-1), 0.6)
	fbit = Solution("FBI Standard", 1M)
	Pt = Powder("Powder test", 1mol/mg, 1mg/mm^2)
	ps = Powder("Standard", 1mmol/mg, 1mg/cm^2)
	pt = Pellet(1mm^2, 1mg)
	tml = Mlayer("test", N_A*mol/cm^2)
	pe = Pellet(1μm^2, 1mg)
	ml = Mlayer("test", 1e+6/μm^2)

	@test unit(fbi.σ) == cm^2
	@test uconvert(mol, fbit.c*1L) == 1mol
	@test nofv(1M, 1L) ≈ N_A*1.0*mol
	@test nofv(fbit, 1L) ≈ N_A*1.0*mol
	@test nofa(Pt, pt, 1mm^2) ≈ N_A*1.0*mol
	@test nofa(Pt,1mm^2) ≈ N_A*1.0*mol
	@test pe.area * ml.σm ≈ 1e+6
	@test ps.cs * pe.mass ≈ 1mmol
	@test nofa(tml, 1cm^2) ≈ N_A*1.0*mol
end

@testset "FbiFluorophores_test" begin

	@test keys(fbiFluo.fbi) == keys(fbiFluo.fbiba)
	#
	for gn in ["g1", "g2"]
		efbi    = string("ϵFbi",titlecase(gn))
		efbiba  = string("ϵFbiBa",titlecase(gn))
		for l in ls
			lfbi = string("l", l, gn)
			@test fbiFluo.fbi[lfbi].ϵ*(M*cm) ≈ select_element(adf,"λ",l,efbi)
			@test fbiFluo.fbiba[lfbi].ϵ*(M*cm) ≈ select_element(adf,"λ",l,efbiba)
	   	end
	end

end


@testset "iruFluorophores_test" begin

	rdf = load_df_from_csv(datadir("fbi"),
	                       "molar_extinction_coefficient_RuIr.csv",
	                       enG)
	#
	iruFluo = iru_fluorophores(rdf,
	                       ["IrF", "Ir", "Ru", "IrF+","Ir++","Ru++"],
				           [325,405],
				           [0.01,0.01,0.01,0.01,0.01, 0.01])

	@test iruFluo["l325Ru"].name == "ϵRu"
	@test iruFluo["l325Ru"].ϵ    == 9339.0 * 1.0/(cm*M)
	@test iruFluo["l405Ru++"].name == "ϵRu++"
	@test iruFluo["l405Ru++"].ϵ    == 7293.0 * 1.0/(cm*M)
end

@testset "fluorescence_test" begin
	efluo = emitted_fluorescence(fbiFluo, lT)
	mfov =  molecules_in_fov(lT, solfbi)
	ffov =  emitted_fluorescence_fov(efluo, mfov)
	ofilter = transmission(obj)
	fobj =  filtered_fluorescence(ffov, ofilter)

	@test keys(fbiFluo.fbi) == keys(efluo.fbi)
	@test keys(fbiFluo.fbiba) == keys(efluo.fbiba)

	for gn in gs
		for l in ls
			klt = string("l",l)
			kff = string(klt,gn)
			@test ffov.fbi[kff] ≈  efluo.fbi[kff] * mfov[klt]
			@test fobj.fbi[kff] ≈  ffov.fbi[kff] * ofilter
			@test ffov.fbiba[kff] ≈  efluo.fbiba[kff] * mfov[klt]
			@test fobj.fbiba[kff] ≈  ffov.fbiba[kff] * ofilter

		end
	end
end

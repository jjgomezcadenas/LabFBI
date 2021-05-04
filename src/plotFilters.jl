using Plots

markercolors = [:green :orange :black :purple :red  :yellow :brown :white]

function plot_filter(filtername, fdf, fint, wl)
	function labels()
		xlabel!("λ (nm)")
		ylabel!("T")
		title!(filtername)
	end

	ifbi = [fint(x*1.0) for x in wl]
	p1 = plot(fdf.λ, fdf.T,
		label = "Filter data",
    	shape = :circle,
    	color = :black,
    	markersize = 3, leg=:topright)
	p2 = plot!(p1,wl, ifbi,
		label = "Filter data interpolation",
    	color = :blue,
    	markersize = 3, leg=:topright)
	labels()
	p3 = plot(wl, ifbi,
		label = "Filter data interpolation",
    	color = :blue,
    	markersize = 3, leg=:topright)
	labels()

	plot(p2, p3, layout = (1, 2), legend = false, fmt = :png)
end

function plot_filters(filterlist, filternames, wl)
	function labels(filtername)
		xlabel!("λ (nm)")
		ylabel!("T")
		title!(filtername)
	end
	PP = []
	for (i, ff) in enumerate(filterlist)
		ifbi = [ff(x*1.0) for x in wl]

		p = plot(wl, ifbi,
				 label = filternames[i],
    			 color = markercolors[i],
    			 lw=2, leg=:topright)
		labels(filternames[i])
		push!(PP,p)
	end

	return PP
end

function plot_filterset(fset, fname, wl)
	function labels()
		xlabel!("λ (nm)")
		ylabel!("T")
		title!(fname)
	end

	ifbi = [fset(x) for x in wl]
	p    = plot(wl, ifbi,
		label = fname,
    	color = :blue,
    	lw=2, leg=:false)
	labels()
end

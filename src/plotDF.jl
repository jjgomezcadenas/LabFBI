using Plots
using DataFrames

markercolors = [:green :orange :black :purple :red  :yellow :brown :white]

"""
	plot_xy(x, y, xlabel, ylabel, title;
				 color = :black,
				 lw = 2,
				 linestyle = :solid,
				 legend=false)

"""
function plot_xy(x, y, xlabel, ylabel, title;
	             color = :black,
				 lw = 2,
				 linestyle = :solid,
				 legend=false)

	p1 = plot(x, y,
		      lw=lw,
		      color=color,
		      linestyle = linestyle,
		      legend=legend,
			  fmt = png)
	xlabel!(xlabel)
	ylabel!(ylabel)
	title!(title)
	return p1
end

"""
	function log_xy(x, y, xlabel, ylabel, title;
	             color = :black,
				 lw = 2,
				 linestyle = :solid,
				 legend=false)

"""
function log_xy(x, y, xlabel, ylabel, title;
	             color = :black,
				 lw = 2,
				 linestyle = :solid,
				 legend=false)

	p1 = plot(x, y,
			  xaxis=:log,
		      lw=lw,
		      color=color,
		      linestyle = linestyle,
		      legend=legend,
			  fmt = png)
	xlabel!(xlabel)
	ylabel!(ylabel)
	title!(title)
	return p1
end

"""
	loglog_xy(x, y, xlabel, ylabel, title;
				 color = :black,
				 lw = 2,
				 linestyle = :solid,
				 legend=false)

"""
function loglog_xy(x, y, xlabel, ylabel, title;
	             color = :black,
				 lw = 2,
				 linestyle = :solid,
				 legend=false)

	p1 = plot(x, y,
			  xaxis=:log,
			  yaxis=:log,
		      lw=lw,
		      color=color,
		      linestyle = linestyle,
		      legend=legend,
			  fmt = png)
	xlabel!(xlabel)
	ylabel!(ylabel)
	title!(title)
	return p1
end

"""
	scatter_xy(x, y, xlabel, ylabel, title;
				 shape = :circle, color = :black,
				 markersize = 2, legend=false)

"""
function scatter_xy(x, y, xlabel, ylabel, title;
	             shape = :circle, color = :black,
				 markersize = 2, legend=false)

	p1 = scatter(x, y,
		      shape=shape,
		      color=color,
		      markersize=markersize,
		      legend=legend,
			  fmt = png)
	xlabel!(xlabel)
	ylabel!(ylabel)
	title!(title)
	return p1
end

"""
	plotdf_gfs(gfs, wl::AbstractRange,
					labels::Array{String}, colors,
					xlabel::String, ylabel::String, title::String;
					pdf::Bool=true,
					lw = 2,
					legend=false)

Plot (x,gy) where gy is a general function. If pdf = true, the pdf of the GF
is used, otherwise the fn (not normalized function)
"""
function plotdf_gfs(gfs, wl::AbstractRange,
	                labels::Array{String}, colors,
	                xlabel::String, ylabel::String, title::String;
					pdf::Bool=true,
					lw = 2,
		            legend=false)

	w = collect(wl)

	if pdf
		ff =[gf.pdf for gf in gfs]
	else
		ff =[gf.f for gf in gfs]
	end

	p1 = plot(w, ff[1].(w),
		      label=labels[1],
		      color=colors[1],
		      lw=lw,
		      legend=legend,
			  fmt = png)

	for (i, gf) in enumerate(ff[2:end])
		p1 = plot!(p1, w, gf.(w),
			       label=labels[i+1],
		           color=colors[i+1],
		           lw=lw,
		           legend=legend,
			       fmt = png)
	end

	xlabel!(xlabel)
	ylabel!(ylabel)
	title!(title)
	return p1
end


"""
	plotdf_xy(df::DataFrame, dfx::String, dfy::String,
		 	  xlabel, ylabel;
		 	  label="DF data", shape = :circle, color = :black,
		 	  markersize = 3, legend=false)

Plot columns (dx,dy) of dataframe df

"""
function plotdf_xy(df::DataFrame, dfx::String, dfy::String,
	              xlabel, ylabel;
	              label="DF", shape = :circle, color = :black,
				  markersize = 2, legend=false)

	p1 = plot(df[!,dfx], df[!,dfy],
		      label=label,
		      shape=shape,
		      color=color,
		      markersize=markersize,
		      legend=legend,
			  fmt = png)
	xlabel!(xlabel)
	ylabel!(ylabel)
	title!(dfy)
	return p1
end

"""
	plotdf_xys(df::DataFrame, dfx::String, dfys::Array{String},
					labels::Array{String}, colors,
					xlabel, ylabel, title;
					shape = :circle,
					markersize = 2,
					legend=false)

Plot  (dx,cols) of dataframe df, where cols is a list of columns

"""
function plotdf_xys(df::DataFrame, dfx::String, dfys::Vector{String},
					norm::Bool,
	                labels::Vector{String}, colors,
	                xlabel, ylabel, title;
	                shape = :circle,
				    markersize = 2,
		            legend=false)

	function csum(dfc)
		if norm == true
			xsum = sum(dfc)
		else
			xsum = 1.0
		end
		return xsum
	end


	p1 = plot(df[!,dfx], df[!,dfys[1]]/csum(df[!,dfys[1]]),
		      label=labels[1],
		      shape=shape,
		      color=colors[1],
		      markersize=markersize,
		      legend=legend,
			  fmt = png)

	for (i, dfy) in enumerate(dfys[2:end])
		p1 = plot!(p1, df[!,dfx], df[!,dfy]/csum(df[!,dfy]),
			      label=labels[i+1],
			      shape=shape,
			      color=colors[i+1],
			      markersize=markersize,
			      legend=legend,
				  fmt = png)
	end

	xlabel!(xlabel)
	ylabel!(ylabel)
	title!(title)
	return p1
end


"""
	plotf(f, X, xlabel, ylabel, title;
		  label="DF", lw = 2, color = :black, legend=false)

Plot f.(X) where X is a vector
"""
function plotf(f, X, xlabel, ylabel, title;
	           label="DF", lw = 2, color = :black, legend=false)
	fx = f.(X)
	p = plot(X, fx,
		     lw = lw,
    	     colour = color,
    	     label  = label,
    	     legend = legend,
	   	     fmt = :png)

	xlabel!(xlabel)
	ylabel!(ylabel)
	title!(title)
	return p
end

function plot_data_and_model(x, y, fx, xlabel, ylabel, title;
	                         lw = 2,
							 color = :black,
							 markersize=2,
	                         legend=false)
	p1 = scatter(x, y,
	             label="data",
				 markersize=markersize,
				 colour = color,
	    	     legend = legend,
		   	     fmt = :png)

	p1 = plot!(p1, x, fx,
			   lw = lw,
			   colour = color,
			   label  = "f(x)",
			   legend = legend,
			   fmt = :png)

	xlabel!(xlabel)
	ylabel!(ylabel)
	title!(title)
	return p1
end

function merge_plots!(sp1::Plots.Subplot, sp2::Plots.Subplot)
    append!(sp1.series_list, sp2.series_list)
    Plots.expand_extrema!(sp1[:xaxis], xlims(sp2))
    Plots.expand_extrema!(sp1[:yaxis], ylims(sp2))
    Plots.expand_extrema!(sp1[:zaxis], zlims(sp2))
    return sp1
end

"""
	merge_plots!(plt, plts...)

Merge plts plts with plt.
"""
function merge_plots!(plt, plts...)
	nplt = deepcopy(plt)
    for (i, sp) in enumerate(nplt.subplots)
        for other_plt in plts
            if i in eachindex(other_plt.subplots)
                merge_plots!(sp, other_plt[i])
            end
        end
    end
    return nplt
end

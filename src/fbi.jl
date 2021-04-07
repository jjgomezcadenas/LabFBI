using CSV
using DataFrames
using Interpolations
using QuadGK
using DrWatson


function load_fbidf(fbidir::String, fbiname::String)
	path       = string(datadir(),"/", fbidir)
    fbifname   = string(path, "/", fbiname)
    fbif       = CSV.File(fbifname; delim=';', decimal=',')
    return DataFrame(fbif)
end


function gfpdf(fi, xmin::Float64, xmax::Float64)
	function fn(x)
	   	if x < xmin || x > xmax
			return 0
		else
			return fi(x)
		end
	end
	return fn
end

function dftof(wl::AbstractRange, df::DataFrame, column::String)
    fbili = LinearInterpolation(wl, df[!,column])
	return gfpdf(fbili, wl[1], wl[end])
end


function ftopdf(wl::AbstractRange, f::Function)
	function pdf(x)
		return f(x) / N
	end
    N = quadgk(f, wl[1], wl[end])[1]
	return N, pdf
end


function qpdf(pdf::Function, λmin::Float64, λmax::Float64)
	return quadgk(pdf, λmin, λmax)[1]
end


function eps_band(fbipdf::Function, fbibapdf::Function,
	              λmin::Float64, λmax::Float64)
	eps_c = qpdf(fbibapdf, λmin, λmax)
	eps_u = qpdf(fbipdf, λmin, λmax)
	return eps_c / sqrt(eps_u)
end

using CSV
using DataFrames
using Interpolations
using QuadGK
using DrWatson

# fbiname   = "EDI044_FBI_round4_em325_375_405.csv"
# fbibaname = "EDI044_FBI_Ba_round4_em325_375_405.csv"
# fbidir    = "fluorimeter"

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

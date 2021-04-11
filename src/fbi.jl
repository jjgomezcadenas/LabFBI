using CSV
using DataFrames
using Interpolations
using QuadGK
using DrWatson

struct gf          # general function
	N  ::Float64
	f  ::Function
	pdf::Function
end

struct Φ
	Φu::Float64
	Φc::Float64
	ϵu::Float64
	ϵc::Float64
	Φuc::Float64
	fuc::Float64
end

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


function f_and_pdf_from_df(wr, df, column)
	fint   = dftof(wr, df, column)
	f      = gfpdf(fint, wr[1], wr[end])
	N, pdf = ftopdf(wr, f)
	return gf(N, f, pdf)
end


function fbigen(wr, df, columns)
	return [f_and_pdf_from_df(wr, df, c) for c in columns]
end


function fom(fbi, fbiba, λmin, λmax)
	Φu = qpdf(fbi.f, λmin, λmax)
	ϵu = Φu / fbi.N
	Φc = qpdf(fbiba.f, λmin, λmax)
	ϵc = Φc / fbiba.N
	Φuc = Φc / Φu
	fuc = Φc / sqrt(Φu)
	return Φ(Φu,Φc,ϵu,ϵc,Φuc,fuc)
end


function test_fbis(wr, fbis)
	Ll = collect(wr)
	test = [fbi.f.(Ll)/fbi.N ≈ fbi.pdf.(Ll) for fbi in fbis]
	return all(test)
end


function test_fbi_pdfs(wr, fbis)
	Ll = collect(wr)
	test = [quadgk(fbi.pdf, wr[1], wr[end])[1] ≈ 1 for fbi in fbis]
	return all(test)
end

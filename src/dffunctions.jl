# """
# Module dffunctions.jl provides tools to read data frames expressing
# data series relevant for LabFbi software
# (e.g, molecular cross sections, filters), and describe them when
# relevant as functions.
# """
using DataFrames
using CSV
using Interpolations
using QuadGK

"""
	CsvG

A struct which describes a simple grammar needed to read CSV files
written with spanish/english characters

# Fields
- `delim::Char`  : Delimiter used (e.g., ',' in english, ';' in spanish)
- `decimal::Char`: Symbol to represent decinal point (',' '.')
"""
struct CsvG  # stands for CSV Grammar
	delim::Char
	decimal::Char
end

spG = CsvG(';',',')  # spanish: decimals represented with ',' delimited with ';'
enG = CsvG(',','.')  # english


"""
	Gf

A struct representing a generalized function:

# Fields
- `N::Number`: Normalization constant: pdf(x) = f(x)/N.
- `f::Any`: Any integrable function.
- `pdf::Any`: Normalized to area 1 (PDF).
"""
struct Gf          # general function
	N  ::Number
	f  ::Any
	pdf::Any
end


"""
	load_df_from_csv(path::String, fname::String, csvg::CsvG)

Load a dataframe from a csv file.

# Arguments
- `path::String`: a path to the data.
- `fname::String`: name of the (csv) file expressing the df.
- `csvg=spG`: CVS coding (spanish, default) or english.
"""
function load_df_from_csv(path::String, fname::String, csvg::CsvG)
    name   = string(path, "/", fname)
    csvf   = CSV.File(name; delim=csvg.delim, decimal=csvg.decimal)
    return DataFrame(csvf)
end


"""
	dftof(wl, df::DataFrame, cname::String,
			   bkgnd::Float64=0.0)

Return an interpolated function, valid in range wl.

# Arguments
- `wl`: interpolation range.
- `df::DataFrame`: data frame holding the data.
- `cname::String`: name of the column holding data to be interplated.
- `bkgnd::Float64`: value of the data outside interpolation range.
"""
function dftof(wl, df::DataFrame, cname::String,
	           bkgnd::Float64=0.0)
    fbili = LinearInterpolation(wl, df[!,cname])
	return gfpdf_(fbili, wl[1], wl[end], bkgnd)
end


function gfpdf_(fi, xmin::Float64, xmax::Float64, bkgnd::Float64=0.0)
	function fn(x)
	   	if x < xmin || x > xmax
			return bkgnd
		else
			return fi(x)
		end
	end
	return fn
end


"""
	qpdf(f, λmin::Number, λmax::Number)

Return the integral of f in the interval (λmin, λmax)
(Syntactic sugar for quadgk)

# Arguments
- `f::Function`: Function to be integrated.
- `λmin::Number`: lower bound of range.
- `λmax::Number`: upper bound of range.
"""
function qpdf(f, λmin::Number, λmax::Number)
	return quadgk(f, λmin, λmax)[1]
end


"""
	ftopdf(wl, f)

Compute the PDF of function f in range wl.

# Arguments
- `wl`: Range of application.
- `f::Function`: Input function.
"""
function ftopdf(wl, f)
	function pdf(x)
		return f(x) / N
	end
    N = qpdf(f, wl[1], wl[end])
	return N, pdf
end


"""
	dftogf(wl, df::DataFrame, cname::String,
	       bkgnd::Float64=0.0)

Return a generalized function from dataframe, column and range

# Arguments
- `wl`: Range of application.
- `df::DataFrame`: data frame holding the data.
- `f::Function`: Input function.
- `cname::String`: name of the column holding data to be interplated.
- `bkgnd::Float64`: value of the data outside interpolation range
                    (zero by default).
"""
function dftogf(wl, df::DataFrame, cname::String,
	            bkgnd::Float64=0.0)
	f      = dftof(wl, df, cname, bkgnd)
	N, pdf = ftopdf(wl, f)
	return Gf(N, f, pdf)
end


"""
	find_max_xy(df::DataFrame, xc::String, yc::String)

Return ymax and x such that ymax = f(x).

# Description:
In a DataFrame one has often "XY" variables, that is, a pair of columns
"X" and "Y" which represent correlated variables (e.g intensity and wavelength).
In such cases one often wants the XY maximum, that is, finding the maximum
of "Y" (call it ymax) and the corresponding value in "X"
(x for which y is ymax). This corresponds, for example, to the wavelength at
which the intensity is maximal.

# Arguments
- `df::DataFrame`: data frame holding the data.
- `xc::String`: the X column.
- `yc::String`: the Y column.

"""
function find_max_xy(df::DataFrame, xc::String, yc::String)
	# function findmax returns value and index of max
	ymax, imax = findmax(df[!, yc])
	x_ymax = df[imax, xc]
	return ymax, x_ymax
end

"""
	select_row_from_column(df::DataFrame, column::String, value::Any)

Return an array of values corresponding to a row in which
column `column` has value `value`

# Arguments
- `df::DataFrame`: data frame holding the data.
- `xc::String`: the X column used for the selection.
- `value::Any`: the value used for the selection.
"""
function select_row_from_column(df::DataFrame, column::String, value::Any)
	return Array(df[df[!,column] .== value, :])
end


function select_element(df::DataFrame,
	                    column::String, value::Any,
						column2::String)
	return df[df[!,column] .== value, :][!, column2][1]
end

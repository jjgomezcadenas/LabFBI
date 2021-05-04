# """
# Module filters.jl includes functions to read the data associated
# to the LabFbi filters and return interpolated functions describing them.
# """
include("dffunctions.jl")

fnames = ["BandPass405", "LongPass425", "BandPass430",
          "LongPass450", "DoubleNotch405_522", "NF405"]
fdoc   = [
"Selects wavelengths in a narrow range around 405 nm",
"Selects wavelengths longer than 425 nm",
"Selects wavelengths in a narrow range around 430 nm",
"Selects wavelengths longer than 450 nm",
"Suppress wavelengths in two bands (405, 422) nm",
"Suppress wavelengths around 405 nm"
]

franges = [354.0:0.2:424.0,
           250.0:1.0:1000.0,
		   200.0:1.0:1000.0,
		   200.0:1.0:1000.0,
		   350.0:1.0:1000.0,
		   200.0:1.0:1000.0,
		   ]

fbkgnd = [1e-6, 0., 0., 0., 0., 0.]

filter_ranges = Dict( fnames[i] => franges[i] for (i, n) in enumerate(fnames))
filter_bkgnd  = Dict( fnames[i] => fbkgnd[i]  for (i, n) in enumerate(fnames))

#fpaths = path_from_names.(fnames, datadir("filters"))
#filter_names = Dict( fnames[i] => fpaths[i] for (i, n) in enumerate(fnames))


"""
	path_from_name(name::String, path::String, ext::String=".csv")
Return a full file path (and extension) fron name

# Arguments
- `name::String`     : name of the filter.
- `path::String`     : a path to the directory with filter data.
- `ext::String=csv`  : file extension

"""
function path_from_name(name::String, path::String, ext::String=".csv")
	return string(path,"/", name, ext)
end

"""
	function load_filter(name::String, path::String, ext::String=".csv")

Load a filter dataframe from a csv file. Returns a data frame and an
interpolated function 

# Arguments
- `name::String`       : name of the filter.
- `path::String`       : path to the filter folder.
- `ext::String=".csv"` : file extension.
"""
function load_filter(name::String, path::String, ext::String=".csv")
	function scaledf(df)
		df[!, :T] = map(t -> t/100.0, df[:, :T])
		return df
	end

	function fixdf(df, ci=3, cj=4)
		df = df[!,ci:cj]
		rename!(df, [Symbol("λ"),  Symbol("T")])
		return scaledf(df)
	end

	file    = path_from_name(name, path, ext)
	if name == "NF405"
		csvf    = CSV.File(file; delim=spG.delim, decimal=spG.decimal)
	else
		csvf    = CSV.File(file; delim=enG.delim, decimal=enG.decimal)
	end

	df = DataFrame(csvf)

	if name == "NF405"
		df = fixdf(df, 1,2)
		sort!(df, rev = false)
	elseif name == "BandPass430"
		df = fixdf(df)
		sort!(df, rev = false)
	elseif name == "DoubleNotch405_522"
		df = scaledf(df)
	elseif name == "LongPass425" || name == "LongPass450"
		df = fixdf(df)
	end

	df[!, :λ] = map(t -> t*1.0, df[:, :λ])
	fint = interpolate((df.λ,), df.T, Gridded(Linear()))
	return df, fint
end

#include("dffunctions.jl")
#include("math.jl")

"""
	mec_at_lamda(df::DataFrame, λ::Float64, scols::AbstractRange)

	Compute the Molar Extinction Coefficient (mec, ϵ) for wavelength λ,
	fitting a set of values (A, C) to the formula:

	A = C ϵ

	Where A is the absorbance (for a specified path length in cm),
	C is the concentration (in M) and ϵ is the mec in 1/(M cm)

	# Fields
	- `df::DataFrame`: A data frame with the following structure

		λ	C1	C2	C3...

		where C1, C2, C2 label the value of A (in 1/cm).
		Each column computes A for a different λ (in nm)
	- `λ::Float64`: ϵ is computed at this wavelength
	- `scols::AbstractRange`: select the number of columns in the df

"""
function mec_at_lamda(df::DataFrame, λ::Float64, scols::AbstractRange)
	eps     = select_row_from_column(df, "λ", λ)[scols]
	cnames  = names(df)[scols]
	cs      = parse.(Float64, replace.(cnames, "," => "."))
	return  cs, eps, lfit(cs, eps)
end


@doc raw"""
	struct Φ

A struct that characterizes the separation between free (u) and chelated (c) molecules.

# Fields
- `Φu::Float64` :
- `Φc::Float64` :
- `Φu::Float64` :
- `Φc::Float64` :
- `ϵu::Float64` :
- `ϵc::Float64` :
- `Φuc::Float64` :
- `fuc::Float64` :

The definition of the fields above is:

```math
\Phi_u = \int_{\lambda_{min}}{\lambda_{max}} f_u(\lambda) d\lambda \\
\Phi_c = \int_{\lambda_{min}}{\lambda_{max}} f_c(\lambda) d\lambda \\
\epsilon_u = \int_{\lambda_{min}}{\lambda_{max}} pdf_u(\lambda) d\lambda \\
\epsilon_c = \int_{\lambda_{min}}{\lambda_{max}} pdf_c(\lambda) d\lambda \\
\Phi_{uc} = \Phi_c / \Phi_u \\
f_{uc} = \Phi_c / \sqrt{Phi_u}
```
"""
struct Φ
	Φu::Float64
	Φc::Float64
	ϵu::Float64
	ϵc::Float64
	Φuc::Float64
	fuc::Float64
end


@doc raw"""
	stn_band(fs::Function, fn::Function, λmin::Number, λmax::Number)

Computes the signal to noise ratio between pdfs (signal) and pdfn (noise).
where:


# Arguments

- `fs::Function`   : input signal function
- `fn::Function`   : input noise function
- `λmin::Float64`  : minimum of range
- `λmax::Float64`  : maximum of range

The signal to noise definition is:


```math
	stn = S /\sqrt{N} \\
	S = \int_{\lambda_{min}}^{\lambda_{max}} fs(\lambda) d\lambda \\
	N = \int_{\lambda_{min}}^{\lambda_{max}} fn(\lambda) d\lambda \\
```
"""
function stn_band(fs::Function, fn::Function, λmin::Number, λmax::Number)
	S = qpdf(fs, λmin, λmax)
	N = qpdf(fn, λmin, λmax)
	return S / sqrt(N)
end


@doc raw"""
	double_band_ratio(fbp::Function, flp::Function,
			     	  λbpmin::Number, λbpmax::Number,
				      λlpmin::Number, λlpmax::Number)


# Arguments

- `fbp::Function`   : input fbp function (e.g, a function that peaks in the interval BP)
- `fl`::Function`   : input flp function (e.g, a function that peaks in the interval LP)
- `λbpmin::Float64`  : minimum of range BP
- `λbpmax::Float64`  : maximum of range BP
- `λlpmin::Float64`  : minimum of range LP
- `λlpmax::Float64`  : maximum of range LP

The function computes the ratio:
```math
	r = BP/LP \\
	BP = \int_{\lambda_{bpmin}}^{\lambda_{bpmax}} fbp(\lambda) d\lambda \\
	LP = \int_{\lambda_{lpmin}}^{\lambda_{lpmax}} flp(\lambda) d\lambda
```
"""
function double_band_ratio(fbp::Function, flp::Function,
	                       λbpmin::Number, λbpmax::Number,
				           λlpmin::Number, λlpmax::Number)
	bp = qpdf(fbp, λbpmin, λbpmax)
	lp = qpdf(flp, λlpmin, λlpmax)
	return bp/lp
end


"""
	fom(fbi::Gf, fbiba::Gf, λmin::Number, λmax::Number)

Return figure of merit Φ (see definition of type)

# Arguments

- `fbi::Gf`   : Generalized function defining the fbi spectrum
- `fbiba::Gf` : Generalized function defining the fbiba spectrum
- `λmin::Float64`  : minimum of range
- `λmax::Float64`  : maximum of range
"""
function fom(fbi::Gf, fbiba::Gf, λmin::Number, λmax::Number)
	Φu = qpdf(fbi.f, λmin, λmax)
	ϵu = Φu / fbi.N
	Φc = qpdf(fbiba.f, λmin, λmax)
	ϵc = Φc / fbiba.N
	Φuc = Φc / Φu
	fuc = Φc / sqrt(Φu)
	return Φ(Φu,Φc,ϵu,ϵc,Φuc,fuc)
end

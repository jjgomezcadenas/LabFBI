push!(LOAD_PATH,"../src/")
using DrWatson
using Test
@quickactivate "LabFBI"
include("../src/fbi.jl")
#include("../src/fluorophores.jl")
#include("../src/setup.jl")
include("../src/fluorescence.jl")
include("../src/utils.jl")
include("../src/filters.jl")

include("test_fbi.jl")

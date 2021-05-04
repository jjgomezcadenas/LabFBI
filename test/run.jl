push!(LOAD_PATH,"../src/")
using DrWatson
using Test
@quickactivate "LabFBI"
include("../src/fbi.jl")
include("../src/fluorophores.jl")
include("../src/setup.jl")
include("../src/utils.jl")
include("../src/filters.jl")
#include("test_dffunctions.jl")
include("test_fbi.jl")
include("test_tools.jl") #includes tests of dffunctions.jl

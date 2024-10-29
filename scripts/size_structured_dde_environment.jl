
# load environment: packages, functions 

pkgs_shared = [
    "Pkg", "DrWatson", "Revise", "Test",  "OhMyREPL", "Logging", 
    "StatsBase", "Statistics", "Distributions", "Random", "Setfield", "Memoization", 
    "MCMCChains", 
    "DataFrames", "JLD2", "CSV", "PlotThemes", "Colors", "ColorSchemes", "RData",  
    "Plots",  "StatsPlots", "MultivariateStats", 
    "ForwardDiff", "ReverseDiff", "Enzyme", "ADTypes", "Zygote", 
    "StaticArrays", "LazyArrays", "FillArrays", "LinearAlgebra", "MKL", "Turing"
]

pkgs = [
    "QuadGK", "ModelingToolkit", "DifferentialEquations", "Interpolations",
]

pkgs = unique!( [pkgs_shared; pkgs] )


println( """

Libraries loading:

If this is the initial run, it will take a while to precompile/download all libraries. 
The variable, 'pkgs', contains the list of required libraries. To (re)-install required packages, run: 

install_required_packages(pkgs)

""" )

# install_required_packages(pkgs)  # in case you need to install packages

for pk in pkgs;  @eval using $(Symbol(pk));  end  # needs to be run at the top level
 

include( srcdir( "shared_functions.jl" ) )   
include( srcdir( "size_structured_dde_functions.jl" ) ) 

println( """

Make sure the paths are correct:

Currently active project is: $(projectname())
Path of active project: $(projectdir())

""" )

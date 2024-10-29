# scripts used to create/initiate the project 
 
if false
    using Pkg
    # to do it manually with Pkg
    Pkg.generate("projects/model_fishery")
    Pkg.instantiate()
    Pkg.activate("projects/model_fishery")  
    cd("projects/model_fishery")
end

 
# or via DrWatson
using DrWatson
initialize_project("projects/model_fishery"; 
    authors="Jae S. Choi",
    template = [
        "data",  
        "docs",
        "media",
        "src",
        "scripts",
        "outputs",
        "papers",
        "test"
    ]
)

quickactivate("projects/model_fishery" )
cd( projectdir() )

# add some key libraries:

using Pkg 


if false
    # to run from a script:
    using DrWatson  # shouls ideally bethe first things to run
    quickactivate("projects/model_fishery" )
    include( srcdir( "model_fishery.jl") )
    # cd( projectdir() )

    # then do important stuff:

end

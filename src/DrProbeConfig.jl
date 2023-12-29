module DrProbeConfig

using Dates
using Comonicon
using DelimitedFiles
using LinearAlgebra
using DataStructures
using Images

export drprobe_config

include("drprobefunctions.jl")

#This docstring is required by Comonicon.jl
"""
Run a multislice simulation using the configuration file supplied. Under the hood, this runs 
the CELSLC and MSA programs by Dr. Juri Barthel. The execution is optionally multithreaded 
by running several (one for each scan line) MSA calculations in parallel. The parameters are
supplied in a configuration file in YAML format. This program does not check the validity 
of the configuration file to make sure it is correct and free of conflicting parameters! 

# Args

- `configuration-file`: configuration file for the calculation
"""
@main function drprobe_config(
    configuration_file::String
    )
    run_drprobe(configuration_file)
end

end # module DrProbeConfig

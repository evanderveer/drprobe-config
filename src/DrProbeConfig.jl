module DrProbeConfig

using Dates
using Comonicon

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

# Flags

- `-d, --debug`: turn on debug mode, does not run any calculations
- `-n, --no-cleanup`: skip cleaning up temporary files
- `-m, --multithreaded`: turn on multithreaded mode

# Options

- `-o, --output-folder=<path>`: folder in which to place the output, pwd by default

"""
@main function drprobe_config(
    configuration_file::String;
    output_folder::String = pwd(),
    debug::Bool = false,
    no_cleanup::Bool = false, 
    multithreaded::Bool = false
    )
    if multithreaded && Threads.nthreads() == 1
        throw(ArgumentError("multithreaded set to true, but nthreads == 1"))
    else
        run_drprobe(
            configuration_file, 
            output_folder=output_folder, 
            no_cleanup=no_cleanup, 
            debug=debug,
            multithreaded=multithreaded
            )
    end

end

end # module DrProbeConfig

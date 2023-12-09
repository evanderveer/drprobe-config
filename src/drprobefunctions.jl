using YAML 

include("celslc.jl")
include("msa.jl")
include("image.jl")
include("cellcalculator.jl")

const FILES_TO_KEEP = r"\.dat|\.sli"

"""
    run_drprobe(
        config_file::String; 
        output_folder::String = pwd(), 
        no_cleanup::Bool = false, 
        debug::Bool = false, 
        multithreaded::Bool = false
        )

Run the drprobe tool using the specified configuration file.

# Arguments
- `config_file::String`: Path to the configuration file in YAML format.
- `output_folder::String`: Path to the output folder (default is the current working directory).
- `no_cleanup::Bool`: Flag indicating whether to skip the cleanup process (default is false).
- `debug::Bool`: Flag indicating whether to enable debugging, no calculations are run (default is false).
- `multithreaded::Bool`: Flag indicating whether to use multithreading (default is false).

# Example
```julia
run_drprobe("config.yml", output_folder="/path/to/output", debug=true)
"""
function run_drprobe(
    config_file::String; 
    output_folder::String = pwd(), 
    no_cleanup::Bool = false,
    debug::Bool = false,
    multithreaded::Bool = false
    )

    start_dir = pwd()
    cd(dirname(config_file))
    config = YAML.load_file(basename(config_file))
    temp_dir = mktempdir(pwd())

    #Copy the input file (.cif/.cel) to the temporary folder
    cp(config["input"], joinpath(temp_dir, basename(config["input"])))
    cd(temp_dir)

    run_celslc(config, debug)

    msa_function = multithreaded ? run_msa_multithread : run_msa
    
    if config["run-msa"]
        make_detector_prm_file(config)
        msa_function(config, debug)
    else
        println("MSA deactivated, only running CELSLC")
    end

    cd(start_dir)
    
    if !(debug || no_cleanup)
        cleanup(config, joinpath(dirname(config_file), temp_dir), output_folder)
    end
end

"""
    cleanup(
        config::Dict, 
        temp_folder::String, 
        output_folder::String
        )

Clean up temporary files generated during the drprobe execution.

# Arguments
- `config::Dict`: Configuration dictionary.
- `temp_folder::String`: Path to the temporary folder.
- `output_folder::String`: Path to the output folder.

# Example
```julia
cleanup(config, "/path/to/temp", "/path/to/output")
"""
function cleanup(
    config::Dict, 
    temp_folder::String, 
    output_folder::String
    )
    println("Cleaning up the mess")
    output_folder = mkdir(string(now()))

    files_to_move = filter!(
                                x -> occursin(FILES_TO_KEEP, x), 
                                readdir(temp_folder)
                            )

    for file in files_to_move
        mv(
            joinpath(temp_folder, file), 
            joinpath(output_folder, file)
           )
    end
    rm(temp_folder, recursive=true)
end

using YAML 

include("celslc.jl")
include("msa.jl")
include("image.jl")
include("cellcalculator.jl")
include("imageconstructor.jl")


const FILES_TO_KEEP = r"\.dat|\.sli|\.tif"

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
    temporary_folder::String = pwd(), 
    slice_file_folder::String = pwd(),
    no_cleanup::Bool = false,
    debug::Bool = false,
    multithreaded::Bool = false
    )

    start_dir = pwd()
    cd(dirname(config_file))
    config = YAML.load_file(basename(config_file))

    #Make paths absolute to prevent problems later
    temp_dir = abspath(mktempdir(temporary_folder))
    slice_file_dir = abspath(slice_file_folder)
    output_dir = abspath(output_folder)

    #Copy the input file (.cif/.cel) to the temporary folder
    cp(config["input"], joinpath(temp_dir, basename(config["input"])))

    if config["run-celslc"]
        cd(temp_dir)
        run_celslc(config, debug)
    else
        println("CELSLC disabled, looking for slice files")
        link_slice_files(config, temp_dir, slice_file_dir)
        cd(temp_dir)
    end

    msa_function = multithreaded ? run_msa_multithread : run_msa
    
    if config["run-msa"]
        make_detector_prm_file(config)
        msa_function(config, debug)
        msa_spatial_convolution(config)
    else
        println("MSA deactivated, only running CELSLC")
    end

    make_images(config)

    cd(start_dir)
    
    if !(debug || no_cleanup)
        cleanup(config, temp_dir, output_dir)
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
    output_folder = mkdir(joinpath(output_folder, string(now())))
    println("Output folder: $output_folder")

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

function link_slice_files(
    config::Dict,
    temp_dir::String,
    slice_file_dir::String
)
    sli_file_regex = Regex(join([".*", config["output"], raw"_[0-9]+.sli$"]))
    full_file_list = readdir(slice_file_dir, join=true)
    sli_files = full_file_list[occursin.(sli_file_regex, full_file_list)]

    if length(sli_files) == 0
        error("No slice files found in $slice_file_dir. Make sure slice file name matches config file")
    end
    
    println("Found $(length(sli_files)) slice files, making symbolic links")

    #Slice files can get big so symlink instead of copy
    for sli_file in sli_files
        symlink(sli_file, joinpath(pwd(), temp_dir, basename(sli_file)))
    end
end


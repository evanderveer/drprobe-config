using YAML 

include("celslc.jl")
include("msa.jl")
include("image.jl")
include("cellcalculator.jl")
include("imageconstructor.jl")


const FILES_TO_KEEP = r"\.dat|\.sli|\.tif"

"""
    run_drprobe(config_file::String)

    Run the drprobe tool using the specified configuration file.

    # Arguments
    - `config_file::String`: Path to the configuration file in YAML format.

    # Example
    ```julia
    run_drprobe("./config.yml")
"""
function run_drprobe(config_file::String)

    #Make paths absolute to prevent problems later
    start_dir = abspath(pwd())
    config_file_dir = dirname(config_file)
    cd(config_file_dir == "" ? "." : config_file_dir)
    config = check_config_file(YAML.load_file(basename(config_file)))

    #Copy the input file (.cif/.cel) to the temporary folder
    cp(config["input"], joinpath(config["temporary-folder"], basename(config["input"])))

    if config["run-celslc"]
        cd(config["temporary-folder"])
        celslc(config)
    else
        println("CELSLC disabled, looking for slice files")
        link_slice_files(config)
        cd(config["temporary-folder"])
    end

    if config["run-msa"]
        make_detector_prm_file(config)
        msa(config)
        msa_spatial_convolution(config)
    else
        println("MSA deactivated, only running CELSLC")
    end

    make_images(config)

    cd(start_dir)
    cleanup(config)

end

function check_config_file(config::Dict)

    config["temporary-folder"] = config["temporary-folder"] == "" ? "." : config["temporary-folder"]
    config["temporary-folder"] = abspath(mktempdir(config["temporary-folder"]))

    config["slice-file-folder"] = config["slice-file-folder"] == "" ? "." : config["slice-file-folder"]
    config["slice-file-folder"] = abspath(config["slice-file-folder"])
end

"""
    cleanup(
        config::Dict, 
        temp_folder::String, 
        config["output-folder"]::String
        )

    Clean up temporary files generated during the drprobe execution.

    # Arguments
    - `config::Dict`: Configuration dictionary.

    # Example
    ```julia
    cleanup(config, "/path/to/temp", "/path/to/output")
"""
function cleanup(config::Dict)

    println("Cleaning up the mess")
    config["output-folder"] = mkdir(joinpath(config["output-folder"], string(now())))
    println("Output folder: $(config["output-folder"])")

    files_to_move = filter!(
                            x -> occursin(FILES_TO_KEEP, x), 
                            readdir(config["temporary-folder"])
                            )

    #Copy all files in debugging mode
    files_to_move = "no-cleanup" ∈ config["debug"] ?  readdir() : files_to_move

    for file in files_to_move
        mv(
           joinpath(config["temporary-folder"], file), 
           joinpath(config["output-folder"], file)
           )
    end 
end

function link_slice_files(config::Dict)
    
    sli_file_regex = Regex(join([".*", config["output"], raw"_[0-9]+.sli$"]))
    full_file_list = readdir(config["slice-file-folder"], join=true)
    sli_files = full_file_list[occursin.(sli_file_regex, full_file_list)]

    if length(sli_files) == 0
        error("No slice files found in $(config["slice-file-folder"]). Make sure slice file name matches config file")
    end
    
    println("Found $(length(sli_files)) slice files, making symbolic links")

    #Slice files can get big so symlink instead of copy
    for sli_file in sli_files
        symlink(sli_file, joinpath(config["temporary-folder"], basename(sli_file)))
    end
end


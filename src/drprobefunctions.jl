using YAML 

include("celslc.jl")
include("msa.jl")
include("image.jl")
include("cellcalculator.jl")
include("imageconstructor.jl")
include("cellbuilder.jl")
include("spatialcoherence.jl")

const FILES_TO_KEEP = r"\.dat|\.tif|\.cel|\.cif"

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
    config = YAML.load_file(basename(config_file)) 
    println("Using config file $(basename(config_file)) in $config_file_dir")
    check_config(config)

    run_config(config)

    cd(start_dir)
    cleanup(config)
end

function run_drprobe(config::Dict)
    #Make paths absolute to prevent problems later
    start_dir = abspath(pwd())
    check_config(config)

    run_config(config)

    cd(start_dir)
    cleanup(config)
end

function run_config(config::Dict)
    
    #Copy the input file (.cif/.cel) to the temporary folder
    cp(config["input"], joinpath(config["temporary-folder"], basename(config["input"])))


    if !config["run-celslc"] && config["run-msa"]
        println("CELSLC disabled, looking for slice files")
        link_slice_files(config)
    end
    cd(config["temporary-folder"])
    println("Running CellBuilder")
    config = cellbuilder(config)
    if config["run-celslc"]; celslc(config); end;

    if config["run-msa"]
        make_detector_prm_file(config)
        msa(config)

        println("Applying spatial convolution")
        msa_spatial_convolution(config)
        try 
            make_images(config)
        catch
            println("Could not make images")
        end
 
    else
        println("MSA deactivated, only running CELSLC")
    end
end

function check_config(config::Dict)

    config["temporary-folder"] = config["temporary-folder"] == "" ? ENV["TMPDIR"] : config["temporary-folder"]
    config["temporary-folder"] = abspath(mktempdir(config["temporary-folder"]))
    println("Using temporary folder $(config["temporary-folder"])")

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
    files_to_move = "no-cleanup" ∈ config["debug"] ?  readdir(config["temporary-folder"]) : files_to_move

    for file in files_to_move
        try
            mv(
                joinpath(config["temporary-folder"], file), 
                joinpath(config["output-folder"], file)
              )
        catch 
            continue
        end
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

function open_data_as_matrix(
    config::Dict,
    filename::String
)
    resolution = (config["scan-frame"]["resolution"]["y"], config["scan-frame"]["resolution"]["x"])
    output = Matrix{Float32}(undef, resolution...)
    f = open(filename)
    read!(f, output)
    close(f)
    output
end

function write_matrix_as_data(
    matrix::AbstractMatrix{<:Real},
    filename::String
)
    f = open(filename, "w")
    write(f, matrix)
    close(f)
end
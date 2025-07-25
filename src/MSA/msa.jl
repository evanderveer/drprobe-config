"""
    run_msa(config::Dict)

    Run the Multislice Algorithm (MSA) with the provided configuration.

    # Arguments
    - `config::Dict`: Dictionary containing configuration parameters.

    # Example
    ```julia
    run_msa(config, true)
"""
function msa(config::Dict)

    resolution_from_pixel_size(config)
    println("Number of scan lines: $(config["scan-frame"]["resolution"]["y"])")
    if Threads.nthreads() == 1 
        run_msa_singlethreaded(config)
    else
        run_msa_multithreaded(config)
    end
end

function run_msa_singlethreaded(config::Dict)

    make_msa_prm_file(config)
    msa_command = make_msa_command(config)
    println("Running MSA with command: ")
    println(msa_command)
    "print-commands" ∈ config["debug"] || run(msa_command)
end

"""
    run_msa_multithreaded(config::Dict)

    Run the Multislice Algorithm (MSA) using multiple threads, iterating over scan lines.

    # Arguments
    - `config::Dict`: Dictionary containing configuration parameters.

    # Example
    ```julia
    run_msa_multithread(config, true)
"""
function run_msa_multithreaded(config::Dict)

    if "test-multithreading" ∈ config["debug"] 
        msa_iterate_lines_debug(config)
    else
        msa_iterate_lines(config)
    end

    if "no-stitch" ∈ config["debug"] || "print-commands" ∈ config["debug"]
         return
    end
    println("Stitching lines")
    stitch_lines(config)
end

"""
    make_msa_prm_file(
        config::Dict
        )

    Generate the parameter file for the Multislice Algorithm (MSA) based on the given configuration.

    # Arguments
    - `config::Dict`: Dictionary containing configuration parameters.

    # Example
    ```julia
    make_msa_prm_file(config)
"""
function make_msa_prm_file(config::Dict)

    f = haskey(config, "scan-line") ? open("msa$(config["scan-line"]).prm", "w") : open("msa.prm", "w")
    write(f, "'[Microscope Parameters]'\n")
    config["aperture"]["radius"]
    write(f, "$(config["aperture"]["radius"]),\
              $(config["aperture"]["asymmetry"]["amount"]),\
              $(config["aperture"]["asymmetry"]["angle"]),\
              $(config["aperture"]["smoothness"])\n")

    write(f, "1\n")#Detector inner radius, not used
    write(f, "1\n")#Detector outer radius, not used

    write(f, "1, 'detector_settings.prm'\n")

    write(f, "$(config["high-tension"])\n")
    write(f, "$(config["source-radius"])\n")
    write(f, "$(config["focus-spread"]["value"])\n")
    write(f, "$(config["focus-spread"]["kernel-size"])\n")
    write(f, "$(config["focus-spread"]["kernel-steps"])\n")

    write(f, "24\n")
    for (i, aber) in enumerate(config["aberrations"])
        write(f, "$(i-1) $aber\n")
    end

    write(f, "'[Multislice Parameters]'\n")
    write(f, "$(config["tilt"]["object"]["x"])\n")
    write(f, "$(config["tilt"]["object"]["y"])\n")

    write(f, "$(config["scan-frame"]["offset"]["x"])\n")
    write(f, "$(config["scan-frame"]["offset"]["y"])\n")

    write(f, "$(config["scan-frame"]["size"]["x"])\n")
    write(f, "$(config["scan-frame"]["size"]["y"])\n")

    write(f, "$(config["scan-frame"]["rotation"])\n")

    write(f, "$(config["scan-frame"]["resolution"]["x"])\n")
    write(f, "$(config["scan-frame"]["resolution"]["y"])\n")

    focus_spread = config["convolutions"]["focus-spread"] ? 1 : 0
    write(f, "$focus_spread\n")

    if config["convolutions"]["source-size"]["activate"]
        if config["convolutions"]["source-size"]["profile"] == "gaussian"
            source_size = 1
        else
            source_size = 2
        end
    else
        source_size = 0
    end
    write(f, "$source_size\n")

    write(f, "$(config["repeat-supercell"]["x"])\n")
    write(f, "$(config["repeat-supercell"]["y"])\n")
    write(f, "1\n") #Supercell repeat factor along z, obsolete

    write(f, "'./$(config["output"])'\n")
    num_sli_files = number_of_slice_files(config["output"])
    write(f, "$num_sli_files\n")

    if config["frozen-lattice"]["activate"]
        write(f, "$(config["frozen-lattice"]["number-of-variants"])\n")
        write(f, "$(config["minimum-frozen-lattice-variants"])\n")
    else
        write(f, "1\n1\n")
    end
    
    write(f, "$(readout_period(config["readout-period"], num_sli_files))\n")

    max_number_of_slices = num_sli_files * config["final-thickness"]
    write(f, "$max_number_of_slices\n")

    for _ in 1:config["final-thickness"]
        for i in 1:num_sli_files
            write(f, "$(i-1)\n")
        end
    end

    close(f)
end

function resolution_from_pixel_size(config)
    frame_size_x = config["scan-frame"]["size"]["x"]
    frame_size_y = config["scan-frame"]["size"]["y"]

    pixel_size_x = config["scan-frame"]["pixel-size"]["x"]
    pixel_size_y = config["scan-frame"]["pixel-size"]["y"]

    yres = frame_size_y / pixel_size_y
    xres = frame_size_x / pixel_size_x

    config["scan-frame"]["resolution"] = Dict("x"=>0, "y"=>0)

    config["scan-frame"]["resolution"]["x"] = clamp(round(Int, xres), 1, Inf)
    config["scan-frame"]["resolution"]["y"] = clamp(round(Int, yres), 1, Inf)
end

"""
    readout_period(
        period_dict::Dict, 
        num_sli_files::Int
        )

    Calculate the readout period based on the provided dictionary of readout settings from 
    the configuration file and the number of slice files.

    # Arguments
    - `period_dict::Dict`: Dictionary containing readout period parameters.
    - `num_sli_files::Int`: Number of slice files.

    # Example
    ```julia
    readout_period(period_dict, num_sli_files)
"""
function readout_period(
    period_dict::Dict, 
    num_sli_files::Int
    )
    if period_dict["full-thickness-only"]
        return 0
    elseif period_dict["one-unit-cell"]
        return num_sli_files
    else
        return period_dict["fixed-number"]
    end 
end

"""
    number_of_slice_files(base_name::String)

    Count the number of slice files present.

    # Arguments
    - `base_name::String`: Base name used to form the slice file pattern.

    # Example
    ```julia
    number_of_slice_files("output")
"""
function number_of_slice_files(
    base_name::String
    )
    all_files = readdir()
    slice_file_pattern = Regex(join([".*", base_name, raw"_[0-9]+.sli$"]))
    sum(occursin.(slice_file_pattern, all_files))
end

"""
    make_msa_command(
        config::Dict;
        input::String = ""
        )

    Generate the MSA command based on the given configuration.

    # Arguments
    - `config::Dict`: Dictionary containing configuration parameters.
    - `input::String`: Filename of an input file for source size convolution (optional)

    # Example
    ```julia
    make_msa_command(config, true)
"""
function make_msa_command(
    config::Dict;
    input::String = ""
    )

    if haskey(config, "scan-line")
        prm_file = `-prm msa$(config["scan-line"]).prm`
        base_output = "$(config["output"])_$(config["scan-line"])"
    else
        prm_file = `-prm msa.prm`
        base_output = config["output"]
    end

    out_file = input != "" ? `-out $input-convoluted.dat` : `-out $base_output.dat`

    in_file = input != "" ? `-in $input` : ``

    beam_tilt = `-tx $(config["tilt"]["beam"]["x"]) -ty $(config["tilt"]["beam"]["y"])`

    #TODO: Figure out absorption and Debye-Waller factors

    ctem = config["ctem"] ? `/ctem` : ``
    text_output = config["text-output"] ? `/txtout` : ``
    three_d_output = config["3d-output"] ? `/3dout` : ``
    if config["wavefunction-output"]["activate"] 
        wave = config["wavefunction-output"]["average-over-frozen-lattice-variants"] ? `/avwave` : `/wave`
        wave = config["wavefunction-output"]["output-in-fourier-space"] ? `/avwaveft` : wave
    else
        wave = ``
    end
    detimg = config["detector-function-ouput"] ? `/detimg` : ``
    lapro = config["use-large-angle-propagators"] ? `/lapro` : ``
    dftest = config["use-FFTW-ESTIMATE"] ? `/dftest` : ``
    verbose = config["verbose"] ? `/verbose` : ``
    debug = config["msa-debug"] ? `/debug` : ``
    silent = config["silent"] ? `/silent` : ``
    load_ondemand = config["load-slices-on-demand"] ? `/slod` : ``

    `$MSA $prm_file $out_file $in_file $beam_tilt $ctem $text_output $three_d_output $wave $detimg $lapro $dftest $verbose $debug $silent $load_ondemand`
end

"""
    make_detector_prm_file(config::Dict)

    Generate the parameter file for the detector settings based on the given configuration.

    # Arguments
    - `config::Dict`: Dictionary containing configuration parameters.

    # Example
    ```julia
    make_detector_prm_file(config)
"""
function make_detector_prm_file(
    config::Dict
    )
    f = open("detector_settings.prm", "w")
    
    write(f, "'[Detector Parameters]'\n")
    write(f, "2016021801\n")
    write(f, string(config["detector"]["number-of-detectors"])*'\n')

    for detector in keys(config["detector"]["detectors"])
        detector_dict = config["detector"]["detectors"][detector]
        sensitivity_file = detector_dict["use-sensitivity"] ? detector_dict["sensitivity-file"] : ""
        detector_string = "$(detector_dict["inner-radius"]),$(detector_dict["outer-radius"]),\
                           $(detector_dict["azimuth-start"]),$(detector_dict["azimuth-end"]),\
                           $(detector_dict["x-center"]),$(detector_dict["y-center"]),\
                           '$detector','$sensitivity_file'\n"
        write(f, detector_string)
    end

    close(f)
end

function msa_spatial_convolution(config::Dict)

    if !config["convolutions"]["source-size"]["activate"]
        return
    end

    all_files = readdir()
    output_regex = Regex(join([config["output"], raw"_.+_sl\d+\.dat"]))
    output_files = all_files[occursin.(output_regex, all_files)]

    for file in output_files

        println("Applying source size convolution to file $file")
        data = open_data_as_matrix(config, file)
        data_conv = apply_spatial_coherence(config, data)
        filename_conv = splitext(file)[1] * "_convoluted.dat"
        write_matrix_as_data(data_conv, filename_conv)
    end
end

function msa_iterate_lines(config)
    Threads.@threads for line_number in 1:config["scan-frame"]["resolution"]["y"]
        msa_line(deepcopy(config), line_number)
    end
end

function msa_iterate_lines_debug(config)
    println("Running multithreading debug mode")
    for line_number in 1:config["scan-frame"]["resolution"]["y"]
        msa_line(deepcopy(config), line_number)
    end
end

"""
    msa_line(config::Dict, line_number::Int)

    Run the Multislice Algorithm (MSA) for a specific scan line.

    # Arguments
    - `config::Dict`: Dictionary containing configuration parameters.
    - `line_number::Int`: Line number for the scan.

    # Example
    ```julia
    msa_line(config, true, 1)
"""
function msa_line(
    config::Dict, 
    line_number::Int
    )

    line_height = config["scan-frame"]["size"]["y"]/config["scan-frame"]["resolution"]["y"]
    config["scan-frame"]["offset"]["y"] += (line_number-1)*line_height
    config["scan-frame"]["size"]["y"] = line_height
    config["scan-frame"]["resolution"]["y"] = 1
    config["scan-line"] = line_number
    config["silent"] = true #Prevents garbled output
    
    run_msa_singlethreaded(config)
end

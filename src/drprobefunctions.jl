using YAML 

function build_celslc_command(config)
    input = config["input"][end] == 'f' ? `-cif $(config["input"])` : `-cel $(config["input"])` 
    output = `-slc $(config["output"])`

    pot_files = config["export-potentials"]["as-pot-files"] ? `-pot` : ``
    sli_files = config["export-potentials"]["as-sli-files"] ? `-pps` : `` 

    if config["projection"]["activate"]
        axis = join(string.(config["projection"]["axis"]), ',')
        vertical = join(string.(config["projection"]["vertical"]), ',')
        supercell = join(string.(config["projection"]["supercell-size"]), ',')
        prj_params = join([axis, vertical, supercell], ',')
        proj = `-prj $prj_params`
    else
        proj = ``
    end

    if config["atom-shifts"]["activate"]
        tla_params = join(config["atom-shifts"]["shifts"], ',')
        tla = `-tla $tla_params`
    else
        tla = ``
    end

    ht = `-ht $(config["high-tension"])`

    nxny = `-nx $(config["samples"]["x-samples"]) -ny $(config["samples"]["y-samples"])`

    if config["samples"]["use-equidistant"]
        nz = `-nz $(config["samples"]["z-samples"])`
    else
        nz = ``
    end

    rev = config["samples"]["reverse"] ? `-rev` : ``

    if config["frozen-lattice"]["activate"]
        nv = `-fl -nv $(config["frozen-lattice"]["number-of-variants"])`
    else
        nv = ``
    end

    if config["single-slice"]["activate"]
        ssc = `-ssc $(config["single-slice"]["slice-number"])`
    else
        ssc = ``
    end

    dwf = config["debye-waller"] ? `-dwf` : ``

    if config["absorption"]["activate"]
        abs = config["absorption"]["use-wk"] ? `-abs` : `-abf $(config["absorption"]["value"])`
    else
        abs = ``
    end

    `celslc $input $output $pot_files $sli_files $proj $tla $ht $nxny $nz $rev $nv $ssc $dwf $abs`
end

function make_msa_prm_file(config)
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
    for (i,aber) in enumerate(config["aberrations"])
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

function readout_period(period_dict, num_sli_files)
    if period_dict["full-thickness-only"]
        return 0
    elseif period_dict["one-unit-cell"]
        return num_sli_files
    else
        return period_dict["fixed-number"]
    end 
end

function number_of_slice_files(base_name)
    all_files = readdir()
    regexs = ["^", base_name, raw"_\d\d\d.sli$"]
    sum(occursin.(Regex(join(regexs, "")), all_files))
end

function make_msa_command(config, input::Bool = false)

    if haskey(config, "scan-line")
        prm_file = `-prm msa$(config["scan-line"]).prm`
        base_output = "$(config["output"])_$(config["scan-line"])"
    else
        prm_file = `-prm msa.prm`
        base_output = config["output"]
    end

    out_file = input ? `-out $base_output-convoluted.dat` : `-out $base_output.dat`

    in_file = input ? `-in $base_output.dat` : ``

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
    debug = config["debug"] ? `/debug` : ``
    silent = config["silent"] ? `/silent` : ``

    `msa $prm_file $out_file $in_file $beam_tilt $ctem $text_output $three_d_output $wave $detimg $lapro $dftest $verbose $debug $silent`
end

function make_detector_prm_file(config)
    f = open("detector_settings.prm", "w")
    
    write(f, "'[Detector Parameters]'\n")
    write(f, "2016021801\n")
    write(f, string(config["detector"]["number-of-detectors"])*'\n')

    for detector in keys(config["detector"]["detectors"])
        detector_dict = config["detector"]["detectors"][detector]
        sensitivity_file = detector_dict["use-sensitivity"] ? detector_dict["sensitivity-file"] : ""
        detector_string = "$(detector_dict["inner-radius"]),$(detector_dict["outer-radius"]),\
                           $(detector_dict["azimuth-start"]),$(detector_dict["azimuth-start"]),\
                           $(detector_dict["x-center"]),$(detector_dict["y-center"]),\
                           '$detector','$sensitivity_file'\n"
        write(f, detector_string)
    end

    close(f)
end

function cleanup(config, temp_folder, output_folder)
    println("Cleaning up the mess")
    output_folder = mkdir(string(now()))
    all_files = readdir(temp_folder)
    files_to_move = all_files[occursin.(r"\.dat|\.sli", all_files)]
    for file in files_to_move
        mv(joinpath(temp_folder, file), joinpath(output_folder, file))
    end
    rm(temp_folder, recursive=true)
end

function run_msa(config, debug)
    
    make_msa_prm_file(config)
    msa_command = make_msa_command(config)
    println("Running MSA with command: ")
    println(msa_command)
    debug || run(msa_command)
end

function run_drprobe(
    config_file; 
    debug=false, 
    output_folder=pwd(), 
    no_cleanup=false
    )
    start_dir = pwd()
    cd(dirname(config_file))
    config = YAML.load_file(basename(config_file))
    temp_dir = mktempdir(pwd())
    cp(config["input"], temp_dir*"/"*config["input"])
    cd(temp_dir)

    celslc_command = build_celslc_command(config)
    println("Running CELSLC with command: ")
    println(celslc_command)
    debug || run(celslc_command)
    
    if config["run-msa"]
        make_detector_prm_file(config)
        run_msa(config, debug)
    else
        println("MSA deactivated, only running CELSLC")
    end

    cd(start_dir)
    
    debug || no_cleanup || cleanup(config, joinpath(dirname(config_file), temp_dir), output_folder);
    println(" ")
end

function run_drprobe_multithreaded(
    config_file; 
    debug=false, 
    output_folder=pwd(), 
    no_cleanup=false
    )
    start_dir = pwd()
    cd(dirname(config_file))
    config = YAML.load_file(basename(config_file))
    temp_dir = mktempdir(pwd())
    cp(config["input"], temp_dir*"/"*config["input"])
    cd(temp_dir)

    celslc_command = build_celslc_command(config)
    println("Running CELSLC with command: ")
    println(celslc_command)
    debug || run(celslc_command)
    
    if config["run-msa"]
        make_detector_prm_file(config)
        msa_multithread(config, debug)
        debug || stitch_lines(config)
    else
        println("MSA deactivated, only running CELSLC")
    end
    
    cd(start_dir)

    debug || no_cleanup || cleanup(config, temp_dir, output_folder);
    println(" ")
end

function msa_multithread(
    config,
    debug
    )
    Threads.@threads for line_number in 1:config["scan-frame"]["resolution"]["y"]
        msa_line(deepcopy(config), debug, line_number)
    end
end

function msa_line(config, debug, line_number)

    line_height = config["scan-frame"]["size"]["y"]/config["scan-frame"]["resolution"]["y"]
    config["scan-frame"]["offset"]["y"] = (line_number-1)*line_height
    config["scan-frame"]["size"]["y"] = line_height
    config["scan-frame"]["resolution"]["y"] = 1
    config["scan-line"] = line_number
    config["silent"] = true #Prevents garbled output
    
    run_msa(config, debug)
end

function stitch_lines(config)
    filelist = readdir()
    match_list = match.(r"test_([0-9])_(.+)_sl([0-9]+)\.dat", filelist)
    match_list = match_list[isnotnothing.(match_list)]

    detectors = Set(getindex.(match_list, 2))
    slice_nums = Set(getindex.(match_list, 3))

    for detector in detectors
        for slice in slice_nums
            stitch_image(config, detector, slice)
        end
    end
end

isnotnothing(value) = isnothing(value) |>  !

function stitch_image(config, detector, slice)
    output_size = (config["scan-frame"]["resolution"]["y"], 
                   config["scan-frame"]["resolution"]["x"])

    output_data = Matrix{Float32}(undef, output_size...)
    for line in 1:config["scan-frame"]["resolution"]["y"]
        line_filename = "$(config["output"])_$(line)_$(detector)_sl$(slice).dat"
        line_data = similar(output_data[:, line])
        read!(line_filename, line_data)
        output_data[:, line] = line_data
        rm(line_filename)
    end
    out_filename = "$(config["output"])_$(detector)_sl$(slice).dat"

    write(out_filename, output_data)
end
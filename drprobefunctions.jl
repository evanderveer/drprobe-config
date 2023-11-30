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
        prj_params = join([axis,vertical,supercell], ',')
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
    f = open("msa.prm", "w")
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
    
    close(f)
end

function make_msa_command(config)

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

function run_drprobe(config)
    temp_dir = mktempdir(pwd())
    cp(config["input"], temp_dir*"/"*config["input"])
    cd(temp_dir)

    celslc_command = build_celslc_command(config)
    make_detector_prm_file(config)
    make_msa_prm_file(config)
    msa_command = make_msa_command(config)
    run(celslc_command)

    cd("..")
end
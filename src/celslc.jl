"""
build_celslc_command(
    config::Dict
    )

Build the CELSLC command from the configuration supplied.

# Arguments
- `config::Dict`: A dictionary containing configuration parameters.

# Returns
- CELSLC command
"""
function build_celslc_command(
    config::Dict
    )
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

"""
function run_celslc(
    config::Dict,
    debug::Bool
    )

Build A CELSLC command from the configuration supplied, then run CELSLC.

# Arguments
- `config::Dict`: A dictionary containing configuration parameters.
- `debug::Bool`: Flag for turning on debug mode

# Returns
- Nothing, but writes phase gratings/projected potentials to disk.
"""
function run_celslc(
    config::Dict,
    debug::Bool
    )
    celslc_command = build_celslc_command(config)
    println("Running CELSLC with command: ")
    println(celslc_command)
    debug || run(celslc_command)
end
function cellbuilder(config)

    if !occursin(r"\.cif$|\.cel$", config["input"])
        error("Incompatible input file. Must be .cel or .cif.")
    end
    if occursin(r"\.cif$", config["input"])
        temp_cel_file = config["output"] * "_buildcell.cel"
        run(`$BUILDCELL --cif=$(config["input"]) --output=$temp_cel_file`)
        config["input"] = temp_cel_file
    end
    temp_cel_file = config["output"] * "_cellmuncher.cel"
    run(`$CELLMUNCHER --input=$(config["input"]) --output=$temp_cel_file --cif`)
    config["input"] = temp_cel_file

    cell_parameters, basis = load_cell(config["input"])
    CoBmatrix = find_orthogonal_cell(
        config["zone-axis"], 
        cell_parameters,
        tolerance=config["tolerance"], 
        maximum_iterations=config["iteration-limit"], 
        maximum_index=config["index-limit"]
        )
    
    new_basis = transform_basis(basis, CoBmatrix)
    new_cp = new_cell_parameters(cell_parameters, CoBmatrix)

    temp_cel_file = config["output"] * "_cellbuilder.cel"
    save_cell(temp_cel_file, new_cp, new_basis)

    config["input"] = temp_cel_file
    temp_cel_file = config["output"] * "_final.cel"

    run(`$CELLMUNCHER --input=$(config["input"]) --output=$temp_cel_file --cif`)
    config["input"] = temp_cel_file

    if config["scan-frame"]["size"]["auto-detect"]
        ysize, xsize = get_size_from_cel(temp_cel_file)
        config["scan-frame"]["size"]["y"] = ysize
        config["scan-frame"]["size"]["x"] = xsize
    end

    config #Return the config 
end

function get_size_from_cel(cel_file)
    f = open(cel_file)
    
    readline(f) #Ignore the first line of the file
    parameters = [s for s in split(readline(f), " ") if (s != "" && s != " ")]
    close(f)
    xsize, ysize = parse.(Float64, parameters[2:3])
    return ysize, xsize
end
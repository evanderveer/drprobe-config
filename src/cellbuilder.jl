function cellbuilder(config)

    if !occursin(r"\.cif$|\.cel$", config["input"])
        error("Incompatible input file. Must be .cel or .cif.")
    end
    if occursin(r"\.cif$", config["input"])
        temp_cel_file = config["output"] * "_buildcell.cel"
        run(`buildcell --cif=$(config["input"]) --output=$temp_cel_file`)
        config["input"] = temp_cel_file
    end
    temp_cel_file = config["output"] * "_cellmuncher.cel"
    run(`cellmuncher --input=$(config["input"]) --output=$temp_cel_file --cif`)
    config["input"] = temp_cel_file

    cell_parameters, basis = load_cell(config["input"])
    CoBmatrix = find_orthogonal_cell(
        config["zone-axis"], 
        cell_parameters,
        tolerance=config["tolerance"], 
        max_iterations=config["iteration-limit"], 
        max_index=config["index-limit"]
        )
    
    new_basis = transform_basis(basis, CoBmatrix)
    new_cp = new_cell_parameters(cell_parameters, CoBmatrix)

    temp_cel_file = config["output"] * "_cellbuilder.cel"
    save_cell(temp_cel_file, new_cp, new_basis)

    run(`cellmuncher --input=$(config["input"]) --output=$temp_cel_file --cif`)
    config["input"] = temp_cel_file
    config #Return the config 
end
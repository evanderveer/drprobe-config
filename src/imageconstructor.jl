function make_images(
        config::Dict
    )
    all_files = readdir()
    files_to_convert = all_files[occursin.(r"sl\d+\.dat", all_files)]
    
    for file in files_to_convert
        image = make_image(file, 
                           (config["scan-frame"]["resolution"]["y"], 
                            config["scan-frame"]["resolution"]["x"]), 
                           true)
        output_filename = splitext(file)[1]*".tif"
        save(output_filename, image)
    end
end

function make_image(
    name::String,
    size::Tuple{<:Int, <:Int}, 
    stretch::Bool
    )
    dat = Matrix{Float32}(undef, size...)
    read!(name, dat)
    if stretch
        dat .-= minimum(dat)
        dat ./= maximum(dat)
    end
    Gray.(dat)
end
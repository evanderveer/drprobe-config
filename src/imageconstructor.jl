function make_images(
        config::Dict
    )
    all_files = readdir()
    files_to_convert = all_files[occursin.(r"sl\d+\.dat", all_files)]

    println("Converting data files to images: ")
    for f in files_to_convert; println(f); end;
    
    for file in files_to_convert
        image = make_image(file, 
                           (config["scan-frame"]["resolution"]["y"], 
                            config["scan-frame"]["resolution"]["x"]), 
                            true)
        output_filename = splitext(file)[1]*".tif"
        save(output_filename, image)
    end
    if check_idpc(config)
        make_idpc_images(config)
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

#Check if iDPC image can be made
function check_idpc(config::Dict)
    detectors = keys(config["detector"]["detectors"])
    for segment in ["A", "B", "C", "D"]
        if !("DF4-"*segment ∈ detectors)
            return false
        end
    end
    return true
end

function image_filename(config, slice_num, segment, conv)
    conv_string = conv ? ".dat-convoluted" : ""
    "$(config["output"])_DF4-$(segment)_sl$slice_num$conv_string.tif"
end

function make_idpc_images(config::Dict)
    all_files = readdir()
    matches = match.(r"DF4-[A-D]_sl(\d+)(\.dat-convoluted)?\.tif", all_files)
    ipdc_to_convert = matches[isnotnothing.(matches)]

    slice_nums = Set(getindex.(ipdc_to_convert, 1))
    is_convoluted = any(isnothing.(getindex.(ipdc_to_convert, 2)))

    for slice_num in slice_nums 
        make_idpc_slice(config, slice_num, false)
        if is_convoluted
            make_idpc_slice(config, slice_num, true)
        end
    end
end

function make_idpc_slice(config, slice_num, convoluted)
    image_size = (config["scan-frame"]["resolution"]["y"], 
                  config["scan-frame"]["resolution"]["x"])
    images = []
    for segment in ["A", "B", "C", "D"]
        image = Matrix{Float32}(undef, image_size...)
        fn = image_filename(config, slice_num, segment, convoluted)
        println("Processing $fn")
        f = open(fn)
        read!(f, image)
        close(f)
        push!(images, image)
    end
    idpc,ddpc = dpc(images...)
    save(image_filename(config, slice_num, "iDPC", convoluted), idpc)
    save(image_filename(config, slice_num, "dDPC", convoluted), ddpc)

end

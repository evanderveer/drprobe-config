function make_images(
        config::Dict
    )
    all_files = readdir()
    files_to_convert = all_files[occursin.(r".dat$", all_files)]

    println("Converting data files to images: ")
    for f in files_to_convert; println(f); end;
    
    for file in files_to_convert
        data = open_data_as_matrix(config, file)
        image = Gray.(normalize_image(data))

        output_filename = splitext(file)[1]*".tif"
        save(output_filename, image)
    end
    if check_idpc(config)
        make_idpc_images(config)
    end
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
    conv_string = conv ? "_convoluted" : ""
    "$(config["output"])_DF4-$(segment)_sl$slice_num$conv_string"
end

function make_idpc_images(config::Dict)
    all_files = readdir()
    matches = match.(r"DF4-[A-D]_sl(\d+)(_convoluted)?\.tif", all_files)
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
    images = []
    for segment in ["A", "B", "C", "D"]
        fn = image_filename(config, slice_num, segment, convoluted) * ".dat"
        println("Processing $fn")
        image_data = open_data_as_matrix(config, fn)
        push!(images, image_data)
    end

    #Not sure why the order needs to be this
    idpc, ddpc = dpc(images..., order=[1,4,3,2])
    save(image_filename(config, slice_num, "iDPC", convoluted) * ".tif", clamp01nan.(idpc))
    save(image_filename(config, slice_num, "dDPC", convoluted) * ".tif", clamp01nan.(ddpc))

end

function normalize_image(
    image::Matrix{<:Real}
    )
    image_shifted = image .- minimum(image)
    image_shifted ./ maximum(image_shifted)
end

function tile_image(
    image::AbstractMatrix{T},
    factor::Tuple{<:Int, <:Int}
) where T

    image_size = size(image)
    new_image = Matrix{T}(undef, (image_size .* factor)...)
    for i in 1:factor[1]
    for j in 1:factor[2]
        im_start = image_size .* (i-1, j-1)
        im_end = image_size .* (i, j)

        new_image[im_start[1]+1:im_end[1], im_start[2]+1:im_end[2]] = image
    end
    end
    new_image
end
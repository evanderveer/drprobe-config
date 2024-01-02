"""
isnotnothing(value) = isnothing(value) |> !

Check if the given value is not nothing.

# Arguments
- `value`: The value to check.

# Returns
`true` if the value is not nothing, `false` otherwise.
"""
isnotnothing(value) = isnothing(value) |> !

"""
data_file_name_pattern(base_name::String)

Generate a regular expression pattern for matching data file names.

# Arguments
- `base_name::String`: The base name for the data file.

# Returns
A regular expression pattern for matching data file names.
"""
function data_file_name_pattern(base_name::String)
    Regex(base_name * raw"_([0-9])_(.+)_sl([0-9]+).dat")
end

"""
stitch_lines(config::Dict)

Stitch together images from multiple data files (one for each scan line)
based on the provided configuration.

# Arguments
- `config::Dict`: A dictionary containing configuration parameters.

# Returns
Nothing, but writes stitched images to output files.
"""
function stitch_lines(
    config::Dict
    )

    output_files = []
    for file in readdir()
        file_match = match(data_file_name_pattern(config["output"]), file)
        if isnotnothing(file_match)
            push!(output_files, file_match)
        end
    end

    detectors = Set(getindex.(output_files, 2))
    slice_nums = Set(getindex.(output_files, 3))

    if "test-multithreading" in config["debug"]
        iterate_detectors_slices_debug(config, detectors, slice_nums)
    else
        iterate_detectors_slices(config, detectors, slice_nums)
    end
end

function iterate_detectors_slices(config, detectors, slice_nums)
    Threads.@threads for detector in collect(detectors)
        for slice in collect(slice_nums)
            stitch_image(config, detector, slice)
        end
    end
end

function iterate_detectors_slices_debug(config, detectors, slice_nums)
    for detector in collect(detectors)
        for slice in collect(slice_nums)
            stitch_image(config, detector, slice)
        end
    end
end


"""
stitch_image(
    config::Dict, 
    detector::String, 
    slice::String
            )

Stitch together an image based on the provided configuration, detector, and slice.

# Arguments
- `config::Dict`: A dictionary containing configuration parameters.
- `detector::String`: The detector name.
- `slice::String`: The slice number.

# Returns
Nothing, but writes the stitched image to an output file.
"""
function stitch_image(
    config::Dict,
    detector,
    slice
    )
    output_size = (config["scan-frame"]["resolution"]["y"], 
                   config["scan-frame"]["resolution"]["x"])
    output_data = Matrix{Float32}(undef, output_size...)

    fill_image!(output_data, config, detector, slice)
 
    out_filename = "$(config["output"])_$(detector)_sl$(slice).dat"

    write(out_filename, output_data)
end

"""
fill_image!(
    output_data::AbstractMatrix{<:Real}, 
    config::Dict,
    detector::String,
    slice::String)

Fill the output image matrix with data from individual lines.

# Arguments
- `output_data::AbstractMatrix{<:Real}`: The matrix to fill with image data.
- `config::Dict`: A dictionary containing configuration parameters.
- `detector::String`: The detector name.
- `slice::String`: The slice number.

# Returns
Nothing, but modifies the `output_data` matrix in place.
"""
function fill_image!(
    output_data::AbstractMatrix{<:Real},
    config::Dict,
    detector,
    slice
)
    for line in 1:config["scan-frame"]["resolution"]["y"]
        line_filename = "$(config["output"])_$(line)_$(detector)_sl$(slice).dat"
        
        f = open(line_filename)
        output_data[line, :] .= readeach(f, Float32)
        close(f)
        rm(line_filename)
    end
end

# Based on this repository:
# https://github.com/neago/lecroy-reader

using Printf

# Define the structure of the Lecroy .trc file format
const lecroy_format = Dict(
    "descriptor_name" => (0, "16s"),
    "template_name" => (16, "16s"),
    "comm_type" => (32, "h"),
    "comm_order" => (34, "h"),
    "wave_descriptor" => (36, "i"),
    "user_text" => (40, "i"),
    "res_desc1" => (44, "i"),
    "trig_time_array" => (48, "i"),
    "ris_time_array" => (52, "i"),
    "res_array1" => (56, "i"),
    "wave_array1" => (60, "i"),
    "wave_array2" => (64, "i"),
    "res_array2" => (68, "i"),
    "res_array3" => (72, "i"),
    "instrument_name" => (76, "16s"),
    "instrument_number" => (92, "i"),
    "trace_label" => (96, "16s"),
    "reserved1" => (112, "h"),
    "reserved2" => (114, "h"),
    "wave_array_count" => (116, "i"),
    "points_per_screen" => (120, "i"),
    "first_valid_point" => (124, "i"),
    "last_valid_point" => (128, "i"),
    "first_point" => (132, "i"),
    "sparsing_factor" => (136, "i"),
    "segment_index" => (140, "i"),
    "subarray_count" => (144, "i"),
    "sweeps_per_acq" => (148, "i"),
    "points_per_pair" => (152, "h"),
    "pair_offset" => (154, "h"),
    "vertical_gain" => (156, "f"),
    "vertical_offset" => (160, "f"),
    "max_value" => (164, "f"),
    "min_value" => (168, "f"),
    "nominal_bits" => (172, "h"),
    "nom_subarray_count" => (174, "h"),
    "horiz_interval" => (176, "f"),
    "horiz_offset" => (180, "d"),
    "pixel_offset" => (188, "d"),
    "vert_unit" => (196, "48s"),
    "horiz_unit" => (244, "48s"),
    "horiz_uncertainty" => (292, "f"),
    "trigger_time" => (296, "dbbbbhh"),
    "acq_duration" => (312, "f"),
    "record_type" => (316, "h"),
    "processing_done" => (318, "h"),
    "reserved5" => (320, "h"),
    "ris_sweeps" => (322, "h"),
    "time_base" => (324, "h"),
    "vert_coupling" => (326, "h"),
    "probe_att" => (328, "f"),
    "fixed_vert_gain" => (332, "h"),
    "bandwidth_limit" => (334, "h"),
    "vertical_vernier" => (336, "f"),
    "acq_vert_offset" => (340, "f"),
    "wave_source" => (344, "h")
)


# Define the mapping of formats to Julia types and their sizes
const format_to_type = Dict(
    "16s" => (String, 16),
    "48s" => (String, 48),
    "h" => (Int16, 2),
    "i" => (Int32, 4),
    "f" => (Float32, 4),
    "d" => (Float64, 8),
    "dbbbbhh" => (Tuple{Float64, UInt8, UInt8, UInt8, UInt8, Int16, Int16}, 16)
)

function read_trc_file(filename::String; readdata::Bool=true, scale::Bool=true)
    metadata = Dict{String, Any}()
    raw = read(filename)
    raw_ = copy(raw)
    startpos = first(findfirst("WAVEDESC", String(raw_)))

    for (field, (offset, fmt)) in lecroy_format
        fieldlen = format_to_type[fmt][2]
        metadata[field] = raw[(startpos + offset) : (startpos + offset + fieldlen - 1)]
    end

    comm_order = reinterpret(Int16, metadata["comm_order"])
    boc = comm_order == 0 ? '>' : '<'

    for (field, rawval) in metadata
        fmt = lecroy_format[field][2]
        if boc == '>' && field != "comm_order" && fmt != "16s" && fmt != "48s"
            rawval = reverse(rawval)
        end
        if fmt == "16s" || fmt == "48s"
            val = strip(String(rawval), '\x00')
        elseif fmt == "h"
            val = reinterpret(Int16, rawval)[1]
        elseif fmt == "i"
            val = reinterpret(Int32, rawval)[1]
        elseif fmt == "f"
            # println(rawval)
            val = reinterpret(Float32, rawval)[1]
        elseif fmt == "d"
            val = reinterpret(Float64, rawval)[1]
        elseif fmt == "dbbbbhh"
            val = reinterpret(Float64, rawval[1:8])[1]
            val = [val; reinterpret(UInt8, rawval[9:14])]
        else
            val = rawval
        end
        metadata[field] = val
    end

    if readdata
        # Read the trigger times into an array
        trigtimes_startpos = startpos + metadata["wave_descriptor"] + metadata["user_text"]
        trigtimes = reinterpret(reshape, Float64, reshape(raw[trigtimes_startpos:Int(metadata["trig_time_array"])÷8], 8, :))
        
        # Number format for data. COMM_TYPE: 0 -> byte, 1 -> short
        number_type = metadata["comm_type"] == 0 ? Int8 : Int16

        # Read the binary data into an array
        data_startpos = trigtimes_startpos + metadata["trig_time_array"]
        data = reinterpret(reshape, number_type, reshape(raw[data_startpos:data_startpos+sizeof(number_type)*(metadata["wave_array_count"])-1], sizeof(number_type), :))

        # Scale, offset and reshaping into segments
        if scale
            data = data .* metadata["vertical_gain"] .+ metadata["vertical_offset"]
        end

        if metadata["subarray_count"] > 1
            data = reshape(data, metadata["subarray_count"], :)
        end

        return metadata, trigtimes, data
    else
        return metadata
    end
end

function pretty_metadata(metadata::Dict{String, Any})
    for (key, value) in metadata
        println(rpad(key, 20), value)
    end
end

function get_time_array(metadata::Dict{String, Any})
    return metadata["horiz_interval"] .* (0:metadata["wave_array_count"]-1) .+ metadata["horiz_offset"], 1 / metadata["horiz_interval"]
end
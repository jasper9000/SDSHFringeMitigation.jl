using Adapt

include("LeCroyTRC.jl")

@enum Backend cpu cuda

function load_data(filename::String, backend::Backend=cpu, float_type::Type=nothing)
    metadata, trigtimes, data = read_trc_file(filename)
    # pretty_metadata(metadata)
    times, f_sampling = get_time_array(metadata)
    # println("Sampling frequency: ", f_sampling)

    # convert to float type of the data
    T = typeof(data[1])

    # @warn println("MANUAL OVERRIDE TO FLOAT64")
    # T = Float64

    # # give a warning if the data is not float32
    # if T != Float32
    #     @warn println("Data is not Float32!!!")
    # end

    data = convert(Array{T}, data)
    f_sampling = convert(T, f_sampling)
    times = convert(Array{T}, times)

    # convert to CUDA
    if BACKEND == cuda
        data = adapt(CuArray{T}, data);
        times = adapt(CuArray{T}, times);
    end
    if !isnothing(float_type)
        data = change_float_precision(data, float_type)
        times = change_float_precision(times, float_type)
        f_sampling = float_type(f_sampling)
        T = float_type
    end
    return times, data, f_sampling, T
end

function change_float_precision(arr::AbstractArray, T::Type)
    ArrWrapper = typeof(arr).name.wrapper
    return adapt(ArrWrapper{T}, arr)
end

function adapt_array(arr::AbstractArray{T}, BACKEND::Backend) where T<:Real
    if BACKEND == cpu
        return arr
    elseif BACKEND == cuda
        return adapt(CuArray{T}, arr)
    end
end
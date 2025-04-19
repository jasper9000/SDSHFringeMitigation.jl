module SDSHFringeMitigation

    using LinearAlgebra

    include("est_exp_parameters.jl")

    include("KRR/kernels.jl")
    include("KRR/krr.jl")
    include("KRR/training_sample_selection.jl")
    include("KRR/cross_validation.jl")

    include("data_handling/load_file.jl")
    include("data_handling/LeCroyTRC.jl")


    include("phase_estimation.jl")
    include("welch.jl")
    include("beta_linewidth.jl")
    include("gen_custom_noise.jl")
    include("spectrum_equalization.jl")

    # plotting
    include("plotting/lttb.jl")
    include("plotting/misc.jl")

    println("SDSHFringeMitigation.jl loaded")
end 

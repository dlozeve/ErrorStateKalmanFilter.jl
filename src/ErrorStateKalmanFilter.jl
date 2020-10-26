module ErrorStateKalmanFilter

using LinearAlgebra
using Quaternions

export
    NominalState,
    ErrorState,
    Measurement,
    Noise,
    skew_matrix,
    predict

include("types.jl")
include("ESKF.jl")

end

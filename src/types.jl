import Base.vec

struct NominalState{T<:Real}
    p::Vector{T}
    v::Vector{T}
    q::Quaternion{T}
    ab::Vector{T}
    ωb::Vector{T}
    g::Vector{T}
end

function NominalState(x::AbstractVector)
    NominalState(
        x[1:3],   # position
        x[4:6],   # velocity
        Quaternion(x[7], x[8], x[9], x[10]), # attitude
        x[11:13], # accelerometer bias
        x[14:16], # gyrometer bias
        x[17:19]  # gravity vector
    )
end

function vec(nominal::NominalState)::AbstractVector
    [
        nominal.p
        nominal.v
        real(nominal.q)
        Quaternions.imag(nominal.q)
        nominal.ab
        nominal.ωb
        nominal.g
    ]
end

struct ErrorState{T<:Real}
    δp::Vector{T}
    δv::Vector{T}
    δθ::Vector{T}
    δab::Vector{T}
    δωb::Vector{T}
    δg::Vector{T}
end

function ErrorState(x::AbstractVector)
    ErrorState(
        x[1:3],   # position error
        x[4:6],   # velocity error
        x[7:9],   # angles error
        x[10:12], # accelerometer bias error
        x[13:15], # gyrometer bias error
        x[16:18]  # gravity error
    )
end

function vec(error::ErrorState)::AbstractVector
    [
        error.δp
        error.δv
        error.δθ
        error.δab
        error.δωb
        error.δg
    ]
end

struct Measurement{T<:Real}
    am::Vector{T}
    ωm::Vector{T}
end

struct Noise{T<:Real}
    σan::T # velocity noise Vi (m^2/s^2) 
    σωn::T # orientation noise Θi (rad^2)
    σaw::T # accelero bias noise Ai (m^2/s^4)
    σωw::T # gyro bias noise Ωi (rad^2/s^2)
end

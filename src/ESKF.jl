function skew_matrix(a::Vector{T})::Matrix{T} where T<:Real
    m = [  0    -a[3]  a[2]
           a[3]  0    -a[1]
          -a[2]  a[1]  0    ]
    return m
end

function predict(nominal::NominalState{T}, error::ErrorState{T}, error_cov::AbstractMatrix{T}, measurement::Measurement{T}, noise::Noise{T}, dt::T) where T<:Real
    R = rotationmatrix(nominal.q)

    # Nominal state kinematics
    new_nominal = NominalState(
        nominal.p + nominal.v * dt + 0.5 * (R * (measurement.am - nominal.ab) + nominal.g) * dt^2,
        nominal.v + (R * (measurement.am - nominal.ab) + nominal.g) * dt^2,
        nominal.q * qrotation((measurement.ωm - nominal.ωb) * dt),
        nominal.ab,
        nominal.ωb,
        nominal.g
    )

    # Jacobian of the error state kinematics
    Fx = zeros(6 * 3, 6 * 3)
    Fx[1:3, 1:3]     = Matrix{T}(I, 3, 3)
    Fx[1:3, 4:6]     = dt * Matrix{T}(I, 3, 3)
    Fx[4:6, 4:6]     = Matrix{T}(I, 3, 3)
    Fx[4:6, 7:9]     = -R * skew_matrix(measurement.am - nominal.ab) * dt
    Fx[4:6, 10:12]   = -R * dt
    Fx[4:6, 16:18]   = dt * Matrix{T}(I, 3, 3)
    Fx[7:9, 7:9]     = rotationmatrix(qrotation(measurement.ωm - nominal.ωb) * dt)'
    Fx[7:9, 13:15]   = -dt * Matrix{T}(I, 3, 3)
    Fx[10:18, 10:18] = Matrix{T}(I, 9, 9)

    # Error-state kinematics
    # new_error = ErrorState(Fx * vec(error))
    # We don't need to compute this, as the mean of the error is
    # always zero, so this equation always returns zero.

    # Jacobian of the random impulses kinematics
    Fi = zeros(6 * 3, 4 * 3)
    Fi[4:15, :] = Matrix{T}(I, 12, 12)

    # Covariance matrices of the random impulses
    Qi = zeros(4 * 3, 4 * 3)
    Qi[1:3, 1:3]     = noise.σan * dt^2 * Matrix{T}(I, 3, 3) # velocity noise Vi (m^2/s^2)
    Qi[4:6, 4:6]     = noise.σωn * dt^2 * Matrix{T}(I, 3, 3) # orientation noise Θi (rad^2)
    Qi[7:9, 7:9]     = noise.σaw * dt * Matrix{T}(I, 3, 3) # accelero bias noise Ai (m^2/s^4)
    Qi[10:12, 10:12] = noise.σωw * dt * Matrix{T}(I, 3, 3) # gyro bias noise Ωi (rad^2/s^2)

    # Error-state covariance
    new_error_cov = Fx * error_cov * Fx' + Fi * Qi * Fi'

    return new_nominal, error, new_error_cov
end

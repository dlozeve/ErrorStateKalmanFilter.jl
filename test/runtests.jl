using ErrorStateKalmanFilter
using Test
using LinearAlgebra
using Quaternions

@testset "types" begin
    x = collect(range(0, length=19, step=0.5))
    @test vec(NominalState(x)) == x
    x = collect(range(0, length=18, step=0.1))
    @test vec(ErrorState(x)) == x
end

@testset "ESKF" begin
    @test skew_matrix([1, 2, 3]) == [0 -3 2; 3 0 -1; -2 1 0]

    dt = 1e-3
    times = 0:dt:5

    nominal = NominalState(
        [0.0, 0, 0],
        [0.0, 0, 0],
        quat(1.0),
        zeros(3),
        zeros(3),
        [0.0, 0.0, -1.0]
    )

    error = ErrorState(zeros(18))

    error_cov = Matrix{Float64}(I, 18, 18)

    measurement = Measurement([0.0, 0, 0], [0.0, 0, 0])

    noise = Noise(
        2e-7  * (1.0/512),
        1e-10 * (1.0/512),
        1e-9,
        1e-9
    )

    new_nominal, new_error, new_error_cov = predict(nominal, error, error_cov, measurement, noise, dt)

    # error state is unchanged
    @test new_error == error
    # error covariance increases
    @test all(abs.(new_error_cov) .>= abs.(error_cov))
end

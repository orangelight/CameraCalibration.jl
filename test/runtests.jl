using CameraCalibration, Test
@testset "general" begin
    @test true
end
"""
using CSV, DataFrames
imgpoints = Array{Float64, 2}[]
for i = 0:10
    push!(imgpoints, vcat(convert(Array{Float64, 2}, CSV.File(string("./test/imgpoints", i, ".csv"), header = false) |> DataFrame)', ones(1, 35)))
end
"""

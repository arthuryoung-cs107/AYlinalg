include("../JULIA_AUXILIARY/AYio.jl")

mat0_check = AYio.read_matrix("mat0_io_test_matlab")

AYio.write_matrix(mat0_check, "mat0_io_test_julia")

println(mat0_check)
println(size(mat0_check))
println(typeof(mat0_check))

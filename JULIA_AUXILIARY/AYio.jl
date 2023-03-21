module AYio

export read_matrix, write_matrix

function read_matrix(name_::String, prefix_::String="../dat_dir/")
    full_name = prefix_*name_*".aydat"

    file = open(full_name,"r")
    header = Vector{Int32}(undef,5)
    read!(file,header)
    mat_raw = Matrix{Float64}(undef,header[5],header[4])
    read!(file,mat_raw)
    close(file)
    mat_out = Matrix(transpose(mat_raw));
    return mat_out
end

function write_matrix(mat_::Matrix, name_::String, prefix_::String="../dat_dir/")
    full_name = prefix_*name_*".aydat"
    dims = size(mat_)
    header = convert(Array{Int32,1}, [2, 2, dims[1]*dims[2]])
    ivec = convert(Array{Int32,1}, [dims[1], dims[2]])
    dvec = convert(Array{Float64,2}, transpose(mat_))
    
    file = open(full_name,"w")
    write(file,header,ivec,dvec)
    close(file)
    return nothing
end

end

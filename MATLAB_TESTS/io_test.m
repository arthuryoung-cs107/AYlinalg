clear
close all

test_path = cd('../MATLAB_AUXILIARY/');

mat0 = reshape(1:12,4,3)';
AYio.write_matrix(mat0,'mat0_io_test_matlab');

mat0_check_matlab = AYio.read_matrix('mat0_io_test_matlab');
mat0_check_julia = AYio.read_matrix('mat0_io_test_julia');

mat0
mat0_check_matlab
mat0_check_julia

cd(test_path);

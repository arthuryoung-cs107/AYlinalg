mat1 = aysml_read('../dat_dir/mat1')
fprintf_AYdata(mat1, '../dat_dir/mat1_back');

tens1 = aysml_read('../dat_dir/tens1')
fprintf_AYdata(tens1, '../dat_dir/tens1_back');

sym = [1.000000 2.000000 3.000000 4.000000; 2.000000 5.000000 6.000000 7.000000; 3.000000 6.000000 8.000000 9.000000; 4.000000 7.000000 9.000000 10.000000];
vec123 = (1:4)';
sym_check = sym*vec123

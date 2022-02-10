clear
close all

mat1 = AYdata.aysml_read('../dat_dir/mat1')
vec1 = AYdata.aysml_read('../dat_dir/vec1')
sym1 = AYdata.aysml_read('../dat_dir/sym1')
tens1 = AYdata.aysml_read('../dat_dir/tens1')
mat2 = AYdata.aysml_read('../dat_dir/mat2')
vec2 = AYdata.aysml_read('../dat_dir/vec2')
sym2 = AYdata.aysml_read('../dat_dir/sym2')
tens2 = AYdata.aysml_read('../dat_dir/tens2')
tens1_split = AYdata.aysml_read('../dat_dir/tens1_split')
tens2_split = AYdata.aysml_read('../dat_dir/tens2_split')

AYdata.aydat_write(mat1, '../dat_dir/mat1_b');
AYdata.aydat_write(vec1, '../dat_dir/vec1_b');
AYdata.aydat_write(tens1, '../dat_dir/tens1_b');
AYdata.aydat_write(mat2, '../dat_dir/mat2_b');
AYdata.aydat_write(vec2, '../dat_dir/vec2_b');
AYdata.aydat_write(tens2, '../dat_dir/tens2_b');
AYdata.aydat_write_split_tens(tens1_split, '../dat_dir/tens1_split_b');
AYdata.aydat_write_split_tens(tens2_split, '../dat_dir/tens2_split_b');

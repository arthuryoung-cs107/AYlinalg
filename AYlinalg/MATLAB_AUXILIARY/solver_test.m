clear
close all

sym05 = AYdata.aysml_read('../dat_dir/sym05')
sym1 = AYdata.aysml_read('../dat_dir/sym1')
b1 = AYdata.aysml_read('../dat_dir/b1')
xtrue = AYdata.aysml_read('../dat_dir/xtrue')

xtrue_mtlb = sym1\b1

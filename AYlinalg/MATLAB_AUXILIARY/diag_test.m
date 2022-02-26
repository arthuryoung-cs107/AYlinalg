clear
close all
mat1 = [0.988853 0.532895 0.885593 0.585933; ...
0.765989 0.845969 0.693015 0.324607; ...
0.125836 0.126754 0.122653 0.463925; ...
0.690825 0.141597 0.395719 0.681022; ...
0.878038 0.599037 0.620764 0.498526]

diag1 = [1.000000 0.000000 0.000000 0.000000 0.000000; ...
0.000000 2.000000 0.000000 0.000000 0.000000; ...
0.000000 0.000000 3.000000 0.000000 0.000000; ...
0.000000 0.000000 0.000000 4.000000 0.000000; ...
0.000000 0.000000 0.000000 0.000000 5.000000]

sym1 = (mat1')*(diag1)*(mat1)

diag1_vec = AYdata.aysml_read('../dat_dir/diag1')
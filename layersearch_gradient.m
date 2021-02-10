clear all; close all
load('HPPC23_H4_train.mat')
dI = gradient(I_data);
dV = gradient(V_cell); 
input_train = [I_data, V_cell', dI, dV'];
output_train = [soc_bulk_n'];
load('UDDS23_H1_train.mat')
dI = gradient(I_data);
dV = gradient(V_cell);
input_train = [input_train; I_data, V_cell', dI, dV'];
output_train = [output_train; soc_bulk_n'];

load('US0623_H1_train.mat')
dI = gradient(I_data);
dV = gradient(V_cell);
input_test = [I_data, V_cell', dI, dV'];
output_test = [soc_bulk_n'];

for i = 1:5 
    [train_error, test_error] = createNNlayer(input_train, output_train, input_test, output_test, i, 5)
    TR(i) = train_error; 
    TE(i) = test_error; 
end

figure
semilogy(total_tr)
hold on
semilogy(total_te)
%This code trains/validates a 1-layer NN to tune the number of neurons.
%Four inputs into the NN: I, V, dI, dV. 
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

for i = 1:30 
    [train_error, test_error] = createNN(input_train, output_train, input_test, output_test, i)
    TR(i) = train_error; 
    TE(i) = test_error; 
end

plot(TR)
hold on
plot(TE)
legend('Training Error','Testing Error') 



%This code looks trains and validates a neural net to tune the number of
%layers, given 5 neurons. 
clear all; close all
load('HPPC23_H4_train.mat')
input_train = [I_data, V_cell'];
output_train = [soc_bulk_n'];
load('UDDS23_H1_train.mat')
input_train = [input_train; I_data, V_cell'];
output_train = [output_train; soc_bulk_n'];

load('US0623_H1_train.mat')
input_test = [I_data, V_cell'];
output_test = [soc_bulk_n'];

for i = 1:5 
    [train_error, test_error] = createNNlayer(input_train, output_train, input_test, output_test, i, 5)
    TR(i) = train_error; 
    TE(i) = test_error; 
end

plot(TR)
hold on
plot(TE)




%Testing 3 Layer, 5 Neurons. 

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

load('test_set_2_spm.mat')
dI = gradient(I_data);
dV = gradient(V_cell);
input_test = [I_data, V_cell', dI, dV'];
output_test = [SOC_est'];

hiddenlayersize = [5, 5, 5];
net = fitnet(hiddenlayersize, 'trainlm');
net.divideFcn = 'divideind'; 
net.divideParam.trainInd = 1:1:length(output_train);
net.divideParam.valInd = [];
net.divideParam.testInd = [];
net.trainParam.epochs = 300; 

[trained_net, tr] = train(net, input_train', output_train');
%predicted_train = trained_net(input_train'); 

%train_err = sum((predicted_train' - output_train).^2)/length(predicted_train); 
view(trained_net)
predicted_test = trained_net(input_test');
%filtered = medfilt1(predicted_test, 400);
%filtered = movmean(predicted_test, [1300 0]); 
filtered = movavg(predicted_test', 'exponential', 500); 
filtered(1) = predicted_test(1);
test_err = sum((filtered - output_test).^2)/length(filtered); 
rmse = sqrt(sum((output_test - filtered).^2)/length(filtered))*100; 

t = (0:1:length(filtered)-1)*0.1;
plot(t, output_test, 'linewidth',2)
hold on
plot(t, filtered, 'linewidth',2)
xlabel('Time (s)', 'fontweight','bold','fontsize',12);
ylabel('SOC','fontweight','bold','fontsize',12);





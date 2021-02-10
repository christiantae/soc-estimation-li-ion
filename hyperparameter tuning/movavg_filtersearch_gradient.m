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

hiddenlayersize = [5, 5, 5];
net = fitnet(hiddenlayersize, 'trainlm');
net.divideFcn = 'divideind'; 
net.divideParam.trainInd = 1:1:length(output_train);
net.divideParam.valInd = [];
net.divideParam.testInd = [];
net.trainParam.epochs = 300; 

[trained_net, tr] = train(net, input_train', output_train');

predicted_test = trained_net(input_test');
n = [5, 10, 50, 100, 200, 250, 300, 340, 360, 380, 400, 420, 440, 460, 480, 500, 550, 600, 650, 700, 725, 750, 800, 1000]; 
for i = 1:length(n)
    filtered = movmean(predicted_test, [n(i) 0]); %kb, kf  
    filtered(1) = predicted_test(1); 
    test_err = sqrt(sum((filtered' - output_test).^2)/length(filtered))*100;
    TE(i) = test_err;      
end

plot(n,TE)

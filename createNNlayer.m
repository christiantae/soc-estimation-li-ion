function [train_err, test_err] = createNNlayer(input_train, output_train, input_test, output_test, layer, neurons)
%Wrapper function that creates a 2 Layer NN, varying the number of neurons
%   
hiddenlayersize = [];
for i = 1:layer
    hiddenlayersize = [hiddenlayersize, neurons];
end
net = fitnet(hiddenlayersize, 'trainlm');
net.divideFcn = 'divideind'; 
net.divideParam.trainInd = 1:1:length(output_train);
net.divideParam.valInd = [];
net.divideParam.testInd = [];
net.trainParam.epochs = 300; 

[trained_net, tr] = train(net, input_train', output_train');
predicted_train = trained_net(input_train'); 

train_err = sum((predicted_train' - output_train).^2)/length(predicted_train); 

predicted_test = trained_net(input_test');
test_err = sum((predicted_test' - output_test).^2)/length(predicted_test) 

end


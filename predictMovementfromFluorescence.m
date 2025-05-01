% learns patterns in speed, angular velocity, and then sees if fluorescence
% data can predict these patterns

% W Scott Conrad 14/2/25 <3 

clear all
tankfolder = '\\vs03.herseninstituut.knaw.nl\VS03-CSF-1\Conrad\Innate_approach\Data_analysis\24.35.01\ZI approach\Prey Laser\';
d = open([tankfolder 'allDatComb.mat']);
 

fluorescence_data = d.allDatComb.trialSignal; % Your fluorescence measurements
speed_data1 = d.allDatComb.trialSpeed; % Your speed measurements
speed_data = zeros(size(speed_data1, 1)*2, size(speed_data1, 2));

for i = 1:size(speed_data1, 1)
    speed_data(2*i-1, :) = speed_data1(i, :);  % Copy the i-th row of A to the 2*i-1-th row of B
    speed_data(2*i, :) = speed_data1(i, :);    % Copy the i-th row of A to the 2*i-th row of B
end


speedDataGC = normalize(reshape(speed_data', 1, []));
fluorDataGC = normalize(reshape(fluorescence_data', 1, []));
ts = 1:length(speedDataGC);
figure;
plot(ts, speedDataGC, '-b');
hold on;
plot(ts, fluorDataGC,'-k');

numseries = 2;
numlags = (1:20)';
nummdls = numel(numlags);
tbl = table(speedDataGC', fluorDataGC');
T = size(tbl,1);

% Partition time base.
maxp = max(numlags); % Maximum number of required presample responses
idxpre = 1:maxp;
idxest = (maxp + 1):T;

% Preallocation
EstMdl(nummdls) = varm(numseries,0);
aic = zeros(nummdls,1);

% Fit VAR models to data.
Y0 = tbl{idxpre,:}; % Presample
Y = tbl{idxest,:};  % Estimation sample
for j = 1:numel(numlags)
    Mdl = varm(numseries,numlags(j));
    Mdl.SeriesNames = tbl.Properties.VariableNames;
    EstMdl(j) = estimate(Mdl,Y,'Y0',Y0);
    results = summarize(EstMdl(j));
    aic(j) = results.AIC;
end

[~,bestidx] = min(aic);
p = numlags(bestidx)

for x = 1:20
    x
    BestMdl = EstMdl(x);
    h = gctest(BestMdl)
end


% Partition the data
numTrials = size(speed_data, 1);

X = normalize(fluorescence_data'); %might need to normalize training data first, then test using same nrmlztn factors as train
T = normalize(speed_data');

% Partition the data
cv = cvpartition(numTrials, 'HoldOut', 0.5);
X_train = X(:, cv.training);
T_train = T(:, cv.training);
X_test = X(:, cv.test);
T_test = T(:, cv.test);

% Create and configure the network
% hiddenLayerSize = 10;
% net = feedforwardnet(hiddenLayerSize);
numChannels = size(X_train,2);


layers = [
    sequenceInputLayer(numChannels)
    lstmLayer(128)
    fullyConnectedLayer(numChannels)];

% Configure network parameters
% net.divideFcn = ''; % Disable automatic data division since we're doing it manually
% net.performFcn = 'mse';
% net.trainFcn = 'trainlm'; % Levenberg-Marquardt backpropagation
% 
% % Configure training parameters
% net.trainParam.epochs = 100;
% net.trainParam.lr = 0.001;
% net.trainParam.min_grad = 1e-7;
% net.trainParam.showWindow = true;

options = trainingOptions("adam", ...
    MaxEpochs=200, ...
    SequencePaddingDirection="left", ...
    Shuffle="every-epoch", ...
    Plots="training-progress", ...
    Verbose=true, ...
    ExecutionEnvironment = "gpu");

% % Configure the network for the training data
% net = configure(net, X_train, T_train);

% Train the network
net = trainnet(X_train, T_train, layers, 'mse', options);

YTest = minibatchpredict(net,X_test, ...
    SequencePaddingDirection="left", ... % maybe right?
    UniformOutput=false);

% % Test the network
% Y_test = netTrained(X_test);

% % Make predictions
% predictions = predict(net, X_test);
% 

% Evaluate the performance by calculating Mean Squared Error (MSE)
mse = mean((YTest{:} - T_test).^2);  % Mean Squared Error between predicted and actual movement speed

% Display the result
disp(['Mean Squared Error on the test set: ', num2str(mse)]);

% Optionally, visualize predictions vs true values for the test set
figure;
plot(mean(T_test,2), 'b', 'LineWidth', 1.5);  % True movement data (1st trial as an example)
hold on;
plot(mean(YTest{:},2), 'r--', 'LineWidth', 1.5);  % Predicted movement data (1st trial as an example)
xlabel('Time Steps');
ylabel('Movement Speed');
legend('True', 'Predicted');
title('True vs Predicted Movement Speed (1st Test Trial)');


X_train = gpuArray(X_train);  % Move to GPU
T_train = gpuArray(T_train);  % Move to GPU
X_test = gpuArray(X_test);    % Move to GPU
T_test = gpuArray(T_test);    % Move to GPU
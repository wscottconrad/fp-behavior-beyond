function semData = sem(data)
% data = traceData;
% stdData = std(data);
% weeee = sqrt(length(data));
semData = std(data,0,2)/sqrt(length(data(1,:)));

% calculates standard error of mean, meant for time series data
% scott conrad 20/12/2024
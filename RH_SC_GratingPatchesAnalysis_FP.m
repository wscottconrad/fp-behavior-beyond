
%GratingPatchesAnalysis
% adapted from Robin Haak's script, Scott Conrad July 2025
% major changes: uses fiber photometry data, filters out trials where animal
% is moving 
%
% Robin Haak's original comments below:
%online analysis for 'GratingPatches' using NI stream times
%you need a GPU to run this script!
%Robin Haak, March 2023
% 
% vecChs = flip(1:385); %channels to plot (1= bottom Ch)
% 
% %% query user for file names & locations
% %ap.bin en ap.meta
% [strAp,strPathIm] = uigetfile('*.ap.bin','Select imec .ap.bin file','MultiSelect','off');
% sMetaIm = DP_ReadMeta(fullpath(strPathIm,[strAp(1:end-4) '.meta']));
% 
% %structEP
% [strLog,strLogPath] = uigetfile('*.mat','Select trial-based log file','MultiSelect','off');
% load(fullpath(strLogPath,strLog)); % %#ok<ows');
% 
% %% detect spikes on each channel
% [vecSpikeCh,vecSpikeT,~] = DP_DetectSpikesInBinaryFile(fullpath(strPathIm,strAp),[],[],'int16'); %strClass='int16'
% vecSpikeSecs = double(vecSpikeT)/str2double(sMetaIm.imSampRate)+...
%     str2double(sMetaIm.firstSample)/str2double(sMetaIm.imSampRate); %convert to seconds+add offset
% intNumChs = 385; %str2num(sMetaIm.nSavedChans); %#ok<ST2NM>

% preprocessing initialization
lowpass_cutoff = 3; %low-pass cut-off in Hz, e.g. 2
filt_steepness = .95; %how steep the filter transition is (0.5-1, default 0.85)
db_atten = 90; %decibel attenuation of filters (default = 60)
sr = 30; % sampling rate (Hz or FPS)
ma = 90*sr; % moving average window, default = 90 seconds
setUp = 30 * sr; % remove everything up to 30 sec before trials begin, get better fitted signal usually 
numChannels = 2; % number of channels you record from


% trace window
pre = 0.2; % 0.2 seconds before trial start
post = 0.5; % 0.5 seconds after
before = pre*sr; 
after = post*sr; 
traceTiming = -pre:(sr/(before+after)):(post-(sr/(before+after)));
trial_length = 3 * sr; % 3 seconds converted to frames
% corfactor = 0.133; % for ttl lag? 0.133 for fiber only

date = "2025_07_15" + ...
    "_";
animal = "111609";
filePath = strcat('\\vs03.herseninstituut.knaw.nl\VS03-CSF-1\Conrad\Innate_approach\Data_collection\24.35.01\', date, animal, '\');

% filePath = '\\vs03.herseninstituut.knaw.nl\VS03-CSF-1\Conrad\Innate_approach\Data_collection\24.35.01\2025_06_12_rfTest_fiberOnly\';
% signal = readtable(strcat(filePath, dates{d}, animalIDs{id}, '\', 'Fluorescence.csv'));
signal = readtable(strcat(filePath, 'Fluorescence.csv'));
events = readtable(strcat(filePath, 'Events.csv'));

stimDataPath = '\\vs03.herseninstituut.knaw.nl\VS03-CSF-1\Conrad\Innate_approach\Data_collection\24.35.01\RF mapping\';
% load(strcat(stimDataPath, '20250612_fiberOnly_02_16_48_16_RH_GratingPatches.mat'))
% Build the wildcard pattern
filePattern = strcat(stimDataPath, erase(date, '_'), '_', animal, '_withspeed*.mat');

% Get a list of matching files
files = dir(filePattern);

% Check if at least one match was found
if ~isempty(files)
    % Load the first matching file
    fullFileName = fullfile(files(1).folder, files(1).name);
    load(fullFileName);
else
    warning('No files matching pattern found: %s', filePattern);
end



timestamps = signal{:,1}/1000; % convert from ms to sec
eventtimestamps = events{:,1}(events{:,3}==0); %0 is ttl on ramp, 1 is off ramp

if strcmp(date, "2025_06_17_") || (strcmp(date, "2025_06_18_") && strcmp(animal, "112741"))
    ttlidx = 1;
elseif strcmp(date, "2025_06_19_") && strcmp(animal, "111609")
    ttlidx = 2;
elseif strcmp(date, "2025_07_15_") && (strcmp(animal, "111609") || strcmp(animal, "112742"))
    ttlidx = 4;
else
    ttlidx = 3;
end 

% clipStart = round((eventtimestamps(ttlidx,1)-setUp)/33.333); %where to clip data, in seconds
clipStart = round((eventtimestamps(ttlidx,1)-30000)/33.333); %where to clip data, in frames

clipEnd = round(((max((eventtimestamps), [], 'all')+ 10000))/33.333); % in frames
eventtimestamps = round(eventtimestamps(ttlidx:720+ttlidx-1)/33.3333)-clipStart;


for ch = 1:numChannels

chGreen = signal{clipStart:clipEnd,((2*ch)+2)};
chIsos = signal{clipStart:clipEnd,((2*ch)+1)};

% plot raw data
figure;
plot(timestamps(clipStart:clipEnd), chIsos, 'Color', '#7b8be8')
hold on;
plot(timestamps(clipStart:clipEnd), chGreen, 'Color', '#64a858')
hold on;
plot([(eventtimestamps'+clipStart)/sr; (eventtimestamps'+clipStart)/sr], repmat(ylim',1,length(eventtimestamps'/sr)), '-k')
box off
title([animal ' Channel ' num2str(ch) ' raw data'])
xlabel('Seconds since RWD recording start')
ylabel('AU')

% 4 lowpassfilter % changed conversion to sr
lp_normDatG = lpFilter(chGreen, sr, lowpass_cutoff,...
    filt_steepness, db_atten);
lp_normDatI = lpFilter(chIsos, sr, lowpass_cutoff,...
    filt_steepness, db_atten);

% 1 fit isosbestic signal to green channel (motion correction). if using controlfit instead of IRLS_dFF, try
% degree higher than 1
[dFoF, ft_iso_signal] = IRLS_dFF(lp_normDatG, lp_normDatI, 3); % suggested constant is 1.4



% plot filtered data % looks trash with fiber only data
figure;
plot(timestamps(clipStart:clipEnd)*sr-clipStart, dFoF)
hold on;
plot([eventtimestamps'; eventtimestamps'], repmat(ylim',1,length(eventtimestamps')), '-k')
title([animal, 'Ch' num2str(ch) ' filtered data'])
xlabel('Seconds (trimmed data)')
ylabel('AU')

%Make traces...
traces = NaN(length(eventtimestamps), before + trial_length);
tracesG = NaN(length(eventtimestamps), before + trial_length);
tracesI = NaN(length(eventtimestamps), before + trial_length);

AUCfitted = NaN(length(eventtimestamps),1);
AUC_G = NaN(length(eventtimestamps),1);
AUC_I = NaN(length(eventtimestamps),1);

% Compute the indices for the window
startIdx = eventtimestamps - before; 
AUCendIdx = eventtimestamps + after-1; 
endIdx = eventtimestamps + trial_length-1; 

% AUC of df/f trial traces
for m = 1:length(eventtimestamps)
    % Extract the window and store it
    % traces(m) = trapz(dFoF(eventtimestamps(m):endIdx(m))); % traces of fitted data
    traces(m) = ((dFoF(round(eventtimestamps(m):endIdx(m)))-median(dFoF(round(startIdx(m):eventtimestamps(m)-1)))))/median(dFoF(round(startIdx(m):eventtimestamps(m)-1))); % traces of fitted data
    tracesG(m) = ((lp_normDatG(round(eventtimestamps(m):endIdx(m)))-median(lp_normDatG(round(startIdx(m):eventtimestamps(m)-1)))))/median(lp_normDatG(round(startIdx(m):eventtimestamps(m)-1))); % traces of filtered green signal
    tracesI(m) = ((lp_normDatI(round(eventtimestamps(m):endIdx(m)))-median(lp_normDatI(round(startIdx(m):eventtimestamps(m)-1)))))/median(lp_normDatI(round(startIdx(m):eventtimestamps(m)-1))); % traces of filtered isosbestic
    
    AUCfitted(m) = trapz(((dFoF(round(eventtimestamps(m):AUCendIdx(m)))-median(dFoF(round(startIdx(m):eventtimestamps(m)-1)))))/median(dFoF(round(startIdx(m):eventtimestamps(m)-1)))); 
    AUC_G(m) = trapz(((lp_normDatG(round(eventtimestamps(m):AUCendIdx(m)))-median(lp_normDatG(round(startIdx(m):eventtimestamps(m)-1)))))/median(lp_normDatG(round(startIdx(m):eventtimestamps(m)-1)))); 
    AUC_I(m) = trapz(((lp_normDatI(round(eventtimestamps(m):AUCendIdx(m)))-median(lp_normDatI(round(startIdx(m):eventtimestamps(m)-1)))))/median(lp_normDatI(round(startIdx(m):eventtimestamps(m)-1)))); 

end

%% get movement data
speed = structEP.speed;
dt = datetime(structEP.speedTS, 'InputFormat', 'HH:mm:ss.SSS');
timeMillis = milliseconds(dt - dt(1)); % speed ms since recording start

 % sampling rate is not consistant
 % 0 is stim computer recording start
 

 % align ttls and 

dtTTL = datetime(structEP.ttlTimes, 'InputFormat', 'HH:mm:ss.SSS');
ttlMillis = milliseconds(dtTTL - dt(1)); % ttl ms from recording start

[dif,speedclip] = min(abs(timeMillis-(ttlMillis(1)-300)));
fprintf('first ttl -300 ms is offset from speed data by %d ms\n', dif)

timeMillis = timeMillis - timeMillis(speedclip);
speedSet = [timeMillis(speedclip:end)', abs(speed(speedclip:end)')];


% need to not hard code this
ts = timestamps*1000;
figure;
plot(timeMillis, abs(speed))
hold on;
plot(ts(1:length(dFoF(891:end))), dFoF(891:end))
hold on;
% plot([ttlMillis; ttlMillis], repmat(ylim',1,length(ttlMillis)), '-k')
plot([(eventtimestamps'-891)*33.3333; (eventtimestamps'-891)*33.3333], repmat(ylim',1,length(eventtimestamps')), '-k')
title([animal, 'Ch', ch, ' absolute wheel speed', 'vs dFoF'])



%% filter ttls where animal moves

ttlPre = 300; % ms before ttl 
ttlPost = 800; % ms after ttl 
speedThreshold = 0.05; % average, for wheelturning (0.01 leads to 300-400 nans)


for i = 1:720  

    [~,Pre_idx]=min(abs(timeMillis-(ttlMillis(i)-ttlPre)));
    [~,Post_idx]=min(abs(timeMillis-(ttlMillis(i)+ttlPost)));

    if mean(speedSet(Pre_idx:Post_idx,2)>speedThreshold) %is average thebest? maybe if any point gets above a certain value?
        traces(i) = NaN;
        tracesG = NaN;
        tracesI = NaN;
        
        AUCfitted(i) = NaN;
        AUC_G(i) = NaN;
        AUC_I(i) = NaN;
    end 


end 



%% get grid data
vecUniqueRects = unique(structEP.vecDstRect','rows'); %unique dst rects
vecUniqueStims = 1:length(vecUniqueRects);
vecStimIdx = zeros(size(structEP.vecDstRect,2),1);
for intStim = 1:length(vecUniqueRects)
    vecStimIdx(ismember(structEP.vecDstRect',vecUniqueRects(intStim,:),'rows')) = vecUniqueStims(intStim);
end

vecX_pix = unique(vecUniqueRects(:,1))+(vecUniqueRects(1,3)-unique(vecUniqueRects(1,1)))/2;
vecY_pix = unique(vecUniqueRects(:,2))+(vecUniqueRects(1,4)-unique(vecUniqueRects(1,2)))/2;



%% loop through data

% Get logical index of valid (non-NaN) trials
validTrials = ~isnan(AUCfitted);  

% Filter vecStimIdx and traces for valid trials only
vecStimIdx_valid = vecStimIdx(validTrials);
traces_valid = AUC_G(validTrials);

matAvgRespAll = NaN(numel(vecY_pix),numel(vecX_pix),1);
vecAUC = traces_valid;

matAvgResp = NaN(numel(vecY_pix), numel(vecX_pix));
for intLoc = vecUniqueStims
    matAvgResp(intLoc) = mean(vecAUC(vecStimIdx_valid == intLoc));
end
matAvgRespAll(:,:,1) = matAvgResp;


%% plot data
%interpolate
vecX_pix_interp = linspace(vecX_pix(1),vecX_pix(end),16);
vecY_pix_interp = linspace(vecY_pix(1),vecY_pix(end),9);

%get colormap(s)
cellColorMaps = RH_ColorMaps; % or colormap parula
%%
%loop through channels
set(0,'DefaultFigureWindowStyle','docked')

    matAvgRespAll_interp = matAvgRespAll(:,:,1);
    figure; hold on; 
    title([animal, 'Ch: ', num2str(ch)]);
    imagesc(vecX_pix_interp,vecY_pix_interp,matAvgRespAll_interp);
    set(gca, 'YDir','reverse');
    colormap(cellColorMaps{2});
    cb=colorbar;
    cb.Label.String='AUC of Z-score or df/F';
    axis image
    fixfig;


end

% %%
% function dblRateSpontaneous = computeRateSpontaneous(vecSpikeTimes,vecStimOnSecs,vecStimOffSecs,sParams)
% %compute unit's spontaneous/baseline rate
% intCount = 0;
% dblPeriod = 0;
% for intTrial=1:length(vecStimOffSecs)-1
%     intCount = intCount + ...
%         length(find(vecSpikeTimes>(vecStimOffSecs(intTrial)+sParams.dblSecsFromPrevStimOff) & ...
%         vecSpikeTimes<vecStimOnSecs(intTrial+1)));
%     dblPeriod = dblPeriod + vecStimOnSecs(intTrial+1) - (vecStimOffSecs(intTrial)+sParams.dblSecsFromPrevStimOff);
% end
% dblRateSpontaneous = intCount / dblPeriod;
% if dblPeriod<5
%     fprintf('Less than 5s to compute spontaneous rate.\n')
% end
% if intCount<10
%     fprintf('Less than 10 spikes to compute spontaneous rate.\n')
% end
% end

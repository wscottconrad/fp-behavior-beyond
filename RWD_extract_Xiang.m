% extracts data from RWD system, processes it, and collates across trials
% scott conrad 02/12/2024, adapted from Isis Alonso-Lozares

clear all
animalIDs = {'24_03761_01_' , '24_03733_08_', '24_03733_07_', '24_03733_01_'};
dates = {'2025_01_13-', '2025_01_14-'};
conditions = {'DCZ', 'noinjection', 'saline'};
filePath = '\\vs03.herseninstituut.knaw.nl\VS03-CSF-1\Conrad\Collaboration_Bosman\Data_collection\';
savePath = '\\vs03.herseninstituut.knaw.nl\VS03-CSF-1\Conrad\Collaboration_Bosman\Data_analysis\';
numChannels = 1; % number of channels you record from

% preprocessing initialization
lowpass_cutoff = 3; %low-pass cut-off in Hz, e.g. 2
filt_steepness = .95; %how steep the filter transition is (0.5-1, default 0.85)
db_atten = 90; %decibel attenuation of filters (default = 60)
sr = 30; % sampling rate (Hz or FPS)
ma = 90*sr; % moving average window, default = 90 seconds
setUp = 30 * sr; % how far before first trace
% trace window
pre = 5; % 5 seconds before
post = 8; % seconds after
before = pre*sr; 
after = post*sr; 
traceTiming = -pre:(sr/(before+after)):(post-(sr/(before+after)));
labelsX = -pre:5:post; 

for d = 1:length(dates)
    
for con = 1:length(conditions) % maybe tries to find animal/condition combos that don't exist
    
for id = 1:length(animalIDs)
    try 
        signal = readtable(strcat(filePath, dates{d}, '24_data\', dates{d}, animalIDs{id}, conditions{con}, '\', 'Fluorescence.csv'));
    catch 
        fprintf(strcat(dates{d}, animalIDs{id}, conditions{con}, ' does not exist'))
        continue 
    end 
    
%     eventimestamps = readtable(strcat(filePath, dates{1}, animalIDs{id}, '\', 'Eventimestamps.csv'));
    timestamps = signal{:,1}/1000;
    eventtimestamps = signal{:,2};
    Go = find(contains(eventtimestamps, 'Input2*2*0'));
    allTrials = find(contains(eventtimestamps, 'Input1*2*0'));
    noGo = zeros(length(Go),1);
    ng=1;
    for t = 1:length(allTrials)
        if ~ismember(0,(allTrials(t) < Go - 10) + (allTrials(t) > Go + 10)) % filter out out of sync TTL inputs,seperate noGo from all trials 
            noGo(ng) = allTrials(t);
            ng=ng+1;
        end 
    end 

    events = {Go, noGo};
    clipStart = min(cell2mat(events), [], 'all')-setUp; %where to clip data
    clipEnd = max(cell2mat(events), [], 'all')+ 15*sr;
    events = {Go-clipStart, noGo-clipStart};

    for ch = 1:numChannels
    
    chIsos = signal{clipStart:clipEnd,((2*ch)+1)};
    chGreen = signal{clipStart:clipEnd,((2*ch)+2)};
     
    % initialize (unnecessary?)
    fittedCh = {};  
    
    
    % plot raw data
%     figure;
%     plot(timestamps(clipStart:clipEnd), chIsos, 'Color', '#7b8be8')
%     hold on;
%     plot(timestamps(clipStart:clipEnd), chGreen, 'Color', '#64a858')
%     % hold on;
% %     plot([(events{1,1}'/+clipStart)/sr; (events{1,1}'+clipStart)/sr], repmat(ylim',1,size(events{1,1}',1)), '-k')
% %     hold on;
% %     plot([(events{1,2}'/+clipStart)/sr; (events{1,2}'+clipStart)/sr], repmat(ylim',1,size(events{1,2}',1)), '-g')
%     box off
%     title([animalIDs{id} ' ' conditions(con) ' raw data'])
    
    % 4 lowpassfilter % changed conversion to sr
    lp_normDatG = lpFilter(chGreen, sr, lowpass_cutoff,...
    filt_steepness, db_atten);
    lp_normDatI = lpFilter(chIsos, sr, lowpass_cutoff,...
    filt_steepness, db_atten);

    % 1 fit isosbestic signal to green channel (motion correction). if using controlfit instead of IRLS_dFF, try
    % degree higher than 1
    [dFoF, ft_iso_signal] = IRLS_dFF(lp_normDatG, lp_normDatI, 4.685); % default constant is 4.685, can go lower
    %     fittedCh = controlfit(lp_normDatG, lp_normDatI); % old fit function

    
    % 2 calc normalized dF/F (only use for control fit, not IRLS_dFF) 
    %     dFoF = ((lp_normDatG-fittedCh)./fittedCh)*100;
    
    % 3 hp and moving average do i need this?
    %     [hp_normDat, mov_normDat] = hpFilterNew(timestamps, dFoF, ma); %i think timestamps is time?
    

  % plot filtered data 
figure;
plot(timestamps(clipStart:clipEnd)*sr-clipStart, dFoF)
hold on;
plot([events{1,1}'; events{1,1}'], repmat(ylim',1,size(events{1,1}',1)), '-k')
hold on;
plot([events{1,2}'; events{1,2}'], repmat(ylim',1,size(events{1,2}',1)), '-g')
title([animalIDs{id}, 'Ch' num2str(ch) conditions{con} ' filtered data'])

%Make traces...
% Compute the indices for the window, per trial type
for ev = 1:length(events)
traces = {};
tracesGraw = {};
tracesIraw = {};
startIdx = events{ev} - before; 
endIdx = events{ev} + after-1; 

for m = 1:length(events{ev})
    % Extract the window and store it
    traces{end + 1} = dFoF(startIdx(m):endIdx(m)); % traces of fitted data
    tracesGraw{end+1} = lp_normDatG(startIdx(m):endIdx(m)); % traces of filtered green signal
    tracesIraw{end+1} = lp_normDatI(startIdx(m):endIdx(m)); % traces of filtered isosbestic

end

%colates
traceData = cell2mat(traces)';
traceDataG = cell2mat(tracesGraw)';
traceDataI = cell2mat(tracesIraw)';

% baseline corrected (z-scored) (collated)
traceDataSD = std(traceData(:,1:pre*sr-2*sr),0,2);
ZdFoF = (traceData - mean(traceData(:,1:pre*sr-2*sr),2))./traceDataSD;

traceDataSDG = std(traceDataG(:,1:pre*sr-2*sr),0,2);
Gdata = (traceDataG - mean(traceDataG(:,1:pre*sr-2*sr),2))./traceDataSDG;

traceDataSDI = std(traceDataI(:,1:pre*sr-2*sr),0,2);
Idata = (traceDataI - mean(traceDataI(:,1:pre*sr-2*sr),2))./traceDataSDI;

% visualize traces by trial and channel (heatmaps)
figure;
imagesc(ZdFoF);  
colorbar;  
xline(150, '--w', 'Linewidth', 1.5)
xlabel('Time');
ylabel('Trial');
title([animalIDs{id}, ' Event ' num2str(ev) ' dF/F']);
ax = gca;
box off
xticks(1:150:900);
xticklabels(labelsX)

sesdat.date = dates{d}; 
sesdat.mouse = animalIDs{id};
sesdat.conversion = sr;
sesdat.lp_normDat = lp_normDatG;
sesdat.ZdFoF = ZdFoF;
sesdat.Gdata = Gdata;
sesdat.Idata = Idata;
sesdat.channel = ch;
sesdat.events = ev; % use for multiple events
sesdat.condition = conditions{con};

save([savePath dates{d} animalIDs{id}  conditions{con} ' Event ' num2str(ev) '.mat'], 'sesdat')
clear sesdat traces tracesGraw tracesIraw signal

end % end of event block

end % end of channel block 

clear Go noGo eventtimestamps

end % end of animal block

end % end of condition block

end % end of date block




% go to combineData.m next


% to do later, seperate trials into approach/avoid/NR, or even stamp
% within laser trials based on action, direction, other behavior

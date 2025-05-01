% extracts data from RWD system, processes it, and collates traces

% default traces made for TTL inputs (trials), but additional scripts can be
% added here to make traces around other events of interest (e.g. movement
% initiation, movement during iter trial intervals, etc)

% scott conrad 02/12/2024, adapted from Isis Alonso-Lozares

clear all
animalIDs = {'105647' , '105648'};
dates = {'2024_12_12_', '2024_12_16_'};
filePath = '\\vs03.herseninstituut.knaw.nl\VS03-CSF-1\Conrad\Innate_approach\Data_collection\24.35.01\ZI approach\Prey Laser\';
savePath = '\\vs03.herseninstituut.knaw.nl\VS03-CSF-1\Conrad\Innate_approach\Data_analysis\24.35.01\ZI approach\Prey Laser\';
numChannels = 2; % number of channels you record from
ntFilePth = '\\vs03.herseninstituut.knaw.nl\VS03-CSF-1\Conrad\Innate_approach\Data_collection\Neurotar\';
trialIndexAll = {[0 0 0 1 1 0 0]', [2 2 2 2 2 0 2]', [0 1 0 0 0 0]', [2 2 0 2 0 2 1]'};


% preprocessing initialization
lowpass_cutoff = 3; %low-pass cut-off in Hz, e.g. 2
filt_steepness = .95; %how steep the filter transition is (0.5-1, default 0.85)
db_atten = 90; %decibel attenuation of filters (default = 60)
sr = 30; % sampling rate (Hz or FPS)
ma = 90*sr; % moving average window, default = 90 seconds
setUp = 120 * sr; % remove first 2 minutes? maybe get better fit

% trace window
pre = 5; % 5 seconds before
post = 25; % 25 seconds after
before = pre*sr; 
after = post*sr; 
traceTiming = -pre:(sr/(before+after)):(post-(sr/(before+after)));
labelsX = -pre:5:post; 

for d = 1:length(dates)
    
for id = 1:length(animalIDs)
    ntFile = strcat(ntFilePth, dates{d}, animalIDs{id}, '\data.mat'); 
    ttlFile = strcat(ntFilePth, animalIDs{id}, '_', regexprep(dates{d}, '_', ''), '_01_FP_ttl_triggers'); 
    signal = readtable(strcat(filePath, dates{d}, animalIDs{id}, '\', 'Fluorescence.csv'));
%     eventimestamps = readtable(strcat(filePath, dates{1}, animalIDs{id}, '\', 'Eventimestamps.csv'));
    timestamps = signal{:,1}/1000;
    eventtimestamps = find(contains(signal{:,2},'Input1*2*0')); 
    clipStart = min((eventtimestamps), [], 'all')-setUp; %where to clip data
    clipEnd = max((eventtimestamps), [], 'all')+ 30*sr;
    eventtimestamps = eventtimestamps-clipStart;
    
    for ch = 1:numChannels
    
    chIsos = signal{clipStart:clipEnd,((2*ch)+1)};
    chGreen = signal{clipStart:clipEnd,((2*ch)+2)};
   
%     time = 1:length(signal{setUp:end,((2*ch)+3)});
%     timestamps = (timestamps/1000); %convert values to seconds? 
    
    % initialize (unnecessary?)
    fittedCh = {}; dFoF = {}; 
    
    
    % plot raw data
%     figure;
%     plot(timestamps(clipStart:clipEnd), chIsos, 'Color', '#7b8be8')
%     hold on;
%     plot(timestamps(clipStart:clipEnd), chGreen, 'Color', '#64a858')
%     hold on;
%     plot([(eventtimestamps'+clipStart)/sr; (eventtimestamps'+clipStart)/sr], repmat(ylim',1,length(eventtimestamps'/sr)), '-k')
%     box off
%     title([animalIDs{id} ' Channel ' num2str(ch) ' raw data'])
    
    % 4 lowpassfilter % changed conversion to sr
    lp_normDatG = lpFilter(chGreen, sr, lowpass_cutoff,...
    filt_steepness, db_atten);
    lp_normDatI = lpFilter(chIsos, sr, lowpass_cutoff,...
    filt_steepness, db_atten);

    % 1 fit isosbestic signal to green channel (motion correction). if using controlfit instead of IRLS_dFF, try
    % degree higher than 1
    [dFoF, ft_iso_signal] = IRLS_dFF(lp_normDatG, lp_normDatI, 3); % suggested constant is 1.4
    %     fittedCh = controlfit(lp_normDatG, lp_normDatI); % old fit function

    
    % 2 calc normalized dF/F (only use for control fit, not IRLS_dFF) 
    %     dFoF = ((lp_normDatG-fittedCh)./fittedCh)*100;
    
    % 3 hp and moving average do i need this?
    %     [hp_normDat, mov_normDat] = hpFilterNew(timestamps, dFoF, ma); %i think timestamps is time?
    

% plot filtered data 
% figure;
% plot(timestamps(clipStart:clipEnd)*sr-clipStart, dFoF)
% hold on;
% plot([eventtimestamps'; eventtimestamps'], repmat(ylim',1,length(eventtimestamps')), '-k')
% title([animalIDs{id}, 'Ch' num2str(ch) ' filtered data'])

%Make traces...
traces = {};
tracesGraw = {};
tracesIraw = {};

% Compute the indices for the window
startIdx = eventtimestamps - before; 
endIdx = eventtimestamps + after-1; 

% tracesITI = cell(3,1);

% if d == 1 && id == 1
%     x = 1;
%     
% elseif d ==1 && id == 2
%     x = 2;
%     
% elseif d == 2 && id == 1
%     x = 3;
%     
% elseif d == 2 && id ==2
%     x = 4;
%     
% else 
%     fprintf('something has gone horribly wrong!')
% end 

[ITIidx, speedITI, trialIdxInit, speedTrials, ttl, movTrialIdx, ...
    speedTrialsMov] = nt_ITI_movement(ntFile, ttlFile, eventtimestamps, sr); 

% trial traces
for m = 1:length(eventtimestamps)
    % Extract the window and store it
    traces{end + 1} = dFoF(startIdx(m):endIdx(m)); % traces of fitted data
    tracesGraw{end+1} = lp_normDatG(startIdx(m):endIdx(m)); % traces of filtered green signal
    tracesIraw{end+1} = lp_normDatI(startIdx(m):endIdx(m)); % traces of filtered isosbestic

end

% trial traces aligned to movement
tracesInit = {};
tracesInitGraw = {};
tracesInitIraw = {};
if ~isempty(trialIdxInit)
    trialIdxInit = rmmissing(trialIdxInit);
    for m = 1:length(trialIdxInit)
        % Extract the window and store it
        tracesInit{end + 1} = dFoF(trialIdxInit(m)-before:trialIdxInit(m)+after-1); % traces of fitted data
        tracesInitGraw{end+1} = lp_normDatG(trialIdxInit(m)-before:trialIdxInit(m)+after-1); % traces of filtered green signal
        tracesInitIraw{end+1} = lp_normDatI(trialIdxInit(m)-before:trialIdxInit(m)+after-1); % traces of filtered isosbestic
        
    end
end 
% ITI traces  
tracesITI = [];
tracesITIGraw = {};
tracesITIIraw = {};
for m = 1:length(ITIidx)
        % Extract the ITI movement window and store it
        if ~isempty(ITIidx(m))
            tracesITI(end + 1, :) = dFoF(ITIidx(m)-before:ITIidx(m)+after-1); % traces of fitted data
            tracesITIGraw{end+1} = lp_normDatG(ITIidx(m)-before:ITIidx(m)+after-1); % traces of filtered green signal
            tracesITIIraw{end+1} = lp_normDatI(ITIidx(m)-before:ITIidx(m)+after-1); % traces of filtered isosbestic

        else
            tracesITI = NaN;
            tracesITIGraw = NaN;
            tracesITIGraw = NaN;
        end
end

%colates
traceData = cell2mat(traces)';
% tracesITI = cell2mat(tracesITI)'; 
traceDataG = cell2mat(tracesGraw)';
traceDataI = cell2mat(tracesIraw)';


% combine channels if desired
% traceDataCr = traceData';

% baseline corrected (z-scored) (collated)
traceDataSD = std(traceData(:,1:pre*sr),0,2);
ZdFoF = (traceData - mean(traceData(:,1:pre*sr),2))./traceDataSD;

traceDataSDG = std(traceDataG(:,1:pre*sr),0,2);
Gdata = (traceDataG - mean(traceDataG(:,1:pre*sr),2))./traceDataSDG;

traceDataSDI = std(traceDataI(:,1:pre*sr),0,2);
Idata = (traceDataI - mean(traceDataI(:,1:pre*sr),2))./traceDataSDI;

% and for trial initiation
if ~isempty(tracesInit)
    tracesInit = cell2mat(tracesInit)';
    ZdFoFinit = (tracesInit - mean(traceData(movTrialIdx,1:pre*sr),2))./traceDataSD(movTrialIdx);
    
    tracesInitGraw = cell2mat(tracesInitGraw)';
    InitGdata = (tracesInitGraw - mean(traceDataG(movTrialIdx,1:pre*sr),2))./traceDataSDG(movTrialIdx);

    tracesInitIraw = cell2mat(tracesInitIraw)';
    InitIdata = (tracesInitIraw - mean(traceDataI(movTrialIdx,1:pre*sr),2))./traceDataSDI(movTrialIdx);
else
    ZdFoFinit = NaN;
    InitGdata = NaN;
    InitIdata = NaN;
end

% and for ITI
if ~isempty(tracesITI)
    traceDataSDITI = std(tracesITI(:,1:pre*sr),0,2);
    ZdFoFITI = (tracesITI - mean(tracesITI(:,1:pre*sr),2))./traceDataSDITI;
    
    tracesITIGraw = cell2mat(tracesITIGraw)';
    traceDataSDG = std(tracesITIGraw(:,1:pre*sr),0,2);
    ITIGdata = (tracesITIGraw - mean(tracesITIGraw(:,1:pre*sr),2))./traceDataSDG;

    tracesITIIraw = cell2mat(tracesITIIraw)';
    traceDataSDI = std(tracesITIIraw(:,1:pre*sr),0,2);
    ITIIdata = (tracesITIIraw - mean(tracesITIIraw(:,1:pre*sr),2))./traceDataSDI;

else
    ZdFoFITI = NaN;
    ITIGdata = NaN;
    ITIIdata = NaN;
end
    
% visualize traces by trial and channel (heatmaps)
figure;
imagesc(ZdFoF);  
colorbar;  
xline(150, '--w', 'Linewidth', 1.5)
% yline(length(ZdFoF(:,1))/2+0.5, '-k', 'Linewidth', 1.5)
xlabel('Time');
ylabel('Trial');
title([animalIDs{id}, ' Channel ' num2str(ch) ' dF/F']);
ax = gca;
box off
xticks(1:150:900);
xticklabels(labelsX)

% figure;
% imagesc(ZdFoFITI);  
% colorbar;  
% xline(150, '--w', 'Linewidth', 1.5)
% % yline(length(ZdFoF(:,1))/2+0.5, '-k', 'Linewidth', 1.5)
% xlabel('Time');
% ylabel('ITI Trial');
% title([animalIDs{id}, ' Channel ' num2str(ch) ' dF/F']);
% ax = gca;
% box off
% xticks(1:150:900);
% xticklabels(labelsX)

% collated
% baselineSD = std(traceData(:,1:pre*sr)); % sd of baseline, for z-scoring later % check this!!!!!!!!!!!!!!!!
% traceDataCol = traceDataCol - mean(traceDataCol(1:pre*sr)); % baseline subtraction

% plot first trial of session
% figure;
% plot(traceTiming, ZdFoF(1,:), 'Color', '#7b8be8')
% hold on;
% xline(0, '--k'); yline(0, '--k');
% box off
% xlabel('Time (s)');
% ylabel('Normalized dF/F');
% title([animalIDs{id}, 'Channel ' num2str(ch) ' averaged dF/F over session']);

% zscoring
% ZdFoF = traceDataCol/baselineSD;

sesdat.session = d; 
sesdat.mouse = animalIDs{id};
sesdat.conversion = sr;
sesdat.lp_normDat = lp_normDatG;

sesdat.ZdFoF = ZdFoF;
sesdat.ZdFoFinit = ZdFoFinit;
sesdat.ZdFoFITI = ZdFoFITI; % ZdFoFITI is cell

sesdat.Gdata = Gdata;
sesdat.Idata = Idata;

sesdat.InitGdata = InitGdata;
sesdat.InitIdata = InitIdata;

sesdat.ITIGdata = ITIGdata;
sesdat.ITIIdata = ITIIdata;

sesdat.channel = ch;
if ch ==1
    sesdat.speedTrials = speedTrials;
    sesdat.speedITI = speedITI;
    sesdat.speedTrialsMov = speedTrialsMov;
    
else 
    sesdat.speedTrials = {};
    sesdat.speedITI = {};
    sesdat.speedTrialsMov = {};
    
end 
% sesdat.events = ev; % use for multiple events

save([savePath dates{d} animalIDs{id}  ' Channel ' num2str(ch) '.mat'], 'sesdat')
clear sesdat

    end % end of channel block 

end % end of animal block

end % end of date block




% go to combineData.m next


% to do later, stamp within laser trials based on action, direction, other behavior

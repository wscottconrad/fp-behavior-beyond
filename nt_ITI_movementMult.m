function [ITIidx, speedITI, speedTrials, ttl] = nt_ITI_movement(ntFile, ttlFile, eventtimestamps, sr)


% gets indices for periods of neurotar movement in intertrial intervals (ITI), in fp rate and time 
% W Scott Conrad 28/01/25


sr_nt = 100; % sample rate of neurotar
corFac = 20 * sr; % correction factor (estimated) in seconds
% trialIndexAll = trialIndexAll{1,x}; % delete after debugging

data= open(ntFile);
% data= open('\\vs03.herseninstituut.knaw.nl\VS03-CSF-1\Ren\Innate_approach\Data_collection\Neurotar\Track_[2024-12-12_13-53-00]_105647_session1\Track_[2024-12-12_13-53-00].mat');
% ts = data.neurotar_data.Time; % use for intrp1
ttl = readtable(ttlFile);
ttl.Time = seconds(ttl.DateTime - ttl.DateTime(1))+(ttl.Milliseconds)/1000-ttl.Milliseconds(1)/1000;

startIdx = 5*sr;
endIdx = 25*sr;
setUp = 120 * sr; % used to clip data later

ttl = (ttl.Time(2:3:end))*sr; % ttl is now matched to fp sr

clipStart = min((ttl), [], 'all')-setUp- corFac; %where to clip data
% clipEnd = max((ttl), [], 'all')+ 30*fp_sr;

% ttl = ttl - clipStart - corFac;
ttl = eventtimestamps;

speed = resample(data.neurotar_data.Speed,sr, sr_nt); %off by ~5 seconds compared to video
speed = speed(clipStart:end);
speedC = speed; % a copy for later

figure;
plot(1:length(speedC), speedC)
hold on;
plot([ttl'; ttl'], repmat(ylim', 1, size(ttl,1)), '-k')
title('speed with ttl')

speedTrials = {};
speedITI = cell(3,1);
ITIidx = cell(3,1);

% speed during prey laser movement
% separate based on approach avoid? i can do this later
for t = 1:length(ttl)
    speedTrials{end+1} = speed((ttl(t)-startIdx):(ttl(t)+endIdx));
    speedC((ttl(t)-startIdx):((ttl(t))+endIdx*1.5)) = 0; %erase trace speed for later analysis, + 10 seconds after trial end
    
end

figure;
plot(1:length(speedC), speedC)
hold on;
plot([ttl'; ttl'], repmat(ylim', 1, size(ttl,1)), '-k')
title('speed ttl removed')


% speed during ITI
threshH = 75;
speedIdxH = find(speedC >= threshH); % high threshold

threshM = 30; % changed from 75
speedIdxM = find(speedC >= threshM); % mid threshold

threshL = 15;
speedIdxL = find(speedC >= threshL); % low threshold

Idxs = {speedIdxH, speedIdxM, speedIdxL};
susThreshes = [45 30 10]; % sustained threshold 

for xx = 1:size(Idxs, 2)
    Idxs{1,xx} = Idxs{1,xx}(Idxs{1,xx}<=eventtimestamps(end));
    
    % speedIdxH = speedIdxH(speedIdxH<=eventtimestamps(end));
    
    for test = 1:length(Idxs{1,xx})
        if mean(speedC(Idxs{1,xx}(test):Idxs{1,xx}(test)+400)) <= susThreshes(xx) % filters out too small
            Idxs{1,xx}(test) = NaN;
        elseif xx > 1 && ((mean(speedC(Idxs{1,xx}(test):Idxs{1,xx}(test)+600)) > susThreshes(xx-1)) ...
                || (mean(speedC(Idxs{1,xx}(test)-150:Idxs{1,xx}(test))) > 10))% filter out too high + select initiation
            Idxs{1,xx}(test) = NaN;
        end
    end
    
    Idxs{1,xx} = rmmissing(Idxs{1,xx});
    
    %finding starting point for ITIs
    for iti = 1:length(Idxs{1,xx})
        if iti == 1 || (iti ~= 1 && Idxs{1,xx}(iti) > (Idxs{1,xx}(iti-1) + 30*sr)) % finds first index that meets threshold, ignores indexs within 30 seconds of each other
            
            speedITI{xx,1}(end+1,:) = speedC((Idxs{1,xx}(iti)-startIdx):(Idxs{1,xx}(iti)+endIdx))'; % traces of speed
            ITIidx{xx,1}(end+1,:) = Idxs{1,xx}(iti);
            
        end
    end


end 



% ITIidx = cell2mat(ITIidx)';
% speedITI = cell2mat(speedITI)';
speedTrials = cell2mat(speedTrials)';
figure;
if ~isempty(speedITI{1,1})
plot((1:size(speedITI{1,1},2)), mean(speedITI{1,1}, 1), 'Color', [.3 .3 .3])
hold on;
end
if ~isempty(speedITI{2,1})
plot((1:size(speedITI{2,1},2)), mean(speedITI{2,1}, 1), 'Color', [.6 .6 .6])
hold on;
end
if ~isempty(speedITI{3,1})
plot((1:size(speedITI{3,1},2)), mean(speedITI{3,1}, 1), 'Color', [.8 .8 .8])
hold on;
end
plot((1:size(speedTrials,2)), mean(speedTrials, 1), 'Color', [0.78 0 0])
title('mean ITI movement above threshold')


% i also want to find the 'start' of this movement initiation, so scan for
% periods of nothing
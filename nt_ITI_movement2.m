function [ITIidx, speedITI, trialIdxInit, speedTrials, ttl, movTrialIdx, ...
    speedTrialsMov] = nt_ITI_movement(ntFile, ttlFile, eventtimestamps, sr)


% gets indices for periods of neurotar movement in intertrial intervals (ITI), in fp rate and time

% this version tries to make the ITI trial criteria to match more closely
% to the trial initiation's

% gonna wait on this until i have more data

% W Scott Conrad 28/01/25


sr_nt = 100; % sample rate of neurotar
corFac = 20 * sr; % correction factor (estimated) in seconds
thresh = 30; % speed considered locomotion start, about half the mean speed during moving trials
colorz = {[0.9 0.48 0.42], [0.94 0.67 0.42], [0.804 0.81 0.31], [0.54 0.87 0.36], [0.37 0.92 0.72], ...
    [0.38 0.77 0.92], [0.53 0.38 0.92], [0.92 0.38 0.91]};

data= open(ntFile);
% data= open('\\vs03.herseninstituut.knaw.nl\VS03-CSF-1\Ren\Innate_approach\Data_collection\Neurotar\Track_[2024-12-12_13-53-00]_105647_session1\Track_[2024-12-12_13-53-00].mat');
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

speedTrials = zeros(length(ttl), startIdx+endIdx);
speedTrialsMov = zeros(length(ttl), startIdx+endIdx);

speedITI = [];
ITIidx = [];
trialIdxInit = NaN(length(ttl),1);

% speed during prey laser movement
% separate based on approach avoid? i can do this later
for t = 1:length(ttl)
    speedTrials(t, :) = speed((ttl(t)-startIdx):(ttl(t)+endIdx)-1)';
    
    % first find area of activity
    if find(speed(ttl(t):ttl(t)+450) >= thresh*2, 1)
        possiblePeak = find(speed(ttl(t):ttl(t)+450 )>= thresh*2);
       for pp = 1:length(possiblePeak)
           if mean(speed(ttl(t)+possiblePeak(pp):ttl(t)+possiblePeak(pp)+150)) >= thresh-5 % check if sustained
               
               % then look for first threshold point 1 second before
               init = find(speed(ttl(t)+possiblePeak(pp)-30:ttl(t)+possiblePeak(pp)) >= thresh, 1); % looks for initiation in the second before bout.
               trialIdxInit(t) = ttl(t) + possiblePeak(pp) - 30 + init;
               if trialIdxInit(t) < ttl(t)+15
                   fprintf('trial movement initiation detected in first 500 ms')
               end
               
               if trialIdxInit(t) < ttl(t)+3
                   fprintf('trial movement initiation also detected in first 100 ms!')
               end
              speedTrialsMov(t, :) = speed(trialIdxInit(t)-startIdx:trialIdxInit(t)+endIdx-1)';  
            break  
           end
           
       end 
    end 
    speedC((ttl(t)-startIdx):((ttl(t))+endIdx*1.5)) = 0; %erase trace speed for later analysis, + 10 seconds after trial end
    
end 

movTrialIdx = find(~isnan(trialIdxInit));
figure;
plot(1:length(speedC), speedC)
hold on;
plot([ttl'; ttl'], repmat(ylim', 1, size(ttl,1)), '-k')
title('speed ttl removed')


% speed during ITI

speedIdxPP = find(speedC >= thresh*2); % mid threshold
speedIdxPP = speedIdxPP(speedIdxPP<=eventtimestamps(end)); 

susThresh = 30; % sustained threshold 

% determine ITI periods with acceptable movement
for xx = 1:length(speedIdxPP)
    
    if mean(speed(ttl(t)+possiblePeak(pp):ttl(t)+possiblePeak(pp)+150)) >= thresh-5 % check if sustained
               
               % then look for first threshold point 1 second before
               init = find(speed(ttl(t)+possiblePeak(pp)-30:ttl(t)+possiblePeak(pp)) >= thresh, 1); % looks for initiation in the second before bout.
               trialIdxInit(t) = ttl(t) + possiblePeak(pp) - 30 + init;
              
              speedTrialsMov(t, :) = speed(trialIdxInit(t)-startIdx:trialIdxInit(t)+endIdx-1)';  
            break  
           end
           
       end 
end 
    
    speedC((ttl(t)-startIdx):((ttl(t))+endIdx*1.5)) = 0; %erase trace speed for later analysis, + 10 seconds after trial end
    
        if mean(speedC(speedIdxPP(xx):speedIdxPP(xx)+150)) >= thresh-5 ...
                || mean(speedC(speedIdx(xx)-150:speedIdx(xx)))  10 % filters out quick bouts and bouts not recently initiated
            
            init = find(speed(ttl(t)+possiblePeak(pp)-30:ttl(t)+possiblePeak(pp)) >= thresh, 1); % looks for initiation in the second before bout.
               trialIdxInit(t) = ttl(t) + possiblePeak(pp) - 30 + init;
            
            speedIdx(xx) = NaN;
      
        end
       
end 

speedIdx = rmmissing(speedIdx);

% ~25 mm/s seems like a reasonable initiation threshold 



% determine start points of ITI
    for iti = 1:length(speedIdx)
        if iti == 1 || (iti ~= 1 && (speedIdx(iti) > (speedIdx(iti-1) + 30*sr))) % finds first index that meets threshold, ignores indexs within 30 seconds of each other
            
            speedITI(end+1,:) = speedC((speedIdx(iti)-startIdx):(speedIdx(iti)+endIdx)-1)'; % ITI speed traces
            ITIidx(end+1,:) = speedIdx(iti); % ITI index
            
        end
    end



% ITIidx = cell2mat(ITIidx)';
% speedITI = cell2mat(speedITI)';

figure;
plot((1:size(speedITI,2)), mean(speedITI, 1), 'Color', [.6 .6 .6])
hold on;
plot((1:size(speedTrials,2)), mean(speedTrials, 1), 'Color', [0.78 0 0])
hold on;
for i = 1:length(ttl)
    plot((1:size(speedTrialsMov,2)), speedTrialsMov(i,:), 'Color', colorz{i})
    hold on;
end 
title('mean ITI movement above threshold')

end


% i also want to find the 'start' of this movement initiation, so scan for
% periods of nothing
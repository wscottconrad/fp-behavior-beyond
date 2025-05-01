function [ITIidx, ITItraces, speedTraces] = nt_ITI_movement(ntFile)

data= open(ntFile);
% data= open('\\vs03.herseninstituut.knaw.nl\VS03-CSF-1\Ren\Innate_approach\Data_collection\Neurotar\Track_[2024-12-12_13-53-00]_105647_session1\Track_[2024-12-12_13-53-00].mat');

ttl = find(double(data.neurotar_data.TTL_inputs));

sr = 100; % sample rate.... probably need to downsample
startIdx = 5*sr;
endIdx = 20*sr;
setUp = 120 * sr; % used to clip data later

clipStart = min((ttl), [], 'all')-setUp; %where to clip data
clipEnd = max((ttl), [], 'all')+ 30*sr;

speed = data.neurotar_data.Speed(clipStart:clipEnd);
speedC = data.neurotar_data.Speed(clipStart:clipEnd); % a copy for later

speedTraces = {};
ITItraces = {};
ITIidx = {};

% average speed during prey laser movement
% separate based on approach avoid?
for t = length(ttl)
    speedTraces{end+1} = speed((ttl(t)-startIdx):(ttl(t)+endIdx));
    speedC((ttl(t)-startIdx):(ttl(t)+endIdx*1.5)) = 0; %erase trace speed for later analysis, + 10 seconds after trial end
    
end

thresh = mean(speedTraces); % i'll need to adjust this
sd = std(speedTraces); % and this

% find intervals where speeds are similar to trials
speedIdx = find(speedC <= thesh + sd && speedC >= thresh - sd);

%finding starting point
for iti = length(speedIdx)
    if iti == 1 || (iti ~= 1 && speedIdx(iti) > speedIdx(iti-1) +15*sr)
        
        ITItraces{end+1} = speedC((speedIdx(iti)-startIdx):(speedIdx(iti)+endIdx));
        ITIidx{end+1} = speedIdX(iti);
        
    end
end




% i also want to find the 'start' of this movement initiation, so scan for
% periods of nothing


test = [1 1 0 1]';
find(test)
%% Scott Conrad 20/12/2024 adapted from Isis Alonso-Lozares 
% script to combine data from different animals
% takes data created from RWD_extract
 
%% define where the stuff is
clear all
tankfolder = '\\vs03.herseninstituut.knaw.nl\VS03-CSF-1\Conrad\Innate_approach\Data_analysis\24.35.01\ZI approach\Prey Laser\';

% type = 'All SA sessions data';
allDat = cell(4,1);
 
files = dir(fullfile(tankfolder)); % need to add .mat filter
files(ismember({files.name}, {'.', '..'})) = [];
for i = 1:length(files) %iterate through experiment folder
   data = open([tankfolder files(i).name]);
   data = data(1,1);
   allDat{1,i} = data.sesdat.ZdFoF;
   allDat{2,i} = data.sesdat.Gdata;
   allDat{3,i} = data.sesdat.Idata;
   allDat{4,i} = data.sesdat.speedTrials; % switched from 5
   allDat{5,i} = data.sesdat.ZdFoFITI; % this is a cell
   allDat{6,i} = data.sesdat.speedITI; % this is a cell
   
end
 % examples for when i need to split data by channel, event type, etc.
%    leftZIapproach = [];
%    rightZIapproach = [];
%    leftZIbehavior2 = [];
%    rightZIbehavior2 = [];   ZIapproach = [];
  allDatComb.trialSignal = vertcat(allDat{1,:});
  allDatComb.Gtrace = vertcat(allDat{2,:});
  allDatComb.Itrace = vertcat(allDat{3,:});
  allDatComb.trialSpeed = cell2mat(vertcat(allDat{4,:}));
  allDatComb.ITIsignal = allDat(5,:);
  allDatComb.ITIspeed = allDat(6,:);
  save([tankfolder 'allDatComb' '.mat'], 'allDatComb');
  
  
  
%   save([tankfolder 'allITIComb' '.mat'], 'allITIComb');
%   save([tankfolder 'allspeedTrials' '.mat'], 'allspeedTrials');
%   save([tankfolder 'allspeedITI' '.mat'], 'allspeedITI');

% traceTiming = 1:size(allGComb,2);
% figure;
% plot(traceTiming, mean(allGComb,1), 'Color', '#64a858');
% hold on;
% jbfill(traceTiming, mean(allGComb,1) + std(allGComb,0,1), mean(allGComb,1) - std(allGComb,0,1), [0.39 0.66 0.35], [0.39 0.66 0.35]);
% hold on;
% plot(traceTiming, mean(allIComb,1), 'Color', '#7b8be8');
% hold on;
% jbfill(traceTiming, mean(allIComb,1) + std(allIComb,0,1), mean(allIComb,1) - std(allIComb,0,1), [0.48 0.55 0.91]);
% box off
% title({'Filtered & Baseline Corrected', 'Green and Isosbestic Signals'});
% ylabel('Z-Score');
% xlabel('Time');

  % go to btstrpCombine next
 
% for l = 1:size(dat, 2)
% %    leftZIapproach = [leftZIapproach; dat{1, l}.approach];
% 
% 
%  
% end     
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
   
   % Z- scored data
   allDat{1,i} = data.sesdat.ZdFoF;
   allDat{5,i} = data.sesdat.ZdFoFITI;
   allDat{7,i} = data.sesdat.ZdFoFinit;
   
   % speed data
   allDat{4,i} = data.sesdat.speedTrials; % switched from 5
   allDat{6,i} = data.sesdat.speedITI; % this is a cell
   if ~isempty(data.sesdat.speedTrialsMov)
       allDat{8,i} = data.sesdat.speedTrialsMov(find(data.sesdat.speedTrialsMov(:,1)),:);
   end 
   
   % traces for green and isosbestic channels seperated 
   allDat{2,i} = data.sesdat.Gdata;
   allDat{3,i} = data.sesdat.Idata;
   allDat{9,i} = data.sesdat.InitGdata;
   allDat{10,i} = data.sesdat.InitIdata;
   allDat{11,i} = data.sesdat.ITIGdata;
   allDat{12,i} = data.sesdat.ITIIdata;
   
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
  allDatComb.ITIsignal = vertcat(allDat{5,:});
  allDatComb.ITIspeed = cell2mat(vertcat(allDat{6,:}));
  allDatComb.ZdFoFinit = vertcat(allDat{7,:});
  allDatComb.speedTrialsMov = vertcat(allDat{8,:});
  
  allDatComb.InitGdata = vertcat(allDat{9,:});
  allDatComb.InitIdata = vertcat(allDat{10,:});
  allDatComb.ITIGdata = vertcat(allDat{11,:});
  allDatComb.ITIIdata = vertcat(allDat{12,:});
  
  save([tankfolder 'allDatComb' '.mat'], 'allDatComb');
  
  
  
%   save([tankfolder 'allITIComb' '.mat'], 'allITIComb');
%   save([tankfolder 'allspeedTrials' '.mat'], 'allspeedTrials');
%   save([tankfolder 'allspeedITI' '.mat'], 'allspeedITI');



traceTiming = 1:size(allDatComb.InitGdata,2);
figure;

%green
plot(traceTiming, mean(allDatComb.InitGdata,1), 'Color', '#64a858');
hold on;
jbfill(traceTiming, mean(allDatComb.InitGdata,1) + std(allDatComb.InitGdata,0,1), ...
    mean(allDatComb.InitGdata,1) - std(allDatComb.InitGdata,0,1), [0.39 0.66 0.35], [0.39 0.66 0.35]);
hold on;

% isos
plot(traceTiming, mean(allDatComb.InitIdata,1), 'Color', '#7b8be8');
hold on;
jbfill(traceTiming, mean(allDatComb.InitIdata,1) + std(allDatComb.InitIdata,0,1), ...
    mean(allDatComb.InitIdata,1) - std(allDatComb.InitIdata,0,1), [0.48 0.55 0.91]);
box off
title({'Filtered & Baseline Corrected', 'Green and Isosbestic Signals', 'For Prey Approach'});
ylabel('Z-Score');
xlabel('Time');

  % go to btstrpCombine next
 
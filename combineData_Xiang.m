%% Scott Conrad 20/12/2024 adapted from Isis Alonso-Lozares 
% script to combine data from different animals
% takes data created from RWD_extract
 
%% define where the stuff is
clear all
tankfolder = '\\vs03.herseninstituut.knaw.nl\VS03-CSF-1\Conrad\Collaboration_Bosman\Data_analysis\';

% type = 'All SA sessions data';
allDat.base = cell(3,1);
allDat.DCZ = cell(3,1);
allDat.saline = cell(3,1);

 
files = dir(fullfile(tankfolder));
files(ismember({files.name}, {'.', '..'})) = [];
for i = 1:length(files) %iterate through experiment folder
   data = open([tankfolder files(i).name]);
   data = data(1,1);
   
   if contains(files(i).name, 'noinjection')
       allDat.base{1,1} = [allDat.base{1,1}; data.sesdat.ZdFoF];
       
   elseif contains(files(i).name, 'DCZ')
       allDat.DCZ{1,1} = [allDat.DCZ{1,1}; data.sesdat.ZdFoF];
       
   elseif contains(files(i).name, 'saline')
       allDat.saline{1,1} = [allDat.saline{1,1}; data.sesdat.ZdFoF];
       
   end

%    allDat{2,i} = data.sesdat.Gdata;
%    allDat{3,i} = data.sesdat.Idata;

   
end
 % examples for when i need to split data by channel, event type, etc.
%    leftZIapproach = [];
%    rightZIapproach = [];
%    leftZIbehavior2 = [];
%    rightZIbehavior2 = [];   ZIapproach = [];

%   allDatComb = vertcat(allDat.base{1,:});
%   allGComb = vertcat(allDat{2,:});
%   allIComb = vertcat(allDat{3,:});
  save([tankfolder 'allDatComb' '.mat'], 'allDat');

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
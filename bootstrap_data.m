function [bootsCI] = bootstrap_data(data,num_boots ,sig)
 
%%Bootstraps data adapted by Isis for conditioning photometry experimens
% 03-10-21 I combined the consec idx function into here so it is easier-
% Isis
 
% If fig_title not empty (i.e. []), plots kernel+bootstraps according to meta and sig - 
% "Trials" if subj_names empty (enter subj names if plotting kernel of subj_means)
 
%%inputs:
%data = photometry data 
%num_boots = number of bootstraps (I do 5000)
%sig = alpha, I put 0.01
 
%%Output:
% bootsCI = LCI+UCI vector
 
%  Copyright 2020 Philip Jean-Richard-dit-Bressel, UNSW Sydney
%  Based on Colin Clifford 2018 bootstrap_CI.m 
 
num_trials = size(data,1);
window = size(data,2);
 
% Minimum 2 trials (otherwise get funky signals due to inevitable crossing oscillations)
if num_trials > 3
 
  % Prep bootstrapping variables (one row for each bootstrap) ...
   data_boots = zeros(num_boots, window);
   bootsCI = zeros(2,window); 
 
 
   for b = 1:num_boots
      % bootstrap data + kernel across all trials ...
      trial_array = ceil((num_trials).*rand(1,num_trials));
      data_boots(b,:) = mean(data(trial_array,:));
   end
   
   %% Calculate bootstrap CI
   data_boots = sort(data_boots,1);
 
   lower_conf_index = ceil(num_boots*(sig/2))+1;
   upper_conf_index = floor(num_boots*(1-sig/2));
   
   bootsCI(1,:) = data_boots(lower_conf_index,:);
   bootsCI(2,:) = data_boots(upper_conf_index,:);
%    tmp = find(bootsCI(1,:)>0);
%    ind(1, :) = tmp(consec_idx(tmp, thresh));
%    clear tmp
%    tmp = find(bootsCI(2,:) <0);
%    ind(2, :) = tmp(consec_idx(tmp, thresh));
 
else
   fprintf('Less than 3 trials - bootstrapping skipped\n');
   bootsCI = NaN;
end
 
% %plot
% ts = linspace(-10, 20, size(data,2));
% plot(ts,mean(data), 'LineWidth', 2, 'Color', 'y')
% hold on
% jbfill(ts, bootsCI(1,:), bootsCI(2,:), 'y', 'y');
% line([ts(1) ts(end)], [0 0], 'LineStyle', '--', 'Color', 'k');
% ylinemin = 1.2*max(bootsCI(2,:));
% ylinemax = 1.3*max(bootsCI(2,:));
% hold on
% plot(ts(bootsCI(2,:)<0), ylinemin*ones(size(ts(bootsCI(2,:)<0),2)), 's', 'MarkerSize', 7, 'MarkerFaceColor','b','Color', 'b')
% plot(ts(bootsCI(1,:)>0), ylinemax*ones(size(ts(bootsCI(1,:)>0),2)), 's', 'MarkerSize', 7, 'MarkerFaceColor','k','Color', 'k')
% yl =get(gca, 'YLim');
% line([0 0], [yl(1) yl(2)], 'Color', 'k', 'LineWidth', 2)
% xlim([ts(1) ts(end)])
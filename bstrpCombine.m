%% Bootstrap analysis and combine plots for experiments
% scott conrad 20/12/2024, adapted from Isis Alonso-Lozares
% takes data created from combineData.m
clear all
tankfolder = '\\vs03.herseninstituut.knaw.nl\VS03-CSF-1\Conrad\Innate_approach\Data_analysis\24.35.01\ZI approach\Prey Laser\';
d = open([tankfolder 'allDatComb.mat']);
trialSignal = d.allDatComb.trialSignal;

trialSignalinit = d.allDatComb.ZdFoFinit;

ITIsignal = d.allDatComb.ITIsignal;

trialSpeed = d.allDatComb.trialSpeed;

speedInit = d.allDatComb.speedTrialsMov;

ITIspeed = d.allDatComb.ITIspeed;
%define variables to save
btsrp = [];
thres = 5; % consecutive threshold length (sr/lp cutoff *0.5) = (30/3 *0.5)

% creating index for trial types 0 = no movement; 1 = avoid prey laser; 
% 2 = approach

trialIndex = [0 0 0 1 1 0 0 0 0 0 1 1 0 0 2 2 2 2 2 0 2 2 2 2 2 2 0 2 0 1 0 0 0 0 0 1 0 0 0 0 2 2 0 2 0 2 1 2 2 0 2 0 2 1]'; % per channel

figure;
plot(1:900, mean(trialSignal([4:5, 30], :)))
title("animal 1 left")

figure;
plot(1:900, mean(trialSignal([11:12, 36], :)))
title("animal 1 right") % this one looks a bit suspicious 

figure;
plot(1:900, mean(trialSignal([15:19, 21, 41:42, 44, 46:47], :)))
title("animal 2 left")

figure;
plot(1:900, mean(trialSignal([22:26, 28, 48:49, 51, 53:54], :)))
title("animal 2 right")

trialIndexSingle = [0 0 0 1 1 0 0 2 2 2 2 2 0 2 0 1 0 0 0 0 2 2 0 2 0 2 1]'; % per animal

trialIndex(trialIndex == 2) = 1; % for any movement
trialIndexSingle(trialIndexSingle == 2) = 1; % for any movement

nothing = trialSignal((trialIndex == 0), :);
% avoid = allDatComb((trialIndex == 1), :);
% approach = allDatComb((trialIndex == 2), :);
movTrial = trialSignal((trialIndex == 1), :);

% for movement
MovTrialSpeed = trialSpeed((trialIndexSingle ==1), :);

%% colours for plotting
gree = [0.47,0.67,0.19];
lgris = [0.8,0.8,0.8];
gris = [0.65,0.65,0.65];
dgris = [0.3,0.3,0.3];
yel = [0.93,0.69,0.13];
dg = [0.33,0.42,0.20];
dy = [0.62,0.50,0.22];
red = [0.78 0 0];

greyIdx = {dgris, gris, lgris};
% trace windows
pre = 5;
post = 25;
ts = linspace(-pre, post, size(ITIspeed ,2));

% movement
[perm.speed, ~] = permTest_array(speedInit, ITIspeed, 1000);

tmp = find(perm.speed(1,:)<0.05);
id = tmp(consec_idx(tmp, thres));
plot(ts(id), (130)*ones(size(ts(id),2), 2), 's', 'MarkerSize', 5, 'MarkerFaceColor',red,'Color', red)
clear tmp id
hold on


% figure;
% plot(ts, mean(MovTrialSpeed), 'Color', gree)
% hold on;
% jbfill(ts, (mean(MovTrialSpeed)-(std(MovTrialSpeed, 0, 1)./sqrt(size(MovTrialSpeed, 1)))),...
%     (std(MovTrialSpeed, 0, 1)./sqrt(size(MovTrialSpeed, 1))+mean(MovTrialSpeed)), gree, gree)
% hold on;

plot(ts, mean(ITIspeed), 'Color', gris)
hold on
jbfill(ts, (mean(ITIspeed)-(std(ITIspeed, 0, 1)./sqrt(size(ITIspeed, 1)))), ...
    (std(ITIspeed, 0, 1)./sqrt(size(ITIspeed, 1))+mean(ITIspeed)), gris, gris)
hold on

plot(ts, mean(speedInit), 'Color', red)
jbfill(ts, (mean(speedInit)-(std(speedInit, 0, 1)./sqrt(size(speedInit, 1)))),...
    (std(speedInit, 0, 1)./sqrt(size(speedInit, 1))+mean(speedInit)), red, red)


legend off
box off
set(gca, 'color', 'none')
xline(0, '--k', 'TTL or Threshold')
 
% btsrp.nothing = bootstrap_data(nothing, 5000, 0.0001);
% btsrp.avoid = bootstrap_data(avoid, 5000, 0.0001);
% btsrp.approach = bootstrap_data(approach, 5000, 0.0001);
% btsrp.all = bootstrap_data(allDatComb, 5000, 0.0001);

btsrp.nothing = bootstrap_data(nothing, 5000, 0.0001);
btsrp.movTrial = bootstrap_data(movTrial, 5000, 0.0001);
btsrp.movITI = bootstrap_data(ITIsignal, 5000, 0.0001); 
btsrp.movInit = bootstrap_data(trialSignalinit, 5000, 0.0001);
 
% [tmp, ~] = permTest_array(right.SArew, right.SAac, 1000);

[perm.trialITI, ~] = permTest_array(ITIsignal, movTrial, 1000);
[perm.alignITI, ~] = permTest_array(ITIsignal, trialSignalinit, 1000);
% [perm.TrITI, ~] = permTest_array(movTrial, ITIsignal, 1000);

%  
% [permPUN.RewAc, ~] = permTest_array(right.PUNrew, right.PUNac, 1000);
% [permPUN.Rewinac, ~] = permTest_array(right.PUNrew, right.PUNinac, 1000);
% [permPUN.AcInac, ~] = permTest_array(right.PUNac, right.PUNinac, 1000);
% [permPUN.PunRew, ~] = permTest_array(right.PUNpun, right.PUNrew, 1000);
% [permPUN.PunAc, ~] = permTest_array(right.PUNpun, right.PUNac, 1000);
% [permPUN.PunInac, ~] = permTest_array(right.PUNpun, right.PUNinac, 1000);
%  
% [permTESTS.RewAc, ~] = permTest_array(right.TESTSrew, right.TESTSac, 1000);
% [permTESTS.Rewinac, ~] = permTest_array(right.TESTSrew, right.TESTSinac, 1000);
% [permTESTS.AcInac, ~] = permTest_array(right.TESTSac, right.TESTSinac, 1000);
 
%% approach combined accross hemispheres
ts = linspace(-pre, post, size(trialSignal ,2));
a = figure;

% plot(ts, mean(allDatComb), 'LineWidth', 2, 'Color', gree)
% hold on

% plot(ts, mean(nothing), 'LineWidth', 2, 'Color', gris)
% hold on
plot(ts,mean(movTrial), 'LineWidth', 2, 'Color', gree)
hold on;
jbfill(ts, (mean(movTrial)-(std(movTrial, 0, 1)./sqrt(size(movTrial, 1)))), ...
    (std(movTrial, 0, 1)./sqrt(size(movTrial, 1))+mean(movTrial)), gree, gree)
hold on;

plot(ts, mean(ITIsignal), 'Color', gris)
hold on;
jbfill(ts, (mean(ITIsignal)-(std(ITIsignal, 0, 1)./sqrt(size(ITIsignal, 1)))), ...
    (std(ITIsignal, 0, 1)./sqrt(size(ITIsignal, 1))+mean(ITIsignal)), gris, gris)
hold on;

plot(ts, mean(trialSignalinit), 'Color', red)
hold on;
jbfill(ts, (mean(trialSignalinit)-(std(trialSignalinit, 0, 1)./sqrt(size(trialSignalinit, 1)))), ...
    (std(trialSignalinit, 0, 1)./sqrt(size(trialSignalinit, 1))+mean(trialSignalinit)), red, red)

hold on;

line([ts(1) ts(end)], [0 0], 'LineStyle', '--', 'Color', 'k');
ylinemin = max(mean(trialSignal, 1));
ylinemax = 1.6*ylinemin;

%significance bars for bootstrap

% all
% tmp = find(btsrp.all(2,:)<0);
% id = tmp(consec_idx(tmp, thres));
% plot(ts(id), ylinemin*ones(size(ts(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',dg,'Color', dg)
% clear tmp id
% tmp = find(btsrp.all(1,:)>0);
% id = tmp(consec_idx(tmp, thres));
% plot(ts(id), ylinemax*ones(size(ts(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',dg,'Color', dg)
% clear tmp id

% tmp = find(btsrp.movInit(2,:)<0);
% id = tmp(consec_idx(tmp, thres));
% plot(ts(id), ylinemax*ones(size(ts(id),2), 2)+2, 's', 'MarkerSize', 7, 'MarkerFaceColor',red,'Color', red)
% clear tmp id
% hold on;
% 
% tmp = find(btsrp.movInit(1,:)>0);
% id = tmp(consec_idx(tmp, thres));
% plot(ts(id), ylinemax*ones(size(ts(id),2), 2)+2, 's', 'MarkerSize', 7, 'MarkerFaceColor',red,'Color', red)
% clear tmp id
% hold on;
% 
% 
% tmp = find(btsrp.movITI(2,:)<0);
% id = tmp(consec_idx(tmp, thres));
% plot(ts(id), ylinemax*ones(size(ts(id),2), 2)+2.3, 's', 'MarkerSize', 7, 'MarkerFaceColor',gris,'Color', gris)
% clear tmp id
% hold on;
% 
% tmp = find(btsrp.movITI(1,:)>0);
% id = tmp(consec_idx(tmp, thres));
% plot(ts(id), ylinemax*ones(size(ts(id),2), 2)+2.3, 's', 'MarkerSize', 7, 'MarkerFaceColor',gris,'Color', gris)
% clear tmp id
% hold on;
% 
% % 
% tmp = find(btsrp.movTrial(2,:)<0);
% id = tmp(consec_idx(tmp, thres));
% plot(ts(id), ylinemax*ones(size(ts(id),2), 2)+2.6, 's', 'MarkerSize', 7, 'MarkerFaceColor',dg,'Color', dg)
% clear tmp id
% hold on;
% 
% tmp = find(btsrp.movTrial(1,:)>0);
% id = tmp(consec_idx(tmp, thres));
% plot(ts(id), ylinemax*ones(size(ts(id),2), 2)+2.6, 's', 'MarkerSize', 7, 'MarkerFaceColor',dg,'Color', dg)
% clear tmp id
% hold on;


% plot permutation tests
tmp = find(perm.trialITI(1,:)<0.025);
id = tmp(consec_idx(tmp, thres));
plot(ts(id), (ylinemax-0.3)*ones(size(ts(id),2), 2)+1.8, 's', 'MarkerSize', 5, 'MarkerFaceColor',gree,'Color', gree)
clear tmp id 

tmp = find(perm.alignITI(1,:)<0.025);
id = tmp(consec_idx(tmp, thres));
plot(ts(id), (ylinemax-0.3)*ones(size(ts(id),2), 2)+1.7, 's', 'MarkerSize', 5, 'MarkerFaceColor',red,'Color', red)
clear tmp id 
% 
% tmp = find(perm.TrITI(1,:)<0.0167);
% id = tmp(consec_idx(tmp, thres));
% plot(ts(id), (ylinemax-0.3)*ones(size(ts(id),2), 2)+1.6, 's', 'MarkerSize', 5, 'MarkerFaceColor',gris,'Color', gris)
% clear tmp id 

hold on
% jbfill(ts, (mean(nothing)-(std(nothing, 0, 1)./sqrt(size(nothing, 1)))),(std(nothing, 0, 1)./sqrt(size(nothing, 1))+mean(nothing)), gris, gris)
% hold on
% jbfill(ts, (mean(ITIsignal)-(std(ITIsignal, 0, 1)./sqrt(size(ITIsignal, 1)))),(std(ITIsignal, 0, 1)./sqrt(size(ITIsignal, 1))+mean(ITIsignal)), yel, yel)
% hold on
% jbfill(ts, (mean(movTrial)-(std(movTrial, 0, 1)./sqrt(size(movTrial, 1)))),(std(movTrial, 0, 1)./sqrt(size(movTrial, 1))+mean(movTrial)), gree, gree)
legend off
box off
set(gca, 'color', 'none')
xline(0, 'k')
% xline(0, 'k', 'Prey Laser')

title('ZI Approach, bootstrap CI with alpha 0.0001')
%     saveas(a, 'C:\Users\Isis\surfdrive\1_Experiments\Nic photometry\Nic14LeftHemiSA.fig')
%     saveas(a, 'C:\Users\Isis\surfdrive\1_Experiments\Nic photometry\Nic14LeftHemiSA.png')
 
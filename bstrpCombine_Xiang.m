%% Bootstrap analysis and combine plots for experiments
% scott conrad 20/12/2024, adapted from Isis Alonso-Lozares
% takes data created from combineData.m
clear all
tankfolder = '\\vs03.herseninstituut.knaw.nl\VS03-CSF-1\Conrad\Collaboration_Bosman\Data_analysis\';
data = open([tankfolder 'allDatComb.mat']);
allDatBase = data.allDat.base{1,1};
allDatDCZ = data.allDat.DCZ{1,1};
allDatSaline = data.allDat.saline{1,1};
%define variables to save
btsrp = [];
thres = 5; % consecutive threshold length (sr/lp cutoff *0.5) = (30/3 *0.5)
 
%% colours for plotting
gree = [0.47,0.67,0.19];
gris = [0.65,0.65,0.65];
yel = [0.93,0.69,0.13];
dg = [0.33,0.42,0.20];
dy = [0.62,0.50,0.22];
dgray = [0.2, 0.2, 0.2];

% trace windows
pre = 5;
post = 8;
 
% btsrp.SArew = bootstrap_data(right.SArew, 5000, 0.0001);
% btsrp.SAac = bootstrap_data(right.SAac, 5000, 0.0001);
% btsrp.SAinac = bootstrap_data(right.SAinac, 5000, 0.0001);
btsrp.base = bootstrap_data(allDatBase, 5000, 0.0001);
btsrp.DCZ = bootstrap_data(allDatDCZ, 5000, 0.0001);
btsrp.saline = bootstrap_data(allDatSaline, 5000, 0.0001);

 
[perm.baseDCZ, ~] = permTest_array(allDatDCZ, allDatBase, 1000);
% [tmp, ~] = permTest_array(btsrp.saline, btsrp.DCZ, 1000);
[perm.baseSaline, ~] = permTest_array(allDatSaline, allDatBase, 1000);
[perm.DCZSaline, ~] = permTest_array(allDatSaline, allDatDCZ, 1000);


% permSA.RewAc = 
% [permSA.Rewinac, ~] = permTest_array(right.SArew, right.SAinac, 1000);
% [permSA.AcInac, ~] = permTest_array(right.SAac, right.SAinac, 1000);
 
%% approach combined accross hemispheres
ts = linspace(-pre, post, size(allDatBase ,2));
a = figure;
plot(ts, mean(allDatBase), 'LineWidth', 2, 'Color', gris)
hold on
plot(ts,mean(allDatDCZ), 'LineWidth', 2, 'Color', gree)
plot(ts,mean(allDatSaline), 'LineWidth', 2, 'Color', yel)
% jbfill(ts, btsrp.CSpOne(1,:), btsrp.CSpOne(2,:), 'k', 'k');
line([ts(1) ts(end)], [0 0], 'LineStyle', '--', 'Color', 'k');
ylinemin = max(mean(allDatDCZ, 1));
ylinemax = 1.1*ylinemin;
box off
xlabel('Time (s)')
ylabel('Z-Score')


%significance bars for bootstrap
% base (no injection)
% tmp = find(btsrp.base(2,:)<0);
% id = tmp(consec_idx(tmp, thres));
% plot(ts(id), ylinemin*ones(size(ts(id),2), 2)-0.3, 's', 'MarkerSize', 7, 'MarkerFaceColor',dgray,'Color', dgray)
% clear tmp id
% tmp = find(btsrp.base(1,:)>0);
% id = tmp(consec_idx(tmp, thres));
% plot(ts(id), ylinemax*ones(size(ts(id),2), 2)-0.3, 's', 'MarkerSize', 7, 'MarkerFaceColor',dgray,'Color', dgray)
% clear tmp id

% DCZ
% tmp = find(btsrp.DCZ(2,:)<0);
% id = tmp(consec_idx(tmp, thres));
% plot(ts(id), ylinemin*ones(size(ts(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',dg,'Color', dg)
% clear tmp id
% tmp = find(btsrp.DCZ(1,:)>0);
% id = tmp(consec_idx(tmp, thres));
% plot(ts(id), ylinemax*ones(size(ts(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',dg,'Color', dg)
% clear tmp id

% saline
% tmp = find(btsrp.saline(2,:)<0);
% id = tmp(consec_idx(tmp, thres));
% plot(ts(id), ylinemin*ones(size(ts(id),2), 2)+0.3, 's', 'MarkerSize', 7, 'MarkerFaceColor',dy,'Color', dy)
% clear tmp id
% tmp = find(btsrp.saline(1,:)>0);
% id = tmp(consec_idx(tmp, thres));
% plot(ts(id), ylinemax*ones(size(ts(id),2), 2)+0.3, 's', 'MarkerSize', 7, 'MarkerFaceColor',dy,'Color', dy)
% clear tmp id

% plot(ts(btsrp.approach(2,:)<0), ylinemin*ones(size(ts(btsrp.approach(2,:)<0),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',gree,'Color', gree)
% plot(ts(btsrp.approach(1,:)>0), ylinemax*ones(size(ts(btsrp.approach(1,:)>0),2),2), 's', 'MarkerSize', 7, 'MarkerFaceColor',gree,'Color', gree)
% hold on

tmp = find(perm.baseDCZ(1,:)<0.0167);
id = tmp(consec_idx(tmp, thres));
plot(ts(id), (ylinemax+0.01)*ones(size(ts(id),2), 2), 's', 'MarkerSize', 5, 'MarkerFaceColor',gris,'Color',  gris)
clear tmp id 

tmp = find(perm.baseSaline(1,:)<0.0167);
id = tmp(consec_idx(tmp, thres));
plot(ts(id), (ylinemax-0.3)*ones(size(ts(id),2), 2), 's', 'MarkerSize', 5, 'MarkerFaceColor',gree,'Color', gree)
clear tmp id 

tmp = find(perm.DCZSaline(1,:)<0.0167);
id = tmp(consec_idx(tmp, thres));
plot(ts(id), (ylinemax-0.05)*ones(size(ts(id),2), 2), 's', 'MarkerSize', 5, 'MarkerFaceColor',yel,'Color', yel)

hold on
jbfill(ts, (mean(allDatBase)-(std(allDatBase, 0, 1)./sqrt(size(allDatBase, 1)))),(std(allDatBase, 0, 1)./sqrt(size(allDatBase, 1))+mean(allDatBase)), gris, gris)
hold on
jbfill(ts, (mean(allDatDCZ)-(std(allDatDCZ, 0, 1)./sqrt(size(allDatDCZ, 1)))),(std(allDatDCZ, 0, 1)./sqrt(size(allDatDCZ, 1))+mean(allDatDCZ)), gree, gree)
hold on
jbfill(ts, (mean(allDatSaline)-(std(allDatSaline, 0, 1)./sqrt(size(allDatSaline, 1)))),(std(allDatSaline, 0, 1)./sqrt(size(allDatSaline, 1))+mean(allDatSaline)), yel, yel)
legend off
box off
set(gca, 'color', 'none')
xline(0, 'k')
% xline(0, 'k', 'Prey Laser')

title('Go/No Go trials combined, bootstrap CI with alpha 0.0001')
%     saveas(a, 'C:\Users\Isis\surfdrive\1_Experiments\Nic photometry\Nic14LeftHemiSA.fig')
%     saveas(a, 'C:\Users\Isis\surfdrive\1_Experiments\Nic photometry\Nic14LeftHemiSA.png')
 
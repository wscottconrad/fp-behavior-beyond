%% Bootstrap analysis and combine plots for Nic experiments
 
%define variables to save
btsrp = [];
 
%% colours for plotting
gree = [0.47,0.67,0.19];
gris = [0.65,0.65,0.65];
yel = [0.93,0.69,0.13];
dg = [0.33,0.42,0.20];
dy = [0.62,0.50,0.22];
 
 
 
% btsrp.SArew = bootstrap_data(right.SArew, 5000, 0.0001);
% btsrp.SAac = bootstrap_data(right.SAac, 5000, 0.0001);
% btsrp.SAinac = bootstrap_data(right.SAinac, 5000, 0.0001);
 btsrp.PUNrew = bootstrap_data(right.PUNrew, 5000, 0.0001);
 btsrp.PUNac   = bootstrap_data(right.PUNac, 5000, 0.0001);
 btsrp.PUNinac = bootstrap_data(right.PUNinac, 5000, 0.0001);
 btsrp.PUNnp = bootstrap_data(right.PUNpun, 5000, 0.0001);
% btsrp.TESTSrew = bootstrap_data(right.TESTSrew, 5000, 0.0001);
% btsrp.TESTSac = bootstrap_data(right.TESTSac, 5000, 0.0001);
% btsrp.TESTSinac = bootstrap_data(right.TESTSinac, 5000, 0.0001);
% right.btsrp = btsrp;
 
[tmp, ~] = permTest_array(right.SArew, right.SAac, 1000);
permSA.RewAc = 
[permSA.Rewinac, ~] = permTest_array(right.SArew, right.SAinac, 1000);
[permSA.AcInac, ~] = permTest_array(right.SAac, right.SAinac, 1000);
 
[permPUN.RewAc, ~] = permTest_array(right.PUNrew, right.PUNac, 1000);
[permPUN.Rewinac, ~] = permTest_array(right.PUNrew, right.PUNinac, 1000);
[permPUN.AcInac, ~] = permTest_array(right.PUNac, right.PUNinac, 1000);
[permPUN.PunRew, ~] = permTest_array(right.PUNpun, right.PUNrew, 1000);
[permPUN.PunAc, ~] = permTest_array(right.PUNpun, right.PUNac, 1000);
[permPUN.PunInac, ~] = permTest_array(right.PUNpun, right.PUNinac, 1000);
 
[permTESTS.RewAc, ~] = permTest_array(right.TESTSrew, right.TESTSac, 1000);
[permTESTS.Rewinac, ~] = permTest_array(right.TESTSrew, right.TESTSinac, 1000);
[permTESTS.AcInac, ~] = permTest_array(right.TESTSac, right.TESTSinac, 1000);
 
%% SA ALL
 
 
%% SA LEFT
ts = linspace(-5, 10, size(left.SArew ,2));
a = figure;
plot(ts, mean(left.SArew), 'LineWidth', 2, 'Color', gree)
hold on
plot(ts,mean(left.SAac), 'LineWidth', 2, 'Color', gris)
plot(ts,mean(left.SAinac), 'LineWidth', 2, 'Color', yel)
% jbfill(ts, btsrp.CSpOne(1,:), btsrp.CSpOne(2,:), 'k', 'k');
line([ts(1) ts(end)], [0 0], 'LineStyle', '--', 'Color', 'k');
ylinemin = max(mean(left.SArew, 1));
ylinemax = 1.6*ylinemin;
plot(ts(left.btsrp.SArew(2,:)<0), ylinemin*ones(size(ts(left.btsrp.SArew(2,:)<0),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',gree,'Color', gree)
plot(ts(left.btsrp.SArew(1,:)>0), ylinemax*ones(size(ts(left.btsrp.SArew(1,:)>0),2),2), 's', 'MarkerSize', 7, 'MarkerFaceColor',gree,'Color', gree)
hold on
plot(ts(left.btsrp.SAac(2,:)<0), (ylinemin-0.1)*ones(size(ts(left.btsrp.SAac(2,:)<0),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',gris,'Color', gris)
plot(ts(left.btsrp.SAac(1,:)>0), (ylinemax-0.1)*ones(size(ts(left.btsrp.SAac(1,:)>0),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',gris,'Color', gris)
plot(ts(left.btsrp.SAinac(2,:)<0), (ylinemin-0.2)*ones(size(ts(left.btsrp.SAinac(2,:)<0),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',yel,'Color', yel)
plot(ts(left.btsrp.SAinac(1,:)>0), (ylinemax-0.2)*ones(size(ts(left.btsrp.SAinac(1,:)>0),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',yel,'Color', yel)
plot(ts(permSA.RewAc(1,:)<0.01), (ylinemax-0.3)*ones(size(ts(permSA.RewAc(1,:)<0.01),2), 2), 's', 'MarkerSize', 5, 'MarkerFaceColor',gree,'Color', gree)
plot(ts(permSA.RewAc(1,:)<0.01), (ylinemax-0.33)*ones(size(ts(permSA.RewAc(1,:)<0.01),2), 2), 's', 'MarkerSize', 5, 'MarkerFaceColor',gris,'Color', gris)
plot(ts(permSA.Rewinac(1,:)<0.01), (ylinemax-0.4)*ones(size(ts(permSA.Rewinac(1,:)<0.01),2), 2), 's', 'MarkerSize', 5, 'MarkerFaceColor',gree,'Color', gree)
plot(ts(permSA.Rewinac(1,:)<0.01), (ylinemax-0.43)*ones(size(ts(permSA.Rewinac(1,:)<0.01),2), 2), 's', 'MarkerSize', 5, 'MarkerFaceColor',yel,'Color', yel)
plot(ts(permSA.AcInac(1,:)<0.01), (ylinemax-0.5)*ones(size(ts(permSA.AcInac(1,:)<0.01),2), 2), 's', 'MarkerSize', 5, 'MarkerFaceColor',gris,'Color', gris)
plot(ts(permSA.AcInac(1,:)<0.01), (ylinemax-0.53)*ones(size(ts(permSA.AcInac(1,:)<0.01),2), 2), 's', 'MarkerSize', 5, 'MarkerFaceColor',yel,'Color', yel)
hold on
jbfill(ts, (mean(left.SArew)-(std(left.SArew, 0, 1)./sqrt(size(left.SArew, 1)))),(std(left.SArew, 0, 1)./sqrt(size(left.SArew, 1))+mean(left.SArew)), gree, gree)
hold on
jbfill(ts, (mean(left.SAac)-(std(left.SAac, 0, 1)./sqrt(size(left.SAac, 1)))),(std(left.SAac, 0, 1)./sqrt(size(left.SAac, 1))+mean(left.SAac)), gris, gris)
hold on
jbfill(ts, (mean(left.SAinac)-(std(left.SAinac, 0, 1)./sqrt(size(left.SAinac, 1)))),(std(left.SAinac, 0, 1)./sqrt(size(left.SAinac, 1))+mean(left.SAinac)), yel, yel)
legend('Rewarded NP', 'Active not rewarded NP', 'Inactive NP')
legend boxoff
box off
set(gca, 'color', 'none')
  vline(0, 'k', 'Nosepoke')
title('SA sessions left hemisphere, 1-sample bootstrap with alpha 0.0001, perm_test alpha 0.01')
    saveas(a, 'C:\Users\Isis\surfdrive\1_Experiments\Nic photometry\Nic14LeftHemiSA.fig')
    saveas(a, 'C:\Users\Isis\surfdrive\1_Experiments\Nic photometry\Nic14LeftHemiSA.png')
 
%% SA right
b = figure;
plot(ts, mean(right.SArew), 'LineWidth', 2, 'Color', gree)
hold on
plot(ts,mean(right.SAac), 'LineWidth', 2, 'Color', gris)
plot(ts,mean(right.SAinac), 'LineWidth', 2, 'Color', yel)
% jbfill(ts, btsrp.CSpOne(1,:), btsrp.CSpOne(2,:), 'k', 'k');
line([ts(1) ts(end)], [0 0], 'LineStyle', '--', 'Color', 'k');
ylinemin = max(mean(right.SArew, 1));
ylinemax = 1.6*ylinemin;
plot(ts(right.btsrp.SArew(2,:)<0), ylinemin*ones(size(ts(right.btsrp.SArew(2,:)<0),2),2), 's', 'MarkerSize', 7, 'MarkerFaceColor',gree,'Color', gree)
plot(ts(right.btsrp.SArew(1,:)>0), ylinemax*ones(size(ts(right.btsrp.SArew(1,:)>0),2),2), 's', 'MarkerSize', 7, 'MarkerFaceColor',gree,'Color', gree)
hold on
plot(ts(right.btsrp.SAac(2,:)<0), (ylinemin-0.05)*ones(size(ts(right.btsrp.SAac(2,:)<0),2),2), 's', 'MarkerSize', 7, 'MarkerFaceColor',gris,'Color', gris)
plot(ts(right.btsrp.SAac(1,:)>0), (ylinemax-0.05)*ones(size(ts(right.btsrp.SAac(1,:)>0),2),2), 's', 'MarkerSize', 7, 'MarkerFaceColor',gris,'Color', gris)
plot(ts(right.btsrp.SAinac(2,:)<0), (ylinemin-0.1)*ones(size(ts(right.btsrp.SAinac(2,:)<0),2),2), 's', 'MarkerSize', 7, 'MarkerFaceColor',yel,'Color', yel)
plot(ts(right.btsrp.SAinac(1,:)>0), (ylinemax-0.1)*ones(size(ts(right.btsrp.SAinac(1,:)>0),2),2), 's', 'MarkerSize', 7, 'MarkerFaceColor',yel,'Color', yel)
plot(ts(permSA.RewAc(1,:)<0.01), (ylinemax-0.15)*ones(size(ts(permSA.RewAc(1,:)<0.01),2), 2), 's', 'MarkerSize', 5, 'MarkerFaceColor',gree,'Color', gree)
plot(ts(permSA.RewAc(1,:)<0.01), (ylinemax-0.17)*ones(size(ts(permSA.RewAc(1,:)<0.01),2), 2), 's', 'MarkerSize', 5, 'MarkerFaceColor',gris,'Color', gris)
plot(ts(permSA.Rewinac(1,:)<0.01), (ylinemax-0.2)*ones(size(ts(permSA.Rewinac(1,:)<0.01),2), 2), 's', 'MarkerSize', 5, 'MarkerFaceColor',gree,'Color', gree)
plot(ts(permSA.Rewinac(1,:)<0.01), (ylinemax-0.22)*ones(size(ts(permSA.Rewinac(1,:)<0.01),2), 2), 's', 'MarkerSize', 5, 'MarkerFaceColor',yel,'Color', yel)
plot(ts(permSA.AcInac(1,:)<0.01), (ylinemax-0.25)*ones(size(ts(permSA.AcInac(1,:)<0.01),2), 2), 's', 'MarkerSize', 5, 'MarkerFaceColor',gris,'Color', gris)
plot(ts(permSA.AcInac(1,:)<0.01), (ylinemax-0.27)*ones(size(ts(permSA.AcInac(1,:)<0.01),2), 2), 's', 'MarkerSize', 5, 'MarkerFaceColor',yel,'Color', yel)
 
hold on
jbfill(ts, (mean(right.SArew)-(std(right.SArew, 0, 1)./sqrt(size(right.SArew, 1)))),(std(right.SArew, 0, 1)./sqrt(size(right.SArew, 1))+mean(right.SArew)), gree, gree)
hold on
jbfill(ts, (mean(right.SAac)-(std(right.SAac, 0, 1)./sqrt(size(right.SAac, 1)))),(std(right.SAac, 0, 1)./sqrt(size(right.SAac, 1))+mean(right.SAac)), gris, gris)
hold on
jbfill(ts, (mean(right.SAinac)-(std(right.SAinac, 0, 1)./sqrt(size(right.SAinac, 1)))),(std(right.SAinac, 0, 1)./sqrt(size(right.SAinac, 1))+mean(right.SAinac)), yel, yel)
legend('Rewarded NP', 'Active not rewarded NP', 'Inactive NP')
legend boxoff
box off
set(gca, 'color', 'none')
  vline(0, 'k', 'Nosepoke')
title('SA sessions right hemisphere, 1-sample bootstrap with alpha 0.0001, Permtest alpha 0.01')
    saveas(b, 'C:\Users\Isis\surfdrive\1_Experiments\Nic photometry\Nic14RightHemiSA.fig')
    saveas(b, 'C:\Users\Isis\surfdrive\1_Experiments\Nic photometry\Nic14RightHemiSA.png')
 
%% PUN LEFT
c = figure;
plot(ts, mean(left.PUNrew), 'LineWidth', 2, 'Color', gree)
hold on
plot(ts,mean(left.PUNac), 'LineWidth', 2, 'Color', gris)
plot(ts,mean(left.PUNinac), 'LineWidth', 2, 'Color', yel)
plot(ts,mean(left.PUNnp), 'LineWidth', 2, 'Color', 'r')
% jbfill(ts, btsrp.CSpOne(1,:), btsrp.CSpOne(2,:), 'k', 'k');
line([ts(1) ts(end)], [0 0], 'LineStyle', '--', 'Color', 'k');
ylinemin = max(mean(left.PUNnp, 1));
ylinemax = 1.6*ylinemin;
plot(ts(left.btsrp.PUNrew(2,:)<0), ylinemin*ones(size(ts(left.btsrp.PUNrew(2,:)<0),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',gree,'Color', gree)
plot(ts(left.btsrp.PUNrew(1,:)>0), ylinemax*ones(size(ts(left.btsrp.PUNrew(1,:)>0),2),2), 's', 'MarkerSize', 7, 'MarkerFaceColor',gree,'Color', gree)
hold on
plot(ts(left.btsrp.PUNac(2,:)<0), (ylinemin-0.1)*ones(size(ts(left.btsrp.PUNac(2,:)<0),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',gris,'Color', gris)
plot(ts(left.btsrp.PUNac(1,:)>0), (ylinemax-0.1)*ones(size(ts(left.btsrp.PUNac(1,:)>0),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',gris,'Color', gris)
plot(ts(left.btsrp.PUNinac(2,:)<0), (ylinemin-0.2)*ones(size(ts(left.btsrp.PUNinac(2,:)<0),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',yel,'Color', yel)
plot(ts(left.btsrp.PUNinac(1,:)>0), (ylinemax-0.2)*ones(size(ts(left.btsrp.PUNinac(1,:)>0),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',yel,'Color', yel)
plot(ts(left.btsrp.PUNnp(2,:)<0), (ylinemin-0.3)*ones(size(ts(left.btsrp.PUNnp(2,:)<0),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor','r','Color', 'r')
plot(ts(left.btsrp.PUNnp(1,:)>0), (ylinemax-0.3)*ones(size(ts(left.btsrp.PUNnp(1,:)>0),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor','r','Color', 'r')
 
plot(ts(permPUN.RewAc(1,:)<0.008), (ylinemax-0.3)*ones(size(ts(permPUN.RewAc(1,:)<0.008),2), 2), 's', 'MarkerSize', 5, 'MarkerFaceColor',gree,'Color', gree)
plot(ts(permPUN.RewAc(1,:)<0.008), (ylinemax-0.33)*ones(size(ts(permPUN.RewAc(1,:)<0.008),2), 2), 's', 'MarkerSize', 5, 'MarkerFaceColor',gris,'Color', gris)
plot(ts(permPUN.Rewinac(1,:)<0.008), (ylinemax-0.4)*ones(size(ts(permPUN.Rewinac(1,:)<0.008),2), 2), 's', 'MarkerSize', 5, 'MarkerFaceColor',gree,'Color', gree)
plot(ts(permPUN.Rewinac(1,:)<0.008), (ylinemax-0.43)*ones(size(ts(permPUN.Rewinac(1,:)<0.008),2), 2), 's', 'MarkerSize', 5, 'MarkerFaceColor',yel,'Color', yel)
plot(ts(permPUN.AcInac(1,:)<0.008), (ylinemax-0.5)*ones(size(ts(permPUN.AcInac(1,:)<0.008),2), 2), 's', 'MarkerSize', 5, 'MarkerFaceColor',gris,'Color', gris)
plot(ts(permPUN.AcInac(1,:)<0.008), (ylinemax-0.53)*ones(size(ts(permPUN.AcInac(1,:)<0.008),2), 2), 's', 'MarkerSize', 5, 'MarkerFaceColor',yel,'Color', yel)
 
plot(ts(permPUN.PunRew(1,:)<0.008), (ylinemax-0.6)*ones(size(ts(permPUN.PunRew(1,:)<0.008),2), 2), 's', 'MarkerSize', 5, 'MarkerFaceColor','r','Color', 'r')
plot(ts(permPUN.PunRew(1,:)<0.008), (ylinemax-0.63)*ones(size(ts(permPUN.PunRew(1,:)<0.008),2), 2), 's', 'MarkerSize', 5, 'MarkerFaceColor',gree,'Color', gree)
plot(ts(permPUN.PunAc(1,:)<0.008), (ylinemax-0.7)*ones(size(ts(permPUN.PunAc(1,:)<0.008),2), 2), 's', 'MarkerSize', 5, 'MarkerFaceColor', 'r' ,'Color', 'r' )
plot(ts(permPUN.PunAc(1,:)<0.008), (ylinemax-0.73)*ones(size(ts(permPUN.PunAc(1,:)<0.008),2), 2), 's', 'MarkerSize', 5, 'MarkerFaceColor',gris,'Color', gris)
plot(ts(permPUN.PunInac(1,:)<0.008), (ylinemax-0.8)*ones(size(ts(permPUN.PunInac(1,:)<0.008),2), 2), 's', 'MarkerSize', 5, 'MarkerFaceColor','r','Color', 'r')
plot(ts(permPUN.PunInac(1,:)<0.008), (ylinemax-0.83)*ones(size(ts(permPUN.PunInac(1,:)<0.008),2), 2), 's', 'MarkerSize', 5, 'MarkerFaceColor',yel,'Color', yel)
hold on
jbfill(ts, (mean(left.PUNrew)-(std(left.PUNrew, 0, 1)./sqrt(size(left.PUNrew, 1)))),(std(left.PUNrew, 0, 1)./sqrt(size(left.PUNrew, 1))+mean(left.PUNrew)), gree, gree)
hold on
jbfill(ts, (mean(left.PUNac)-(std(left.PUNac, 0, 1)./sqrt(size(left.PUNac, 1)))),(std(left.PUNac, 0, 1)./sqrt(size(left.PUNac, 1))+mean(left.PUNac)), gris, gris)
hold on
jbfill(ts, (mean(left.PUNinac)-(std(left.PUNinac, 0, 1)./sqrt(size(left.PUNinac, 1)))),(std(left.PUNinac, 0, 1)./sqrt(size(left.PUNinac, 1))+mean(left.PUNinac)), yel, yel)
hold on
jbfill(ts, (mean(left.PUNnp)-(std(left.PUNnp, 0, 1)./sqrt(size(left.PUNnp, 1)))),(std(left.PUNnp, 0, 1)./sqrt(size(left.PUNnp, 1))+mean(left.PUNnp)), 'r', 'r')
legend('Rewarded NP', 'Active not rewarded NP', 'Inactive NP', 'Punished NP')
legend boxoff
box off
set(gca, 'color', 'none')
  vline(0, 'k', 'Nosepoke')
title('PUN sessions left hemisphere, 1-sample bootstrap with alpha 0.0001, perm test alpha 0.008')
    saveas(c, 'C:\Users\Isis\surfdrive\1_Experiments\Nic photometry\Nic14aLeftHemiPUN.fig')
    saveas(c, 'C:\Users\Isis\surfdrive\1_Experiments\Nic photometry\Nic14aLeftHemiPUN.png')
 
%% PUN right
d = figure;
plot(ts, mean(right.PUNrew), 'LineWidth', 2, 'Color', gree)
hold on
plot(ts,mean(right.PUNac), 'LineWidth', 2, 'Color', gris)
plot(ts,mean(right.PUNinac), 'LineWidth', 2, 'Color', yel)
plot(ts,mean(right.PUNpun), 'LineWidth', 2, 'Color', 'r')
% jbfill(ts, btsrp.CSpOne(1,:), btsrp.CSpOne(2,:), 'k', 'k');
line([ts(1) ts(end)], [0 0], 'LineStyle', '--', 'Color', 'k');
ylinemin = max(mean(right.PUNpun, 1));
ylinemax = 1.6*ylinemin;
plot(ts(right.btsrp.PUNrew(2,:)<0), ylinemin*ones(size(ts(right.btsrp.PUNrew(2,:)<0),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',gree,'Color', gree)
plot(ts(right.btsrp.PUNrew(1,:)>0), ylinemax*ones(size(ts(right.btsrp.PUNrew(1,:)>0),2),2), 's', 'MarkerSize', 7, 'MarkerFaceColor',gree,'Color', gree)
hold on
plot(ts(right.btsrp.PUNac(2,:)<0), (ylinemin-0.2)*ones(size(ts(right.btsrp.PUNac(2,:)<0),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',gris,'Color', gris)
plot(ts(right.btsrp.PUNac(1,:)>0), (ylinemax-0.2)*ones(size(ts(right.btsrp.PUNac(1,:)>0),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',gris,'Color', gris)
plot(ts(right.btsrp.PUNinac(2,:)<0), (ylinemin-0.4)*ones(size(ts(right.btsrp.PUNinac(2,:)<0),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',yel,'Color', yel)
plot(ts(right.btsrp.PUNinac(1,:)>0), (ylinemax-0.4)*ones(size(ts(right.btsrp.PUNinac(1,:)>0),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',yel,'Color', yel)
plot(ts(right.btsrp.PUNnp(2,:)<0), (ylinemin-0.6)*ones(size(ts(right.btsrp.PUNnp(2,:)<0),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor','r','Color', 'r')
plot(ts(right.btsrp.PUNnp(1,:)>0), (ylinemax-0.6)*ones(size(ts(right.btsrp.PUNnp(1,:)>0),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor','r','Color', 'r')
plot(ts(permPUN.RewAc(1,:)<0.008), (ylinemax-0.8)*ones(size(ts(permPUN.RewAc(1,:)<0.008),2), 2), 's', 'MarkerSize', 5, 'MarkerFaceColor',gree,'Color', gree)
plot(ts(permPUN.RewAc(1,:)<0.008), (ylinemax-0.83)*ones(size(ts(permPUN.RewAc(1,:)<0.008),2), 2), 's', 'MarkerSize', 5, 'MarkerFaceColor',gris,'Color', gris)
plot(ts(permPUN.Rewinac(1,:)<0.008), (ylinemax-1)*ones(size(ts(permPUN.Rewinac(1,:)<0.008),2), 2), 's', 'MarkerSize', 5, 'MarkerFaceColor',gree,'Color', gree)
plot(ts(permPUN.Rewinac(1,:)<0.008), (ylinemax-1.03)*ones(size(ts(permPUN.Rewinac(1,:)<0.008),2), 2), 's', 'MarkerSize', 5, 'MarkerFaceColor',yel,'Color', yel)
plot(ts(permPUN.AcInac(1,:)<0.008), (ylinemax-1.2)*ones(size(ts(permPUN.AcInac(1,:)<0.008),2), 2), 's', 'MarkerSize', 5, 'MarkerFaceColor',gris,'Color', gris)
plot(ts(permPUN.AcInac(1,:)<0.008), (ylinemax-1.23)*ones(size(ts(permPUN.AcInac(1,:)<0.008),2), 2), 's', 'MarkerSize', 5, 'MarkerFaceColor',yel,'Color', yel)
 
plot(ts(permPUN.PunRew(1,:)<0.008), (ylinemax-1.4)*ones(size(ts(permPUN.PunRew(1,:)<0.008),2), 2), 's', 'MarkerSize', 5, 'MarkerFaceColor','r','Color', 'r')
plot(ts(permPUN.PunRew(1,:)<0.008), (ylinemax-1.43)*ones(size(ts(permPUN.PunRew(1,:)<0.008),2), 2), 's', 'MarkerSize', 5, 'MarkerFaceColor',gree,'Color', gree)
plot(ts(permPUN.PunAc(1,:)<0.008), (ylinemax-1.6)*ones(size(ts(permPUN.PunAc(1,:)<0.008),2), 2), 's', 'MarkerSize', 5, 'MarkerFaceColor', 'r' ,'Color', 'r' )
plot(ts(permPUN.PunAc(1,:)<0.008), (ylinemax-1.63)*ones(size(ts(permPUN.PunAc(1,:)<0.008),2), 2), 's', 'MarkerSize', 5, 'MarkerFaceColor',gris,'Color', gris)
plot(ts(permPUN.PunInac(1,:)<0.008), (ylinemax-1.8)*ones(size(ts(permPUN.PunInac(1,:)<0.008),2), 2), 's', 'MarkerSize', 5, 'MarkerFaceColor','r','Color', 'r')
plot(ts(permPUN.PunInac(1,:)<0.008), (ylinemax-1.83)*ones(size(ts(permPUN.PunInac(1,:)<0.008),2), 2), 's', 'MarkerSize', 5, 'MarkerFaceColor',yel,'Color', yel)
hold on
jbfill(ts, (mean(right.PUNrew)-(std(right.PUNrew, 0, 1)./sqrt(size(right.PUNrew, 1)))),(std(right.PUNrew, 0, 1)./sqrt(size(right.PUNrew, 1))+mean(right.PUNrew)), gree, gree)
hold on
jbfill(ts, (mean(right.PUNac)-(std(right.PUNac, 0, 1)./sqrt(size(right.PUNac, 1)))),(std(right.PUNac, 0, 1)./sqrt(size(right.PUNac, 1))+mean(right.PUNac)), gris, gris)
hold on
jbfill(ts, (mean(right.PUNinac)-(std(right.PUNinac, 0, 1)./sqrt(size(right.PUNinac, 1)))),(std(right.PUNinac, 0, 1)./sqrt(size(right.PUNinac, 1))+mean(right.PUNinac)), yel, yel)
hold on
jbfill(ts, (mean(right.PUNpun)-(std(right.PUNpun, 0, 1)./sqrt(size(right.PUNpun, 1)))),(std(right.PUNpun, 0, 1)./sqrt(size(right.PUNpun, 1))+mean(right.PUNpun)), 'r', 'r')
legend('Rewarded NP', 'Active not rewarded NP', 'Inactive NP', 'Punished NP')
legend boxoff
box off
set(gca, 'color', 'none')
  vline(0, 'k', 'Nosepoke')
title('PUN sessions right hemisphere, 1-sample bootstrap with alpha 0.0001, perm test alpha 0.008')
    saveas(d, 'C:\Users\Isis\surfdrive\1_Experiments\Nic photometry\Nic14aRightHemiPUN.fig')
    saveas(d, 'C:\Users\Isis\surfdrive\1_Experiments\Nic photometry\Nic14aRightHemiPUN.png')
    
    %% SA LEFT
e = figure;
plot(ts, mean(left.TESTSrew), 'LineWidth', 2, 'Color', gree)
hold on
plot(ts,mean(left.TESTSac), 'LineWidth', 2, 'Color', gris)
plot(ts,mean(left.TESTSinac), 'LineWidth', 2, 'Color', yel)
% jbfill(ts, btsrp.CSpOne(1,:), btsrp.CSpOne(2,:), 'k', 'k');
line([ts(1) ts(end)], [0 0], 'LineStyle', '--', 'Color', 'k');
ylinemin = max(mean(left.TESTSrew, 1));
ylinemax = 1.6*ylinemin;
plot(ts(left.btsrp.TESTSrew(2,:)<0), ylinemin*ones(size(ts(left.btsrp.TESTSrew(2,:)<0),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',gree,'Color', gree)
plot(ts(left.btsrp.TESTSrew(1,:)>0), ylinemax*ones(size(ts(left.btsrp.TESTSrew(1,:)>0),2),2), 's', 'MarkerSize', 7, 'MarkerFaceColor',gree,'Color', gree)
hold on
plot(ts(left.btsrp.TESTSac(2,:)<0), (ylinemin-0.1)*ones(size(ts(left.btsrp.TESTSac(2,:)<0),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',gris,'Color', gris)
plot(ts(left.btsrp.TESTSac(1,:)>0), (ylinemax-0.1)*ones(size(ts(left.btsrp.TESTSac(1,:)>0),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',gris,'Color', gris)
plot(ts(left.btsrp.TESTSinac(2,:)<0), (ylinemin-0.2)*ones(size(ts(left.btsrp.TESTSinac(2,:)<0),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',yel,'Color', yel)
plot(ts(left.btsrp.TESTSinac(1,:)>0), (ylinemax-0.2)*ones(size(ts(left.btsrp.TESTSinac(1,:)>0),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',yel,'Color', yel)
 
plot(ts(permTESTS.RewAc(1,:)<0.01), (ylinemax-0.3)*ones(size(ts(permTESTS.RewAc(1,:)<0.01),2), 2), 's', 'MarkerSize', 5, 'MarkerFaceColor',gree,'Color', gree)
plot(ts(permTESTS.RewAc(1,:)<0.01), (ylinemax-0.33)*ones(size(ts(permTESTS.RewAc(1,:)<0.01),2), 2), 's', 'MarkerSize', 5, 'MarkerFaceColor',gris,'Color', gris)
plot(ts(permTESTS.Rewinac(1,:)<0.01), (ylinemax-0.4)*ones(size(ts(permTESTS.Rewinac(1,:)<0.01),2), 2), 's', 'MarkerSize', 5, 'MarkerFaceColor',gree,'Color', gree)
plot(ts(permTESTS.Rewinac(1,:)<0.01), (ylinemax-0.43)*ones(size(ts(permTESTS.Rewinac(1,:)<0.01),2), 2), 's', 'MarkerSize', 5, 'MarkerFaceColor',yel,'Color', yel)
plot(ts(permTESTS.AcInac(1,:)<0.01), (ylinemax-0.5)*ones(size(ts(permTESTS.AcInac(1,:)<0.01),2), 2), 's', 'MarkerSize', 5, 'MarkerFaceColor',gris,'Color', gris)
plot(ts(permTESTS.AcInac(1,:)<0.01), (ylinemax-0.53)*ones(size(ts(permTESTS.AcInac(1,:)<0.01),2), 2), 's', 'MarkerSize', 5, 'MarkerFaceColor',yel,'Color', yel)
hold on
jbfill(ts, (mean(left.TESTSrew)-(std(left.TESTSrew, 0, 1)./sqrt(size(left.TESTSrew, 1)))),(std(left.TESTSrew, 0, 1)./sqrt(size(left.TESTSrew, 1))+mean(left.TESTSrew)), gree, gree)
hold on
jbfill(ts, (mean(left.TESTSac)-(std(left.TESTSac, 0, 1)./sqrt(size(left.TESTSac, 1)))),(std(left.TESTSac, 0, 1)./sqrt(size(left.TESTSac, 1))+mean(left.TESTSac)), gris, gris)
hold on
jbfill(ts, (mean(left.TESTSinac)-(std(left.TESTSinac, 0, 1)./sqrt(size(left.TESTSinac, 1)))),(std(left.TESTSinac, 0, 1)./sqrt(size(left.TESTSinac, 1))+mean(left.TESTSinac)), yel, yel)
legend('Rewarded NP', 'Active not rewarded NP', 'Inactive NP')
legend boxoff
box off
set(gca, 'color', 'none')
  vline(0, 'k', 'Nosepoke')
title('TEST sessions left hemisphere, 1-sample bootstrap with alpha 0.0001, perm test alpha 0.01')
    saveas(e, 'C:\Users\Isis\surfdrive\1_Experiments\Nic photometry\Nic14aLeftHemiTESTS.fig')
    saveas(e, 'C:\Users\Isis\surfdrive\1_Experiments\Nic photometry\Nic14aLeftHemiTESTS.png')
 
%% SA right
f = figure;
plot(ts, mean(right.TESTSrew), 'LineWidth', 2, 'Color', gree)
hold on
plot(ts,mean(right.TESTSac), 'LineWidth', 2, 'Color', gris)
plot(ts,mean(right.TESTSinac), 'LineWidth', 2, 'Color', yel)
% jbfill(ts, btsrp.CSpOne(1,:), btsrp.CSpOne(2,:), 'k', 'k');
line([ts(1) ts(end)], [0 0], 'LineStyle', '--', 'Color', 'k');
ylinemin = max(mean(right.TESTSrew, 1));
ylinemax = 1.6*ylinemin;
plot(ts(right.btsrp.TESTSrew(2,:)<0), ylinemin*ones(size(ts(right.btsrp.TESTSrew(2,:)<0),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',gree,'Color', gree)
plot(ts(right.btsrp.TESTSrew(1,:)>0), ylinemax*ones(size(ts(right.btsrp.TESTSrew(1,:)>0),2),2), 's', 'MarkerSize', 7, 'MarkerFaceColor',gree,'Color', gree)
hold on
plot(ts(right.btsrp.TESTSac(2,:)<0), (ylinemin-0.1)*ones(size(ts(right.btsrp.TESTSac(2,:)<0),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',gris,'Color', gris)
plot(ts(right.btsrp.TESTSac(1,:)>0), (ylinemax-0.1)*ones(size(ts(right.btsrp.TESTSac(1,:)>0),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',gris,'Color', gris)
plot(ts(right.btsrp.TESTSinac(2,:)<0), (ylinemin-0.2)*ones(size(ts(right.btsrp.TESTSinac(2,:)<0),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',yel,'Color', yel)
plot(ts(right.btsrp.TESTSinac(1,:)>0), (ylinemax-0.2)*ones(size(ts(right.btsrp.TESTSinac(1,:)>0),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',yel,'Color', yel)
plot(ts(permTESTS.RewAc(1,:)<0.01), (ylinemax-0.3)*ones(size(ts(permTESTS.RewAc(1,:)<0.01),2), 2), 's', 'MarkerSize', 5, 'MarkerFaceColor',gree,'Color', gree)
plot(ts(permTESTS.RewAc(1,:)<0.01), (ylinemax-0.33)*ones(size(ts(permTESTS.RewAc(1,:)<0.01),2), 2), 's', 'MarkerSize', 5, 'MarkerFaceColor',gris,'Color', gris)
plot(ts(permTESTS.Rewinac(1,:)<0.01), (ylinemax-0.4)*ones(size(ts(permTESTS.Rewinac(1,:)<0.01),2), 2), 's', 'MarkerSize', 5, 'MarkerFaceColor',gree,'Color', gree)
plot(ts(permTESTS.Rewinac(1,:)<0.01), (ylinemax-0.43)*ones(size(ts(permTESTS.Rewinac(1,:)<0.01),2), 2), 's', 'MarkerSize', 5, 'MarkerFaceColor',yel,'Color', yel)
plot(ts(permTESTS.AcInac(1,:)<0.01), (ylinemax-0.5)*ones(size(ts(permTESTS.AcInac(1,:)<0.01),2), 2), 's', 'MarkerSize', 5, 'MarkerFaceColor',gris,'Color', gris)
plot(ts(permTESTS.AcInac(1,:)<0.01), (ylinemax-0.53)*ones(size(ts(permTESTS.AcInac(1,:)<0.01),2), 2), 's', 'MarkerSize', 5, 'MarkerFaceColor',yel,'Color', yel)
hold on
jbfill(ts, (mean(right.TESTSrew)-(std(right.TESTSrew, 0, 1)./sqrt(size(right.TESTSrew, 1)))),(std(right.TESTSrew, 0, 1)./sqrt(size(right.TESTSrew, 1))+mean(right.TESTSrew)), gree, gree)
hold on
jbfill(ts, (mean(right.TESTSac)-(std(right.TESTSac, 0, 1)./sqrt(size(right.TESTSac, 1)))),(std(right.TESTSac, 0, 1)./sqrt(size(right.TESTSac, 1))+mean(right.TESTSac)), gris, gris)
hold on
jbfill(ts, (mean(right.TESTSinac)-(std(right.TESTSinac, 0, 1)./sqrt(size(right.TESTSinac, 1)))),(std(right.TESTSinac, 0, 1)./sqrt(size(right.TESTSinac, 1))+mean(right.TESTSinac)), yel, yel)
legend('Rewarded NP', 'Active not rewarded NP', 'Inactive NP')
legend boxoff
box off
set(gca, 'color', 'none')
  vline(0, 'k', 'Nosepoke')
title('TEST sessions right hemisphere, 1-sample bootstrap with alpha 0.0001, perm test alpha 0.01')
    saveas(f, 'C:\Users\Isis\surfdrive\1_Experiments\Nic photometry\Nic14aRightHemiTESTS.fig')
    saveas(f, 'C:\Users\Isis\surfdrive\1_Experiments\Nic photometry\Nic14aRightHemiTESTS.png')
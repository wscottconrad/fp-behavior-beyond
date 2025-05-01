clear all
filePath1 = 'C:\Users\conrad\Desktop\2024_12_18-03_14_54';
fileName1 = 'Fluorescence.csv';
myData = readtable(strcat(filePath1, '\', fileName1));
timestamps = myData{:,1};
greensignal = myData{:,4};figure;
p1= plot(myData{:,1}, greensignal)

doricdata = open('\\vs03.herseninstituut.knaw.nl\VS03-CSF-1\Heimel\Photometry\Experiments\233508\fluorescentslide\2024-12-18\Fiberphoto\t00001\fiberphotometry.mat');
data = doricdata.data;

expStart = 16000;
expEnd = 25000;
signalD = data(expStart:expEnd,1);

figure;
plot((1:length(data(:,1))),data(:,1));
hpThreshold = 0.0904;
lpThreshold = 0.0597;
figure;
p2= plot(1:length(signalD),signalD(:))
slideSignal = signalD(signalD>hpThreshold);
offSignal = signalD(signalD<lpThreshold); 


b1 = bar(1, mean(offSignal))
hold on;
% sem1 = std(offSignal)./sqrt(length(offSignal(:,1)));
errorbar(1, mean(offSignal), std(offSignal))
b2 = bar(2, mean(slideSignal))
hold on;
% sem2 = std(slideSignal)./sqrt(length(slideSignal(:,1)));
errorbar(2, mean(slideSignal), std(slideSignal))

function [controlfit] = controlfit (signal, isos)

reg = polyfit(isos, signal, 1);
 
a = reg(1);
b = reg(2);
 
if a < 0
   warning('Regression coefficient is negative - changed to 1');
   a = 1;
   b = mean(signal) - mean(isos);
end
 
controlfit = a.*isos(:) + b;
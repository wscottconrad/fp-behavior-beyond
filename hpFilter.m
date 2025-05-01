function [hpf_dats, mov_ave] = hpFilter(ts, ddf, ma)
%% highpass filtering
% from Isis Alonso-Lozares maybe
hp_window = ma; % Convolution window size in secs
plot_hp = [0 0.6 1]; % color vector [# # #])

ts_beg = ts(1); ts_end = ts(end);

% for hp = 1:length(ddf) %why loop? hp is never used
[~,i_wide] = max(ts_end - ts_beg); % gives me 1?
ref_midpoint = max((ts_end(i_wide) + ts_beg(i_wide))/2);
hp_dpts = nb_datapoints(hp_window, ts, ref_midpoint);
%hp_window is in seconds
% end
 
      % Make sure window size is even (necessary for multi_hpf)
      if mod(hp_dpts(1),2) ~= 0
         hp_dpts(1) = hp_dpts(1)+1;
         fprintf('Window size odd, increased by 1.\n');
      end
% dat.Tr_on = ts_beg; dat.Tr_off = ts_end; dat.timestamps = tsFilt;
win_size = hp_dpts(1);
half_wind = win_size/2;
data_leng = length(ddf);
 
      % linear fit to ends of df_f
 
         %if section > than window, use tails of section to extrapolate
         p_begin = polyfit((1:half_wind)', ddf(1:half_wind), 1);
         p_end = polyfit((1:half_wind)', ddf(end-half_wind+1:end), 1);
         
         % linear extrapolation 
         extrap_begin = p_begin(1).*((1-half_wind):0)'+p_begin(2);
         extrap_end = p_end(1).*(half_wind+1:(-1+half_wind*2))'+p_end(2);
         
         % implement moving average window using extrapolations
         mov_ave = zeros(size(ddf));
         mov_ave = conv([extrap_begin;ddf;extrap_end],ones(1,win_size)./win_size,'valid');
 if length(mov_ave) < length(ddf)
    mov_ave(end+1) = mov_ave(end);
end
         hpf_dats = ddf - mov_ave;
 
end
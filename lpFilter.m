function lp_data = lpFilter(data,samp_rate, lowpass_cutoff, filt_steepness, db_atten)
%lpFilter lowpass filter function and spectral analysis, adapted from
%PJRB by Isis Alonso-Lozares 02-2021
 
%    %% DeltaF/F spectral density
%    fft_dff = fft(data);
%    amp = abs(fft_dff);
%    n = length(amp);
%    f = (0:(n-1))/(n/samp_rate);
%    idx = f < max_hz;
 
 %% Low-pass 
      lp_data = lowpass(data,lowpass_cutoff, samp_rate,...
         'ImpulseResponse','iir','Steepness',filt_steepness, 'StopbandAttenuation',db_atten);
end
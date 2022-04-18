function [smoothedData] = butterLowZero(order, cutoff, sampleRate, x)
% [smoothedData] = butterLowZero(order, cutoff, sampleRate, x)
% Applies a zero-lag lowpass buttworth filter with a cutoff at the specified frequency ('cutoff', Hz)
% to data ('x', may be multiple columns/rows. Will filter along longest axis, so orientation should not matter)
% at specified order ('order') for specified sampling rate ('sampleRate')
[rows, cols] = size(x);

x = squeeze(x);

if cols > rows
    x = x'; %transpose into column format
end


% fill gaps, if any exist 
fillmissing(x, 'linear');

[B,A] = butter(order,cutoff/(sampleRate*0.5),'low');


smoothedData = filtfilt(B,A,x);

smoothedData = smoothedData'; %put everything back where you found it <3 
function [ax h] = createHist(data, hBins)

iBin = [1:hBins];
h = hist(data, hBins);
% value of the difference between bins
delta = (max(data) - min(data))/(hBins);
ax = ((min(data) + iBin.*delta) + (min(data) + (iBin-1).*delta))/2;

%skip all values, which are not presented in histogram
ind = find(h);
h = h(ind);
ax = ax(ind);
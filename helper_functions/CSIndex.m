% function [CSI, CumN] = CSIndex(spk1, spk2, MaxT, Nb)
function Logic = CSIndex(spk1, spk2, MaxT, Nb)

%CSIndex calculates the MU synchronization with reference MU.
% a correlation corrologram is constructed and the cumsum of it is then
% calculated to comput a correlation index
%INPUT:
% Spk1 is the time of firing event of ref MU
% Spk2 is time of firing event of a MU
% maxT is the max window for the sync calculation, e.g. 20 ms
% Nb is the number of bins
% By Xiaogang at PennState, June 16, 2009.

Logic = 0;
Syn = [];
for i = 1:length(spk1)
    RefT = spk1(i);
    NS2 = spk2 - RefT;% timing of spike 2 relative to a spike in ref MU 
    In1 = find(NS2 > -1 * MaxT & NS2 < 1 * MaxT); % check if it falls in the window
    if isempty(In1)~=1
    Syn = [Syn, NS2(In1(1))];% if does, catch that firing event relative timing
    end
end

[N, xbin] = hist(Syn, Nb); % calculate histgram of relative timing
maxN = max(N);
PeakCenter = xbin(N==maxN);
In2 = find(Syn<PeakCenter(1)+0.001&Syn>PeakCenter(1)-0.001);
CommonPercent = length(In2)/length(spk1);
 if CommonPercent>0.5
     Logic = 1;
 end
%MN = 
% CumN = cumsum(N - mean(N));




function Logic = CSIndex(spk1, spk2, MaxT, Nb)
%CSINDEX Calculates the motor unit synchronization with reference MU.
%   A correlation correlogram is constructed and the cumsum of it is then
%   calculated to compute a correlation index.
%
%   INPUT:
%   spk1 - Time of firing events of reference MU
%   spk2 - Time of firing events of comparison MU
%   MaxT - Maximum window for the sync calculation (e.g., 20 ms)
%   Nb   - Number of bins
%
%   Author: Xiaogang at PennState
%   Date: June 16, 2009

Logic = 0;
Syn = [];
for i = 1:length(spk1)
    RefT = spk1(i);
    NS2 = spk2 - RefT;  % Timing of spike 2 relative to a spike in ref MU 
    In1 = find(NS2 > -1*MaxT & NS2 < 1*MaxT);  % Check if it falls in the window
    if ~isempty(In1)
        Syn = [Syn, NS2(In1(1))];  % If it does, record the relative timing
    end
end

[N, xbin] = hist(Syn, Nb);  % Calculate histogram of relative timing
maxN = max(N);
PeakCenter = xbin(N == maxN);
In2 = find(Syn < PeakCenter(1) + 0.001 & Syn > PeakCenter(1) - 0.001);
CommonPercent = length(In2)/length(spk1);
if CommonPercent > 0.5
    Logic = 1;
end

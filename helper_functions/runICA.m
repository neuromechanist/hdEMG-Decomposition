function [source,B,spkie_train] = runICA(emg, M, max_iter)
%runICA Rubs (fast) ICA on the hdEMG data and produces the spike trains.
%
%   INPUT:
%   'emg' : The extended hd-EMG that is the out
%
%   'M' : Number of sources.
%
%   'max_iter': The maximum iterations per source to converge.
%
%   OUTPUT:
%   'spike_train':
%
%   'B': unmixing matrix
%
%   'source': The source signal
%
%   REV:
%   v0 @ 09/14/2022
%
%   Copyright (c) 2022 Seyed Yahya Shirazi, shirazi@ieee.org
%% initialize
% the default numbers should be provided, no need to initialize.
tolerance = 10e-4;
% For the matrix manipulation it is easier to have the channels as rows
[num_chan, frames] = size(emg');
B = zeros(num_chan,1);
spkie_train = zeros(frames, M);

for i = 1:M
    w = [];
    w(:,1) = randn(num_chan,1);
    w(:,2) = randn(num_chan,1);
    for n = 2:max_iter
        if abs(w(:,n)'*w(:,n-1)-1)>tolerance
            A = mean(2*w(:,n)'*EMG);
            w(:,n+1) = EMG*(((w(:,n)'*EMG)').^2)-A*w(:,n);
            w(:,n+1) = w(:,n+1) - B*B'*w(:,n+1);
            w(:,n+1) = w(:,n+1)/norm(w(:,n+1));
        else
            break;
        end
    end
    source = w(:,n)'*EMG; % This is the source signals.

    [pks,loc] = findpeaks(source(:,i).^2);
    [idx,~] = kmeansplus(pks',2);
    if sum(idx==1)<=sum(idx==2)
        spike_loc = loc(idx==1);
    else
        spike_loc = loc(idx==2);
    end
    spkie_train(spike_loc,i) = 1;
    B(:,i) = w(:,end);
end



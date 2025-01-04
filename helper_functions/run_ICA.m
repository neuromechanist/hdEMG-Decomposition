function [source,B,spike_train,score] = run_ICA(emg, M, max_iter)
%RUN_ICA Runs (fast) ICA on the hdEMG data and produces spike trains.
%
%   INPUT:
%   'emg' : The extended hd-EMG data matrix
%
%   'M' : Number of sources to extract
%
%   'max_iter': Maximum number of iterations per source for convergence
%
%   OUTPUT:
%   'spike_train': The spike train from the identified peaks for each
%   source (i.e., motor units)
%
%   'B': Unmixing matrix
%
%   'source': The source signal matrix
%
%   'score': Silhouette score measuring the effectiveness of the
%   clustering performed after running ICA to find the peaks
%
%   REV:
%   v0 @ 09/14/2022
%   v1 @ 09/19/2022: Added silhouette scores here 
%
%   Copyright (c) 2022 Seyed Yahya Shirazi, shirazi@ieee.org
%% initialize
% the default numbers should be provided, no need to initialize.
tolerance = 10e-4;
% For the matrix manipulation it is easier to have the channels as rows
emg = emg';
[num_chan, frames] = size(emg);
B = zeros(num_chan,1);
spike_train = zeros(frames, M);
source = zeros(frames,M);
score = zeros(1,M);
fprintf('running ICA for %d sources \n',M)
for i = 1:M
    w = [];
    w(:,1) = randn(num_chan,1);
    w(:,2) = randn(num_chan,1);
    for n = 2:max_iter
        if abs(w(:,n)'*w(:,n-1)-1)>tolerance
            A = mean(2*w(:,n)'*emg);
            w(:,n+1) = emg*(((w(:,n)'*emg)').^2)-A*w(:,n);
            w(:,n+1) = w(:,n+1) - B*B'*w(:,n+1);
            w(:,n+1) = w(:,n+1)/norm(w(:,n+1));
        else
            break;
        end
    end
    source(:,i) = w(:,n)'*emg; % Extract the source signal

    [pks,loc] = findpeaks(source(:,i).^2);
    [idx,~] = kmeansplus(pks',2);
    sil_score = silhouette(pks,idx);
    score(i) = (mean(sil_score(idx==1))+mean(sil_score(idx==2)))/2;
    if sum(idx==1)<=sum(idx==2)
        spike_loc = loc(idx==1);
    else
        spike_loc = loc(idx==2);
    end
    spike_train(spike_loc,i) = 1;
    B(:,i) = w(:,end);
    fprintf('.')
end
spike_train = sparse(spike_train);
disp("\n ICA decomposition is completed")

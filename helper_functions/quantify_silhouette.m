function score = quantify_silhouette(source, freq)
%quantify_silhouette() Rubs (fast) ICA on the hdEMG data and produces the spike trains.
%
%
%   Similar to the silhouette scoring for all kmeans clustering problems,
%   the aim is to understand how good the clustering was performed in the
%   runICA() method
%   INPUT:
%   'source' : The source signa, T x N, where T is the frames.
%
%   'freq' : signal frequency.
%
%   OUTPUT:
%   'score': The Silhouette score for the source.
%
%   REV:
%   v0 @ 09/19/2022
%
%   Copyright (c) 2022 Seyed Yahya Shirazi, shirazi@ieee.org
%% initialize
[b,a] = butter(4,500/(freq/2),'low'); % low-pass filter for the source signals
score = zeros(1,size(source,2)); % motor units are at the columns, frames are at the rows.


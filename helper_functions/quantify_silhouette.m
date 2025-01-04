function score = quantify_silhouette(source, freq)
%QUANTIFY_SILHOUETTE Calculate silhouette scores for motor unit clusters.
%   Similar to silhouette scoring for k-means clustering problems,
%   this function evaluates the quality of clustering performed by the
%   run_ICA() method.
%
%   INPUT:
%   'source' : The source signal, T x N matrix, where T is the number of frames
%             and N is the number of sources.
%
%   'freq'   : Sampling frequency in Hz.
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
source = filtfilt(b,a,source);
score = zeros(1,size(source,2)); % Initialize scores (columns = motor units, rows = frames)

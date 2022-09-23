function [spike_train,source, good_idx] = remove_motorUnit_duplicates(spike_train, source, freq)
%REMOVE_MOTORUNIT_DUPLICATES removes the duplicate sources based identifed
%by the ICA algoritm.
%
%
%   INPUTS:
%   'spike_train' :  the signal containing the spikes detected by the
%   runICA() for each motor unit.
%
%   'source' : the source singal, also the output of the runICA() function
%
%   'freq' :  The EMG recording frequency. Default = 2048
%
%   OUTPUTS:
%   'spike_train' : This time the duplicates are removed.
%
%   'source' : The duplicates are removed from the source as well.
%
%   'good_idx':  The source indecies that are good, ie, pass this test
%
%   REV:
%   v0 @ 09/15/2022
%
%   Copyright (c) 2022 Seyed Yahya Shirazi, shirazi@ieee.org
%% initialize
if ~exist("freq","var") || isempty(freq), freq = 2048; end % default value for the recoding frequency 
% The minimum of muscle firing is set to 4 Hz and the max should not be
% more than 35 Hz. This insight is form the physiological stnadpoint. For
% refrence look at Winter's biomechanics book chapter 9 and 10.
min_firing = 4; % in Hz
max_firing = 35; % in Hz
min_firing_interval = 1/max_firing; % in miliseconds. This equals to max_firning = 50Hz
time_stamp = linspace(1/freq,length(spike_train)/freq,length(spike_train));

%% step 1, select physiologically plausible sources
firings = sum(spike_train,1);
lowerBound_cond = find(firings>min_firing*time_stamp(end));
upperBound_cond = find(firings<max_firing*time_stamp(end));
plausible_firings = intersect(upperBound_cond,lowerBound_cond); % these are the sources with plausible firings

%% step2, select only one spike from twin spikes
% It so happens that the spikes in a good source with "playsible firgins"
% are very close to each other (say < 20ms, i.e., 50Hz spike rate), Such
% spike rates are not physiological, therefore, only one can exist. We only
% take the one with greater peak
for k = plausible_firings
    spike_timeDiff = diff(time_stamp(spike_train(:,k)==1));
    for t = 1:length(spike_timeDiff)
        if spike_timeDiff(t) < min_firing_interval
            if source(t,k) < source(t+1,k)
                spike_train(t,k) = 0;
            else
                spike_train(t+1,k) = 0;
            end
        end
    end
end

%% step3, finally, the duplicates motor units should be removed
% CSIndex is a fucntion to find the simialrity of source's (i.e., motor
% unit) spike trains. It works based on the overlapping histograms.
max_timeDiff = 0.01; % in seconds (as the time_stamp is in seconds) 
num_bins = 10; % number of bins for the histogram
duplicate_sources = [];
for k = plausible_firings
    if ~ismember(k,duplicate_sources)
        for j = setdiff(plausible_firings(plausible_firings~=k),duplicate_sources)
            isduplicate  = CSIndex(time_stamp(spike_train(:,k)==1),time_stamp(spike_train(:,j)==1),...
                max_timeDiff, num_bins);
            if isduplicate
                duplicate_sources = [duplicate_sources j]; %#ok<AGROW> 
            end
        end
    end
end
good_idx = setdiff(plausible_firings,duplicate_sources);
spike_train = sparse(spike_train(:,good_idx));
source = source(:,good_idx);




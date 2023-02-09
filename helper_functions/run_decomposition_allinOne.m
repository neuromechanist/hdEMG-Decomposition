function motor_unit = run_decomposition(varargin)
%RUN_DECOMPOSITION      motor-unit decomposition on hdEMG datasets
%
%   This function and the helper files are mainly a reorganization of the
%   code accompanied with Hyser Dataset by Jian et. al. <a href=
%   'https://doi.org/10.1109/TNSRE.2021.3082551'>Hyser Dataset</a>
%   The original code and the dataset is also <a href=
%   'https://doi.org/10.13026/ym7v-bh53'>available at PhysioNet</a>
%
%   Requirements:
%           Matlab R2021b+, Statistics and Machine Learning Toolbox, DSP
%           Toolbox, and Communications Toolbox (to inject white noise).
%           This function uses tables, with features available from
%           R2017b+.
%
%   INPUT:
%   Inputs are name pairs:
%
%   'data' : (The path to the) hdEMG data. if 'data' is pointing to the
%   lcation of the data file, the data file must be a MAT array.
%   Default is the sample file included in the toolbox.
%
%   'data_mode' : EMG can be recorded in the 'monopolar' or 'bipolar'
%   mode. Default = 'monopolar'
%
%   'frq','sampling_frequency' : Sampling frequency of the data, default = 2048 Hz
%   
%   'R', 'extension_parameter' : The number of times to repeat the data
%   blocks, see the Hyser paper for more detail. Default = 4
%   
%   'M', 'max_iter' : Maximum iterations for the (FAST) ICA decompsition.
%   Default = 300
%
%   'whiten_flag' : Whether to whiten the data prior to the ICA. Default =
%   1
%
%   'SNR', 'inject_noise' : Adding white noise to the EMG mixutre. Uses
%   Communication Toolbox AWGN fucntion. Default = Inf for not injecting
%   any artificial noise.
%
%   'SIL_thresh' : The silhouette threshold to detect the good motor units.
%   Default = 0.6
%
%   'save_path' : The path that the files should be saved there. The
%       function does not create the path, rather uses it. Default is the
%       is the 'sample' path of the toolbox.
%
%   'save_flag' : Whether the files are saved or not, default is 0, so it is
%       NOT saving your output.
%
%   OUTPUT:
%   'motor_unit': The structure inclduing the follwing fields:
%       spike_train
%       waveform
%       ica_weight
%       whiten_matrix
%       SIL
%
%   EXAMPLE:
%   sample_motorUnit = RUN_DECOMPOSITION()
%
%   REV:
%   v0 @ 09/12/2022
%
%
%   Copyright (c) 2022 Seyed Yahya Shirazi, shirazi@ieee.org

%% initialize
fs = filesep;  % a shortcut to the filesep 
addpath('./helper_functions')
opts = arg_define(varargin, ...
        arg({'data'},['sample_data' fs 'sample1'] ,[],'The file is also from the Hyser dataset'), ...
        arg({'data_mode'}, 'monopolar',['monopolar','bipolar'],'Data mode of the hd-EMG data'), ...
        arg({'frq','sampling_frequency'}, 2048,[],'Sampling frequency of the imported dataset'), ...
        arg({'R', 'extension_parameter'}, 4,[],'The number of times to repeat the data blocks'), ...
        arg({'M', 'max_sources'}, 300,[],'Maximum number of sources being decomposed by (FAST) ICA.'), ...
        arg({'whiten_flag'}, 1,[0,1],'Whether to whiten the data prior to the ICA.'), ...
        arg({'SNR', 'inject_noise'}, Inf,[],'WAdding white noise to the EMG mixutre.'), ...
        arg({'SIL_thresh'}, 0.6,[],'The silhouette threshold to detect the good motor units.'), ...
        arg({'output_file'},['sample_data' fs 'sample1_decomposed'] ,[],'The path to the output file.'), ...
        arg({'save_flag'}, 1,[0,1],'Whether the files are saved or not.'),...
        arg({'plot_spikeTrain', 'plot_fig'}, 1,[0,1],'Whether to plot the resulting spike trains.'),...
        arg({'load_ICA'}, 0,[0,1],'Whether to load precomputed ICA results for debugging.'));
data = opts.data;
if ischar(data)
    emg_file = load(data); % Based on the strcutre of the exproted MAT files from OTBiolab+ v1.5.8
    data = emg_file.Data{1}; % Assuming that the file ONLY contains one hdEMG array!
    opts.frq = emg_file.SamplingFrequency;
end
max_iter = 200;

%% run the decomposition
[extended_emg,~] = emg_preprocess(data,'R',opts.R,'WhitenFlag',opts.whiten_flag,'SNR',opts.SNR);
if ~opts.load_ICA
    [uncleaned_source,B,uncleaned_spkieTrain,score] = run_ICA(extended_emg, opts.M, max_iter); % B is the unmixing matrix
else
    warning("Loading ICA results from a saved file. Change the 'load_ICA' flag if you want to run ICA.")
    ica_results = load([opts.data '_ica_results']);
    uncleaned_source = ica_results.uncleaned_source; uncleaned_spkieTrain = ica_results.uncleaned_spkieTrain;
    score = ica_results.score; B = ica_results.B;
end
[spike_train, source, good_idx] = remove_motorUnit_duplicates(uncleaned_spkieTrain, uncleaned_source, opts.frq);
silhouette_score = score(good_idx);

%% save the results
% motor_unit.raw.B = B; motor_unit.raw.extended_emg = extended_emg;  % this creates a very large file thay would not fit on github.   
% motor_unit.raw.uncleaned_source = uncleaned_source; motor_unit.raw.score = score;
motor_unit.spike_train = spike_train; motor_unit.source = source; motor_unit.good_idx = good_idx;
motor_unit.silhouette_score = silhouette_score;
if opts.save_flag, save(opts.output_file,"-struct","motor_unit"); end

%% plot the motor units in a nice way
minScore_toPlot = 0.9; % The minimum silhouette score of the motor units to be included in the plot
if opts.plot_spikeTrain, plot_spikeTrain(spike_train,opts.frq,silhouette_score,minScore_toPlot); end
end

%% Auxiliary fucntions
function [preprocessed_data,W] = emg_preprocess(data,varargin)
%EMG_PREPROCESS Prepares the emg array for decomposition.
%
%   INPUT:
%   'data_mode' : EMG can be recorded in the 'monopolar' or 'bipolar'
%   mode. Default = 'monopolar'
%
%   'whiten_flag' : Whether to whiten the data prior to the ICA. Default =
%   1
%
%   'SNR': Adding white noise to the EMG mixutre. Uses Communication Toolbox
%    AWGN fucntion. Default = Inf for not injecting any artificial noise.
%
%   'R': The number of times to repeat the data blocks, see the Hyser paper
%    for more detail. Default = 4
%
%   'array_shape': The first element will be used to calcualte the bipolar
%   activitiy if the bipolar flag is on for the 'data_mode'.
%
%   REV:
%   v0 @ 09/12/2022
%
%
%   Copyright (c) 2022 Seyed Yahya Shirazi, shirazi@ieee.org

%% initialize
data_mode = 'monopolar';
whiten_flag = 1;
R = 4;
SNR = Inf;
array_shape = [8,8];

for i = 1:2:length(varargin)
    switch varargin{i}
        case 'SNR'
            SNR = varargin{i+1};
        case 'whiten_flag'
            whiten_flag = varargin{i+1};
        case 'R'
            R = varargin{i+1};
        case 'data_mode'
            data_mode = varargin{i+1};
        case 'array_shape'
            array_shape = varargin{i+1};
    end
end
% data can come in coloumn or row format, but needs to become the coloumn
% format where the 
num_chan = min(size(data));  % Let's assume that we have more than 64 frames
if num_chan ~= size(data,2), data = data'; end

%% preprocessing
% Add white noise
if ~isinf(SNR), data = awgn(data,SNR,'dB'); end
% create a bipolar setting from the monopolar data
if string(data_mode) == "bipolar"
    for i = 1:num_chan-array_shape(1)
        data(:,i) = data(:,i) - data(:,i+array_shape(1));
    end
    data(:,end-8:end) = [];
end

extended_data = zeros(size(data,1),size(data,2)*(R+1));
extended_data(:,1:size(data,2)) = data; % to make a consistent downstream
if R~=0
    for i = 1:R
        % This basically shift the data on the replications one step
        % forward. This is pretty standard in the ICA, as ICA should be
        % able to parse the sources out pretty well. Also, it introduces
        % small delay, R/freq, which with R=64, delay = 64/2048= 31ms.
        % This addition reinforces finding MUAPs, despite having
        % duplicates. Later on, the duplicates will be removed.
        extended_data(1+i:end,size(data,2)*i+1:size(data,2)*i+size(data,2)) = ...
            data(1:end-i,:);
    end
end

if whiten_flag
    [whitened_data, ~, ~, W] = whiten(extended_data);
    preprocessed_data = whitened_data;
else
    preprocessed_data = data;
end
disp("EMG preprocessing is completed")
end

function [Xwh, mu, invMat, whMat] = whiten(X,epsilon)
%function [X,mu,invMat] = whiten(X,epsilon)
%
% ZCA whitening of a data matrix (make the covariance matrix an identity matrix)
%
% WARNING
% This form of whitening performs poorly if the number of dimensions are
% much greater than the number of instances
%
%
% INPUT
% X: rows are the instances, columns are the features
% epsilon: small number to compensate for nearly 0 eigenvalue [DEFAULT =
% 0.0001]
%
% OUTPUT
% Xwh: whitened data, rows are instances, columns are features
% mu: mean of each feature of the orginal data
% invMat: the inverse data whitening matrix
% whMat: the whitening matrix
%
% EXAMPLE
%
% X = rand(100,20); % 100 instance with 20 features
% 
% figure;
% imagesc(cov(X)); colorbar; title('original covariance matrix');
% 
% [Xwh, mu, invMat, whMat] = whiten(X,0.0001);
% 
% figure;
% imagesc(cov(Xwh)); colorbar; title('whitened covariance matrix');
% 
% Xwh2 = (X-repmat(mean(X), size(X,1),1))*whMat; 
% figure;
% plot(sum((Xwh-Xwh2).^2),'-rx'); title('reconstructed whitening error (should be 0)');
% 
% Xrec = Xwh*invMat + repmat(mu, size(X,1),1);
% figure;
% plot(sum((X-Xrec).^2),'-rx'); ylim([-1,1]); title('reconstructed data error (should be zero)');
%
% Author: Colorado Reed colorado-reed@uiowa.edu
%
% Copyright (c) 2012, Colorado Reed
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

if ~exist('epsilon','var')
    epsilon = 0.0001;
end

mu = mean(X); 
X = bsxfun(@minus, X, mu);
A = X'*X;
[V,D,~] = svd(A);
whMat = sqrt(size(X,1)-1)*V*sqrtm(inv(D + eye(size(D))*epsilon))*V';
Xwh = X*whMat;  
invMat = pinv(whMat);
end


function [source,B,spike_train,score] = run_ICA(emg, M, max_iter)
%run_ICA() Runs (fast) ICA on the hdEMG data and produces the spike trains.
%
%   INPUT:
%   'emg' : The extended hd-EMG that is the out
%
%   'M' : Number of sources.
%
%   'max_iter': The maximum iterations per source to converge.
%
%   OUTPUT:
%   'spike_train': The spike train from the identified peaks for each
%   source (ie, motor units)
%
%   'B': unmixing matrix
%
%   'source': The source signal
%
%   'score': This is the silhouette score for the effectiveness of the
%   clusterting performed after the running the ICA to find the peakes.
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
    source(:,i) = w(:,n)'*emg; % This is the source signals.

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
end

function [L,C] = kmeansplus(X,k)
%KMEANS Cluster multivariate data using the k-means++ algorithm.
%   [L,C] = kmeans(X,k) produces a 1-by-size(X,2) vector L with one class
%   label per column in X and a size(X,1)-by-k matrix C containing the
%   centers corresponding to each class.

%   Version: 2013-02-08
%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%
%   References:
%   [1] J. B. MacQueen, "Some Methods for Classification and Analysis of 
%       MultiVariate Observations", in Proc. of the fifth Berkeley
%       Symposium on Mathematical Statistics and Probability, L. M. L. Cam
%       and J. Neyman, eds., vol. 1, UC Press, 1967, pp. 281-297.
%   [2] D. Arthur and S. Vassilvitskii, "k-means++: The Advantages of
%       Careful Seeding", Technical Report 2006-13, Stanford InfoLab, 2006.

L = [];
L1 = 0;

while length(unique(L)) ~= k
    
    % The k-means++ initialization.
    C = X(:,1+round(rand*(size(X,2)-1)));
    L = ones(1,size(X,2));
    for i = 2:k
        D = X-C(:,L);
        D = cumsum(sqrt(dot(D,D,1)));
        if D(end) == 0, C(:,i:k) = X(:,ones(1,k-i+1)); return; end
        C(:,i) = X(:,find(rand < D/D(end),1));
        [~,L] = max(bsxfun(@minus,2*real(C'*X),dot(C,C,1).'));
    end
    
    % The k-means algorithm.
    while any(L ~= L1)
        L1 = L;
        for i = 1:k, l = L==i; C(:,i) = sum(X(:,l),2)/sum(l); end
        [~,L] = max(bsxfun(@minus,2*real(C'*X),dot(C,C,1).'),[],1);
    end
    
end
end

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
end
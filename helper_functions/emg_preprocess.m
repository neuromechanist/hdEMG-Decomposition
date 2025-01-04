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
%   'SNR': Adding white noise to the EMG mixture. Uses Communication Toolbox
%    AWGN function. Default = Inf for not injecting any artificial noise.
%
%   'R': The number of times to repeat the data blocks, see the Hyser paper
%    for more detail. Default = 4
%
%   'array_shape': The first element will be used to calculate the bipolar
%   activity if the bipolar flag is on for the 'data_mode'.
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
            if length(array_shape) ~= 2
                error("array_shape should be strictly a 1 x 2 vector")
            elseif mod(array_shape(1) * array_shape(2), 2) == 1
                error("array_shape elements should be even numbers. Please choose another array shape"+...
                    "\n if [13, 5] is used, please rerun the script with [8, 8].")
            end
    end
end
% Data can come in column or row format, but needs to be in column
% format for processing.
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
extended_data(:,1:size(data,2)) = data; % Initialize first block for consistent downstream processing
if R~=0
    for i = 1:R
        % This shifts the data on the replications one step
        % forward. This is standard practice in ICA, as ICA should be
        % able to parse the sources out effectively. It introduces a
        % small delay of R/freq (e.g., with R=64, delay = 64/2048 = 31ms).
        % This addition reinforces finding MUAPs, despite having
        % duplicates which will be removed later.
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

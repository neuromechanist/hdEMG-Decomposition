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
%   'recording_mode' : EMG can be recorded in the 'monopolar' or 'bipolar'
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
%   Copyright 2022 Seyed Yahya Shirazi, shirazi@ieee.org

%% initialize
fs = filesep;  % a shortcut to the filesep 
addpath('./helper_functions')
opts = arg_define(varargin, ...
        arg({'data'},[fs 'sample_data' fs 'sample.mat'] ,[],'The file is also from the Hyser dataset'), ...
        arg({'recording_mode'}, 'monopolar',['monopolar','bipolar'],'Recroding mode of the hd-EMG data'), ...
        arg({'frq','sampling_frequency'}, 2048,[],'Sampling frequency of the imported dataset'), ...
        arg({'R', 'extension_parameter'}, 4,[],'The number of times to repeat the data blocks'), ...
        arg({'M', 'max_iter'}, 300,[],'Maximum iterations for the (FAST) ICA decompsition.'), ...
        arg({'whiten_flag'}, 1,[0,1],'Whether to whiten the data prior to the ICA.'), ...
        arg({'SNR', 'inject_noise'}, Inf,[],'WAdding white noise to the EMG mixutre.'), ...
        arg({'SIL_thresh'}, 0.6,[],'The silhouette threshold to detect the good motor units.'), ...
        arg({'save_path'},[fs 'sample_data' fs] ,[],'The path that the files should be saved there.'), ...
        arg({'save_flag'}, 0,[0,1],'Whether the files are saved or not.'));
data = opts.data;
if ischar(data), data = load(data); end

%% run the decomposition
[emg_extend,W] = SimEMGProcessing(data,'R',opts.R,'WhitenFlag',opts.whiten_flag,'SNR',opts.SNR);

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

function plot_spikeTrain(spike_train,frq,sil_score,minScore_toPlot)
%PLOT_SPIKETRAIN plots the spike trains in a sorted manner
%   The physiology of the motor unit activations suggests that the motor
%   units with more spikes should have been activated earlier than the motor
%   units with fewer spikes. This fact, however, should be verified by
%   finding the thresholds of each motor unit.
%   PLOT_SPIKETRAIN plots the spike trains above the min silhouette
%   threshold over time in a vertically sorted manner.
%
%   INPUTS:
%   'spike_train': The array of the spikes, it is a FRAME x MOTOR UNIT
%   matrix.
%
%   'frq': The signal frequency.
%
%   'sil_score': The silhouette score for each motor unit.
%
%   'minScore_toPlot': The minimum silhouette score of the motor units that
%   should be plotted.
%
%   REV:
%   v0 @ 09/21/2022
%
%
%   Copyright (c) 2022 Seyed Yahya Shirazi, shirazi@ieee.org
%% initialize
if ~exist("minScore_toPlot","var") || isempty(minScore_toPlot), minScore_toPlot = 0.7; end % Default value for the recording frequency 
selected_spikeTrain = spike_train(:,sil_score>minScore_toPlot);
[~,order] = sort(sum(selected_spikeTrain,1),"descend");
bar_height = 0.2;
figure("Renderer","painters","Name","Spike trains of the motor units");
sub_handle = subplot(1,1,1);
ylim(sub_handle,[0,10])
hold on
%% plot as a rug plot
% Plotting will be very similar to the rug plots I used to have
% (https://github.com/neuromechanist/add_rug_plot). However, here the rugs
% are the plot.
rugVal = num2cell(selected_spikeTrain(:,order)',2);
rugCount = length(rugVal);
ylim(sub_handle,[0,rugCount*bar_height])
hold on
for r = 1: rugCount
    signRugVal = rugVal{r}; % only one row of rug plot at a time.
    signRugVal = abs(sign(signRugVal)); % rug plots only accepts 0 and 1.
    if exist('barColor', 'var') && size(barColor,1)==rugCount
        bCol = barColor(r,:);
    else
        barColor = turbo(rugCount); bCol = barColor(r,:);
    end
    if r == 1
        plot(sub_handle,[find(signRugVal);find(signRugVal)],[zeros(1,length(find(signRugVal)));...
            ones(1,length(find(signRugVal)))* bar_height*0.7],'Color',bCol,"LineWidth",0.01)
    elseif ~isempty(find(signRugVal)) %#ok<EFIND>
        plot(sub_handle,[find(signRugVal);find(signRugVal)],[ones(1,length(find(signRugVal)))*bar_height*(r-1);...
            ones(1,length(find(signRugVal)))*(bar_height*(r-1)+bar_height*0.7)],'Color',bCol,"LineWidth",0.01)
    end
end
xlim(sub_handle,[0,length(selected_spikeTrain)])
xticks(linspace(0,length(selected_spikeTrain),10));
xticklabels(round(linspace(0,length(selected_spikeTrain)/frq,10)));
xlabel("time (sec)")

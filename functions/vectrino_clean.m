
function [u_i,w2] = vectrino_clean(u_i_raw,t,fs,SNR,corr,amp,opts,show)

% PURPOSE - Post-processing code for Vectrino ADV velocity profiles
%
% This script removes points with "bad" correlation, either every point
% with correlation below 70% or the 10% of points with the worst
% correlation, whichever results in the removal of the fewest points. A
% similar approach is applied with SNR <= 15dB.
%
% In addition, despiking is performed on the beam velocities.  However,
% a high pass filter applied prior to despiking, reducing variance due to
% low frequency oscillations and increasing the effectiveness of the
% despiking algorithm.
%
% Arguments IN:: 
% u_i_raw = raw velocity profiles
% t = time vector
% fs = sampling freq. (Hz)
% SNR = signal to noise ratio array
% corr = correlation array
% amp = amplitude array
% opts = options
%   0: (default) replaces bad points with NaN
%   1: Removes bad points
%   2: Replaces bad points with polynomial interpolation via neighbours
% show = show velocity plot or not (1 = Yes)
%
% Arguments OUT:
% u_i = a cell array of length n, where n is the number of profiler bins.
% Each cell is a t by 3 array, where t is the number of time samples.
% The second dimension is the velocity component
% w2 = vector of 2nd vertical velocity
%
% Version control:
% 03/06/2022 - SNR histogram removed. Figures renumbered. Figure formatting settings removed.

%% SETTINGS

% Baseline raw data
u_i = u_i_raw;
u_i(:,4) = [];
w2 = u_i_raw(:,4);

% User prompt for how to deal with 'bad' points
% prompt = ['Which option (0,1,2)? \n' ...
%     '0: (default) replaces bad points with NaN \n' ...
%     '1: Removes bad points \n' ...
%     '2: Replaces bad points with polynomial interpolation via neighbours \n'];
% opts = input(prompt);
opts = 2;

%% FIND BAD POINTS BASED ON CORRELATION

% NOTE: If a point is bad in one beam, all beams are considered bad, e.g. @ t = x, corr logic = 001 (1 in Z direction)... time step is denoted as 'bad'
corr(:,end) = []; % remove w2
SNR_Logic = ones(size(corr));
SNR_Logic(corr<70) = 0;
corr_Logic2 = prod(SNR_Logic,2); % product of array elements. 2 returns column vector containing products of each row
bad_corr = zeros(length(corr_Logic2),1);

if sum(corr_Logic2 == 0) < 0.1*length(corr_Logic2)
    % Remove data points less than 70
    bad_corr(corr_Logic2 == 0,:) = 1;
else
    % Remove the worst 10% of data points
    [~,idx_corr] = sort(corr); % sort vector in ascending order. Indices at the top represent the lowest corr values
    idx_corr = reshape(idx_corr',numel(idx_corr),1); % reshape(X,[M,N]) returns the M-by-N matrix whose elements are taken columnwise from X.
    % sort_corr = reshape(sort_corr',numel(idx_corr),1);
    unique_idx_corr = unique(idx_corr,'stable'); % returns the unique values in a specific order w/o repetitions (noting original array had three columns of data, which were then sorted and indices returned, hence the duplication. Array length = original array / 3). 'Stable' implies no sorting
    % unique_sort_corr = unique(sort_corr,'stable');
    remove_Idx = unique_idx_corr(1:round(length(corr_Logic2)*0.1)); % removes indices equal to first 10% of dataset, noting that the indices represent the corr values sorted in ascending order, i.e. lowest corr values first
    bad_corr(remove_Idx) = 1; % indexing lowest 10% of corr values, and set to 1 where values should be removed
end

%% FIND BAD POINTS BASED ON SNR

% Same logic as corr
SNR(:,end) = []; % remove w2
SNR_Logic = ones(size(SNR));
SNR_Logic(SNR<15) = 0;
SNR_Logic2 = prod(SNR_Logic,2);
bad_SNR = zeros(length(SNR_Logic2),1);

if sum(SNR_Logic2 == 0) < 0.1*length(SNR_Logic2)
    % Remove data points less than 15
    bad_SNR(SNR_Logic2 == 0,:) = 1;
else
    % Remove the worst 10% of data points
    [~,idx_SNR] = sort(SNR); 
    idx_SNR = reshape(idx_SNR',numel(idx_SNR),1); 
    % sort_SNR = reshape(sort_SNR',numel(idx_SNR),1);
    unique_idx_SNR = unique(idx_SNR,'stable');
    % unique_sort_SNR = unique(sort_SNR,'stable');
    remove_Idx = unique_idx_SNR(1:round(length(SNR_Logic2)*0.1));
    bad_SNR(remove_Idx) = 1;
end

%% HIGH-PASS FILTER

% 'designfilt' generates a filter based on frequency-response specifications 
cutoffFreq = 3; % A high-pass filter passes signals with a freq > cutoff freq and attenuates signals with freq lower than the cutoff frequency. May need to adjust
dif = designfilt('highpassfir', 'FilterOrder', 20, 'CutoffFrequency', cutoffFreq, 'SampleRate', fs);
% highpassfir = Finite Impulse Response (FIR) is a filter whose impulse response is of finite period, as a result, it settles to zero in finite time.
% High pass filter to remove long-period fluctuations.
% The rate at which a filter's response falls in the transition band is determined by the filter's order; the higher the order, the faster its rolloff rate.

%% FIND BAD POINTS BASED ON DESPIKING

% Preallocation
bad_Spike = zeros(length(corr_Logic2),1); % one column whose length indicates times series of data
bad_Spike_w2 = zeros(length(corr_Logic2),1); % reserved for w2

% NOTE: If a point is bad in one beam, all beams are considered bad.
% y = filtfilt(b,a,x) applies zero-phase digital filter 'filtfilt by processing input data x in both the forward and reverse directions
% After filtering the data in the forward direction, the function reverses the filtered sequence and runs it back through the filter

for i = 1:width(u_i_raw)
    f = filtfilt(dif,detrend(u_i_raw(:,i))); % apply filter to velocity dataset. Velocity detrended by removing mean value
    [~,spike_idx] = func_despike_phasespace3d_wahl(f); % returns spike indices
    if i == 4
        bad_Spike_w2(spike_idx) = 1; % accumulates indices
    else
        bad_Spike(spike_idx) = 1;
    end
end

%% MERGE INDICES IDENTIFIED FROM 'BAD' POINTS DUE TO CORRELATION, SNR OR DESPIKING

% User input required to define criteria for removing 'bad' points
% prompt = 'Choose removal of data points due to correlation, SNR, and/or despiking in concatenated binary format XYZ, noting 0 = do not apply filter; 1 = apply filter: \n';
% removal = input(prompt,'s');
removal = '001';
if removal(1) == '0'
    bad_corr = zeros(size(bad_corr));
    if removal(2) == '0'
        bad_SNR = zeros(size(bad_SNR));
        if removal(3) == '0'
            bad_Spike = zeros(size(bad_Spike));
        end
    end
elseif removal(2) == '0'
    bad_SNR = zeros(size(bad_SNR));
    if removal(3) == '0'
        bad_Spike = zeros(size(bad_Spike));
    end
elseif removal(3) == '0'
    bad_Spike = zeros(size(bad_Spike));
    bad_Spike_w2 = zeros(size(bad_Spike_w2));
end

% Applies to Vx, Vy and Vz1
bad_idx = zeros(size(bad_Spike));
bad_idx(bad_Spike | bad_corr | bad_SNR) = 1;
bad_idx = logical(bad_idx);

% Applies to Vx, Vy and Vz2
bad_idx_w2 = zeros(size(bad_Spike_w2));
bad_idx_w2(bad_Spike_w2 | bad_corr | bad_SNR) = 1;
bad_idx_w2 = logical(bad_idx_w2);

%% APPLY OPTIONS TO MANAGING THE BAD DATA POINTS

% Replace with NaN (0), remove bad points (1) or replace with polynomial interpolation (2)
    
if opts == 1
    u_i{bad_idx,:} = [];
elseif opts == 0
    u_i{bad_idx,:} = NaN;
elseif opts == 2
    % Applies to Vx, Vy, Vz1 (u_i) array
    if bad_idx(1)
        % first good index
        fgind = find(~bad_idx,1,'first')+1; % Look for the last nonzero element in array
        u_i(1,:) = u_i(fgind,:);
        bad_idx(1) = 0;
    end
    if bad_idx(end)
        % last good index
        bgind = find(~bad_idx,1,'last')+1; % Look for the last non-zero element in array
        u_i(end,:) = u_i(bgind,:);
        bad_idx(end) = 0;
    end
    x = find(~bad_idx); % returns indices of each non-zero element of good points, i.e. not outliers
    xq = find(bad_idx); % returns indices of each non-zero element of outliers
    for i = 1:3
        u_i(xq,i) = interp1(x,u_i(x,i),xq,'PCHIP'); % 1-D linear interpolation where x = good indices, u_i = corresponding good values, xq = coordinates of query points, and 'pchip' is shape-preserving piecewise cubic interpolation (third-order polynomial, through 12 points on either side of the spike.
    end

    fprintf('\n>> Bad points due to corr. < 70, SNR < 15, and spikes = %d + %d + %d = %0.2f%%\n',sum(bad_corr),sum(bad_SNR),sum(bad_Spike),(100*(sum(bad_corr)+sum(bad_SNR)+sum(bad_Spike))/length(u_i)))

    % Applies to w2 array
    if bad_idx_w2(1)
        % first good index
        fgind = find(~bad_idx_w2,1,'first')+1; % Look for the last nonzero element in array
        w2(1,:) = w2(fgind,:);
        bad_idx_w2(1) = 0;
    end
    if bad_idx_w2(end)
        % last good index
        bgind = find(~bad_idx_w2,1,'last')+1; % Look for the last non-zero element in array
        w2(end,:) = w2(bgind,:);
        bad_idx_w2(end) = 0;
    end
    x = find(~bad_idx_w2); % returns indices of each non-zero element of good points, i.e. not outliers
    xq = find(bad_idx_w2); % returns indices of each non-zero element of outliers
    w2(xq) = interp1(x,w2(x),xq,'PCHIP'); % 1D linear interpolation where x = good indices, u_i = corresponding good values, xq = coordinates of query points, and 'pchip' is shape-preserving piecewise cubic interpolation (third-order polynomial, through 12 points on either side of the spike.

    % fprintf('>> 1D velocity (w2) - Bad points due to corr. < 70, SNR < 15, and despiking = %d + %d + %d = %d ...... this equates to %0.2f%% of all points \n\n',sum(bad_corr),sum(bad_SNR),sum(bad_Spike_w2),sum(bad_corr)+sum(bad_SNR)+sum(bad_Spike_w2),(100*(sum(bad_corr)+sum(bad_SNR)+sum(bad_Spike_w2))/length(w2)))
   
end

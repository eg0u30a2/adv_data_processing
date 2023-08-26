
function  [u_i,t,t_ref,order,w2,SNR,corr,amp] = merging(fs,ques_profile,filepath,filelist,params,test_id)

% PURPOSE - Merge profiles along an axis and define the order of concatenation.
%
% Arguments IN:
% fs = sampling freq. (Hz)
% ques_profile = definition of profile axis
% filepath = filepaths for each of the datsets that comprise the profile
% filelist = filenames for each of the datsets that comprise the profile
% params = test parameters
% test_id = test ID
%
% Arguments OUT:
% u_i = 2D velocity array including w1
% t = 1D time array
% t_ref = cell array of individual datasets (used for data windowing)
% order = order of concatenation
% w2 = 1D velocity array of w2
% SNR = 2D SNR array 
% corr = 2D correlation array 
% amp = 2D amplitude array 
% 
% Version control:
% 03/06/2022 - params added as an input
% 22/09/2022 - throw error is no. datasets < required
% 19/01/2023 - Run through vectrino_clean.m twice

%% LOAD DATA 

% Preallocation
data_cell = num2cell(zeros(1,height(test_id)));

for i = 1:height(test_id)
    data_cell{i} = importdata(filepath(i,:),' '); % import data using space delimiter
    data_cell{i} = array2table(data_cell{i}); % convert data to a table
end

%% CONCATENATE DATA FILES AND CONVERT

% Reorder dataset to reflect profile
if ques_profile == 'R' | ques_profile == 'Z'
    [~,order] = sort(params.Z1);
elseif ques_profile == 'X'
    [~,order] = sort(params.X);
end

data = cell(1,height(test_id)); % preallocation
for i = 1:length(data)
    data{i} = data_cell{order(i)};
end

%% ASSIGN VARIABLES AND TABLE HEADERS

% Assign headers to each dataset
for i = 1:width(data) % headers
    data{i}.Properties.VariableNames = {'X1','Time','Step','X2','Vx','Vy','Vz1','Vz2','Amp1','Amp2','Amp3','Amp4','SNR1','SNR2','SNR3','SNR4','Corr1','Corr2','Corr3','Corr4'};
end

% Assign variables of interest
[t,u_i_raw,SNR,corr,amp] = cellfun(@variables,data,'UniformOutput', 0);

%% REVIEW ADV METRICS AND APPLY DESPIKING

% Preallocation
u_i = cell(1,width(data)); % for despiked velocities

% Do you wish to view the metrics for each segment of a profile/transect?
prompt = ['\nDo you wish to view the ADV metrics for each segment? Yes (Y) or No (N) \n>> '];
show = input(prompt,'s');
if show == 'Y' || show == 'y'
    show = 1;
elseif show == 'N' || show == 'n'
    show = 0;
else
end
opts = 2;

for i = 1:width(data)
    a = u_i_raw{i};
    b = t{i};

    % Sampling frequency (Hz) - in case time series does not indicate 200 Hz
    if ( 1 / (b(2,1) - b(1,1)) ~= 200)
        b(:) = [0:(1/fs):((length(b)-1)/fs)]'; % replace spurious time series data with fs
    else
    end

    fprintf('\nDataset %d \n------------',i)
    
    [u_i{i},w2{i}] = vectrino_clean(a,b,fs,SNR{i},corr{i},amp{i},opts,show);

    % Run through cleaning for a second time
    fprintf('\nDataset %d for a second time\n----------------------------',i)
    a = u_i{i};
    a(:,end+1) = w2{i};
    [u_i{i},w2{i}] = vectrino_clean(a,b,fs,SNR{i},corr{i},amp{i},opts,show);

end

%% TERMINATE ANALYSIS IF PLAUSIBILITY CHECKS FAIL

if params.Z1 == -0.5 & params.Z2 <= 0.3
    fprintf('\nThis is likely to be a near-boundary profiling exercise. Please check.\n')
elseif width(t) < 4
    prompt = [ '\nThis is less than the required no. datasets for a profile... Is this:\n' ...
    '- a calibration point, sanity check, or boundary layer analysis (1)?\n' ...
    '- an error (2)?\n\n>> '];
    ques = str2double(input(prompt,'s'));
    if ques == 1
        fprintf('\nOK proceed.\n')
    elseif ques == 2
        error('** You need more time series files. A minimum of four is required. **')
    end
else
end

clear prompt ques

%% CONCATENATE CELLS INTO DOUBLE ARRAYS

corr = cat(1, corr{:});
SNR = cat(1, SNR{:});
amp = cat(1, amp{:});

%% DEVELOP TIME SERIES

len = cellfun(@height, t, 'UniformOutput', false); % length of time series of each dataset
t_steps = sum([len{:}]); % total time series length for traverse
t_ref = t;
t = (0:1/fs:(t_steps-1)/fs)'; % time vector converted to table

%% CONCATENATE VELOCITY DATA INTO ONE DATASET, READY FOR EXPORT

% Clean data
u_i = cat(1, u_i{:});

% w2
w2 = cat(1, w2{:});

close all
end
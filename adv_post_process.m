
% PURPOSE: Main script for ADV data post-processing and preparing the datasets for analysis

% INSTRUCTIONS
% User to add own filepaths, otherwise errors will result
% Press Run and follow instructions in the Command Window
% User defines the nature of the measurement, i.e. single point only (S), streamwise profile of single points (X), vertical profile (Z), cross-stream profile (Y), or transect (R), and calls the linked functions for ADV data post-processing
% The script requires ADV input data in ASCII format (the Nortek) Vectrino software has the option to output .vna files.

% Version control:
% 01/02/2022 - Initial release
% 07/10/2022 - Update static pitch offset on ADV mount rev1 (JMC) based on inclinometer measurement, not calipers
% 12/10/2022 - Update static pitch offsets for ADV mounts rev1 and rev2 based on inclinometer measurements, and accounting for both orientations
% 24/10/2022 - Accommodate transect analysis
% 23/12/2022 - Add provision for streamwise profiling (single points only)
% 26/01/2023 - Pitch offset updated following investigation (0169-0176)
% 27/01/2023 - Removal of pitch rotation tolerancing
% 21/03/2023 - Removal of measured pitch offsets - this is now calculated parametrically

clc ; clearvars; close all;

%% SETUP

warning('off')

%% LOAD TEST MATRIX DOCUMENT

data_dir = '...enter filepath here';
file = 'Test Matrix.xlsx';
path = fullfile(data_dir,file);
opts = detectImportOptions(path,'DataRange','A2','VariableNamesRange','A1');
params = readtable(path,opts);

%% DEFINE PROPERTIES

% Flume
w_f = 1.1; % width of the flume (m)
h_f = 1; % water depth (m)

% Fluid
mu = 0.0010518; % dynamic viscosity (Pa.s) @ 18 degC
nu = 1.0533E-6; % kinematic viscosity (m^2/s) @ 18 degC

% Measurement
fs = 200; % sampling frequency (Hz)

%% DEFINE THE NATURE OF THE ANALYSIS

% Screen 1
prompt = ['Are you processing:\n' ...
    '- a single point (1)\n' ...
    '- a profile (2)\n' ...
    '- a transect (3)? \n>> '];
ques_anal = str2double(input(prompt,'s'));

%% IDENTIFY FILE LOCATIONS

subfolderInfo = dir(cd); % set directory
subfolderNames = {subfolderInfo.name}; % list all folder names within the directory
subfolderNames(ismember(subfolderNames,{'.','..'}))=[]; % remove fields named '.' or '..'
filelist = cell(1,numel(subfolderNames)); % prellocation for a list of files
filepath = cell(1,numel(subfolderNames)); % preallocation for a list of filepaths

%% IDENTIFY FILE(S) AND PARAMETERS OF INTEREST

% Screen 2
if ques_anal == 1
    % Isolate file of interest
    for i = 1:numel(subfolderNames)
        fileInfo = dir(subfolderNames{i}); % file info within the directory
        fileNames = {fileInfo.name};
        fileNames(ismember(fileNames,{'.','..'}))=[]; % remove cells named '.' or '..'
        folder = fileInfo.folder;
        for j = 1:numel(fileNames)
            if (fileNames{j}(end-2:end) == 'vna') & (fileNames{j}(1) == 'S')
                filelist{i}{j,1} = fileNames{j}; % create list of file names
                filepath{i}{j,1} = strcat(folder,'\',fileNames{j}); % create list of filepaths
            else
            end
        end
    end
    filelist = cat(1,filelist{:}); % concatenate the arrays into a X x 1 cell array
    filepath = cat(1,filepath{:}); % concatenate the arrays into a X x 1 cell array
    filelist = filelist(~cellfun('isempty',filelist)); % remove blank cells
        % creates table from cell array (purely to see the indices)    
        list = cell2table(filelist);
        list.idx = (1:1:length(filelist))';
        list.Properties.VariableNames = {'File' 'idx'};
        list = movevars(list,'File','After','idx')
    filepath = filepath(~cellfun('isempty',filepath)); % remove blank cells
    
    % Define file to analyse
    prompt = 'Which file would you like to process? \n>> ';
    id = input(prompt); % index
    path = filepath{id};
    test_id = filelist{id};
    data = importdata(path,' '); % import data using space delimiter
    data = array2table(data); % convert data to a table

    % Extract parameters of interest
    idx = find(strcmp(params.Filename,test_id(1:10))); % compare strings to identify parameters, then find idx of non-zero value of a logical array
    params = params(idx,:); %
    ques_profile = [];

elseif ques_anal == 2
    prompt = '\nDefine the nature of the traverse: Is it streamwise (X) cross-stream (Y) or vertical (Z)? \n>> ';
    ques_profile = input(prompt,'s');

    % Streamwise profiling is only applicable for single points
    if ques_profile == 'X'
        fprintf('\nStreamwise profiling at this stage of post-processing is only applicable for single point measurements.\n\n')
    else
    end

    % Scan folder structure for 'profiling' datasets
    for i = 1:numel(subfolderNames)
        fileInfo = dir(subfolderNames{i}); % file info within the directory
        fileNames = {fileInfo.name};
        fileNames(ismember(fileNames,{'.','..'}))=[]; % remove cells named '.' or '..'
        folder = fileInfo.folder;
        for j = 1:numel(fileNames)
            if (fileNames{j}(end-2:end) == 'vna') & (fileNames{j}(1) == 'P') % filter Profiling datasets only
                filelist{i}{j,1} = fileNames{j}; % create list of file names
                filepath{i}{j,1} = strcat(folder,'\',fileNames{j}); % create list of filepaths
            else
            end
        end
    end
    % Produce list of files and filepaths that relate to profiling
    filelist = cat(1,filelist{:}); % concatenate the arrays into a X x 1 cell array
    filepath = cat(1,filepath{:}); % concatenate the arrays into a X x 1 cell array
    filelist = filelist(~cellfun('isempty',filelist)); % remove blank cells
    filepath = filepath(~cellfun('isempty',filepath)); % remove blank cells

    % Identify indices of profiles from params table (includes archived data)
    idx = find(strcmp(params.Profile_Axis,ques_profile)); % compare strings to identify indices of profiles, then find idx of non-zero value of a logical array
    filenames = params.Filename(idx);

    % Identify indices of profiles within 'filelist' and 'filepath'
    idx = zeros(size(filelist,1),size(filenames,1));
    a = char(filelist); % random variable for use in for loop
    a = cellstr(a(:,1:10)); % create cell array from first 10 characters of character array for use in strcmp below
    for i = 1:size(filenames,1)
        idx(:,i) = strcmp(char(filenames(i,:)),a);
    end
    
    idx = logical(sum(idx,2)); % indices of profiles
    filelist = filelist(idx); % ignores files not related to profile of interest
        % creates table from cell array (purely to see the indices)    
        list = cell2table(filelist);
        list.idx = (1:1:length(filelist))';
        list.Properties.VariableNames = {'File' 'idx'};
        list = movevars(list,'File','After','idx')
    filepath = filepath(idx); % ignores files not related to profile of interest

    if ques_profile == 'Y'

        % Define file to analyse
        prompt = 'Which file would you like to process? \n>> ';
        id = input(prompt); % idx = index
        path = filepath{id};
        test_id = filelist{id};
        data = importdata(path,' '); % import data using space delimiter
        data = array2table(data); % convert data to a table
    
        % Extract parameters of interest
        idx = find(strcmp(params.Filename,test_id(1:10))); % compare strings to identify parameters, then find idx of non-zero value of a logical array
        params = params(idx,:); %
        order = [];
    
        % Sampling frequency (Hz) - in case time series does not indicate 200 Hz
        if ( 1 / (data.data2(2) - data.data2(1))) ~= 200
            data.data2 = [0:(1/fs):((length(data.data2)-1)/fs)]'; % replace spurious time series data with fs
        else
        end

    elseif ques_profile == 'X' | ques_profile == 'Z'
    
        % Define file to analyse
        prompt = 'Enter the indices that constitute a complete traverse:\n>> ';
        id = input(prompt); % idx = index
        filepath = char(filepath(id,:));
        filelist = char(filelist(id,:))
        test_id = char(filelist(:,1:10));
    
        % Extract parameters of interest
        clear idx
        for i = 1:height(test_id)
            idx(i) = find(strcmp(params.Filename,test_id(i,1:10))); % compare strings to identify parameters, then find idx of non-zero value of a logical array
        end
        params = params(idx,:); %
    
        % Concatenated time series' are merged and the order of the concatenation is
        [u_i,t,t_ref,order,w2,SNR,corr,amp] = merging(fs,ques_profile,filepath,filelist,params,test_id);

    end

elseif ques_anal == 3
    ques_profile = 'R';
    % Scan folder structure for 'transect' datasets
    for i = 1:numel(subfolderNames)
        fileInfo = dir(subfolderNames{i}); % file info within the directory
        fileNames = {fileInfo.name};
        fileNames(ismember(fileNames,{'.','..'}))=[]; % remove cells named '.' or '..'
        folder = fileInfo.folder;
        for j = 1:numel(fileNames)
            if (fileNames{j}(end-2:end) == 'vna') & (fileNames{j}(1) == 'R') % filter transect datasets only
                filelist{i}{j,1} = fileNames{j}; % create list of file names
                filepath{i}{j,1} = strcat(folder,'\',fileNames{j}); % create list of filepaths
            else
            end
        end
    end
    % Produce list of files and filepaths that relate to profiling (Y or Z)
    filelist = cat(1,filelist{:}); % concatenate the arrays into a X x 1 cell array
    filepath = cat(1,filepath{:}); % concatenate the arrays into a X x 1 cell array
    filelist = filelist(~cellfun('isempty',filelist)); % remove blank cells
    filepath = filepath(~cellfun('isempty',filepath)); % remove blank cells

    % Identify indices of transects from params table (includes archived data)
    temp = params.Filename;
    temp = cellfun (@(x) x(1),temp,'un',0); % 'un' is short for 'uniformoutput'... wraps the result of each function call into a cell
    idx = find(strcmp(temp,ques_profile)); % compare strings to identify indices of transects, then find idx of non-zero value of a logical array
    filenames = params.Filename(idx);
    clear temp

    % Identify indices of transects within 'filelist' and 'filepath'
    idx = zeros(size(filelist,1),size(filenames,1));
    a = char(filelist); % random variable for use in for loop
    a = cellstr(a(:,1:10)); % create cell array from first 10 characters of character array for use in strcmp below
    for i = 1:size(filenames,1)
        idx(:,i) = strcmp(char(filenames(i,:)),a);
    end
    
    idx = logical(sum(idx,2)); % indices of profiles
    filelist = filelist(idx); % ignores files not related to transects of interest
        % creates table from cell array (purely to see the indices)    
        list = cell2table(filelist);
        list.idx = (1:1:length(filelist))';
        list.Properties.VariableNames = {'File' 'idx'};
        list = movevars(list,'File','After','idx')
    filepath = filepath(idx); % ignores files not related to transects of interest
 
    % Define file to analyse
    prompt = 'Enter the indices that constitute a complete transect:\n>> ';
    id = input(prompt); % idx = index
    filepath = char(filepath(id,:));
    filelist = char(filelist(id,:))
    test_id = char(filelist(:,1:10));

    % Extract parameters of interest
    clear idx
    for i = 1:height(test_id)
        idx(i) = find(strcmp(params.Filename,test_id(i,1:10))); % compare strings to identify parameters, then find idx of non-zero value of a logical array
    end
    params = params(idx,:);

    % Concatenated time series' are merged and the order of the concatenation is established
    [u_i,t,t_ref,order,w2,SNR,corr,amp] = merging(fs,ques_profile,filepath,filelist,params,test_id);

else
    error('Try again.')
end

%% ASSIGN VARIABLES AND HEADERS

if ques_anal == 1 || (ques_anal == 2 && ques_profile == 'Y')
    [t,u_i_raw,SNR,corr,amp] = variables(data);
else
end

%% REVIEW ADV METRICS AND APPLY DESPIKING

if ques_anal == 1 || (ques_anal == 2 && ques_profile == 'Y')

    if ques_anal == 1
        opts = 2; show = 1;
    elseif ques_anal == 2 && ques_profile == 'Y'
        opts = 2; show = []; t_ref = [];
    end
    [u_i,w2] = vectrino_clean(u_i_raw,t,fs,SNR,corr,amp,opts,show);

else
end

%% HARMONISE ADV AND FLUME COORDINATE SYSTEMS

% Define transforms, e.g. *-1
if mean(u_i(:,1)) < 0
    v = -1; % identifier for pitch offset later
    u_i(:,1) = u_i(:,1) *-1; % 
    u_i(:,2) = u_i(:,2) *-1; %
else
    v = 1;
end

%% ROTATION

% Apply rotational corrections via pitch and yaw offsets
% DoR = degrees of rotational correction

if char(params.Probe) == 'D' % down-looking probe

    % Evaluate DoR
    % Pitch correction based on centreline velocities == 0 cannot be applied for unsteady shear flows, i.e. the vertical component should remain true

    if isempty(char(params.Unsteady)) | contains(params.Unsteady(1),'-X-')
        DoR = 2; % corrections assume plug flow, i.e. uniform velocity profile with zero yaw and zero pitch offsets along the flume centreline
    elseif contains(params.Unsteady(1),'G-Y')
        DoR = 0; % above assumption does not hold for turbulent flows as the profile is highly non-uniform
    elseif contains(params.Unsteady(1),'S-Y') 
        DoR = 1;
    end
elseif char(params.Probe) == 'S' % side-looking probe
    DoR = 0; % i.e. zero pitch offset

end

% Single measurement points only, and must be along the centreline, with a downwards-looking probe
if ques_anal == 1 && params.Yi == 0 && params.Zi == 0 && test_id(1) ~= 'F' && char(params.Probe) == 'D' && (params.Y1 == 0) && (params.Z1 == 0)
    
    % 2D velocity array including w1
    u_i1_pre_rot = u_i;
    u_i1 = rotate_axes(u_i1_pre_rot,fs,t,DoR);
    
    if char(params.Probe) == 'D' 
        % 2D velocity array including w2
        u_i2_pre_rot = u_i1_pre_rot; % copy array
        u_i2_pre_rot(:,3) = [];
        u_i2_pre_rot = [u_i2_pre_rot w2];
        u_i2 = rotate_axes(u_i2_pre_rot,fs,t,DoR);
        u_i_pre_rot = [u_i1_pre_rot u_i2_pre_rot(:,3)]; % pre-rotation including w1 and w2
        u_i = [u_i1 u_i2(:,3)]; % contains 4 columns, each representing a velocity component with the appropriate rotational correction applied
    else
        u_i_pre_rot = u_i1_pre_rot; % i.e. side probe
        u_i = u_i1;
    end
    clear u_i1 u_i2 u_i1_pre_rot u_i2_pre_rot

elseif test_id(1,1) == 'P' || test_id(1,1) == 'R'
    % Rotational corrections of profiles are managed in the respective post-processing function
else
    % We do not rotate the side-looking probe or data points that are not along the flume centreline
    u_i = [u_i w2];
    DoR = 0;

    clear w2
end

clear data_dir file fileInfo filelist filenames fileNames folder i idx j list opts subfolderInfo subfolderNames v u_i_raw

%% PROFILING - POST-PROCESSING

% Apply rotational corrections, data windowing, and noise floor analysis
if ques_anal == 2 | ques_anal == 3
     [~,Post_process] = profiling_post_process(DoR,fs,params,ques_profile,t,t_ref,u_i,w2);
     clear t u_i w2
else
end

%% NOISE FLOOR

% Single point measurement only
path = []; pos = []; % dummy variables

% Reference velocity array
if ques_anal == 1 && params.Yi == 0 && params.Zi == 0 && ~isnan(params.Mot_percent)

    % Not necessarily a centreline modulus
    U0_mod = mean(sqrt(u_i(:,1).^2 + u_i(:,2).^2 + u_i(:,3).^2));

    if length(u_i) > 3*60*fs % i.e. if the time series is > 3 minutes
        [Spectra] = noise_floor_analysis(u_i,U0_mod,fs,params,path,pos,ques_anal);
    else
    end

else
    % The noise floor analysis of the profiling is included in the respective post-processing function
end

clear U0_mod_centre_bar

%% SAVE PROCESSED DATA FILES

% Repository for saving dataset
parent_path = '...add your own path here';
if ques_profile == 'Z' | ques_profile == 'R' | ques_profile == 'X'
    target = split(filepath(1,:),'\');
else
    target = split(filepath{id,:},'\');
end
target = char(target{13});

path = strcat(parent_path,'\',target,'\','Processed');

% Create character string
if size(params,1) == 1
    file = sprintf('%04.0f.mat',params.ID(1));
elseif size(params,1) > 1
    file = sprintf('%04.0f-%04.0f.mat',params.ID(1),params.ID(end));
end

% Save
clear a amp corr data data_dir extrap_offset filepath id parent_path pos prompt ques_anal show static_offset SNR switchover t_ref target w2 zpos
save(fullfile(path,file))

close all
clear; clearvars -global; clc; close('all', 'hidden')

fprintf('Analysis complete.\n\n')

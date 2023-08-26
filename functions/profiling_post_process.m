
function [coords,Post_process] = profiling_post_process(DoR,fs,params,ques_profile,t,t_ref,u_i,w2)

% PURPOSE: post-processing code for profiled measurements. Functionality includes data windowing, rotational correction, and removal of noise floor.
% and rotational correction

% Arguments IN:
% DoR = degrees of rotational correction
% fs = sampling frequency (Hz)
% params = parameters of interest
% ques_profile = Defines axis of traverse
% t = time array (secs)
% t_ref = cell array of individual time series'
% u_i = 2D velocity array with w1 (m/s)
% w2 = 1D velocity array with w2 (m/s)

% Arguments OUT:
% coords = measurement coordinates
% Post_process = data structure featuring initial data analysis for noise evaluation
%
% Version control
% 01/05/2022 - DoR = 2 introduced
% 08/10/2022 - Pitch and yaw corrections combined into one figure. Selected pitch offset reflects static offset, rather than calculated pitch offset from central data window. Section removed and archived in 'Notes'
% 10/10/2022 - Provision added into pre-rotation data windowing ('Y' profile) to enable data window size to be modified in case time series was of insufficient length
% 20/10/2022 - Provision added for 10 s sampling close to boundary. Determination of calibration point in Z profile made more robust using 'any' function
% 24/10/2022 - Provision added for transect profiling
% 13/01/2023 - Provision added for (single-point) streamwise profiling 
% 27/01/2023 - Removal of pitch rotation tolerancing
% 21/05/2023 - U0 reference velocity section added to encompass different options

%% DEFINE COORDINATES FOR TRAVERSING

if ques_profile == 'X'
    a = min(params.X);
    b = abs(mean(diff(params.X)));
    c = max(params.X);
    coord = "[" + a + "-" + c + "]" + "-" + params.Y1(1) + "-" + params.Z1(1); % string for filename
    coords = sort(params.X)';
elseif ques_profile == 'Y'
    a = params.Y1;
    b = params.Yi;
    c = params.Y2;
    coord = params.X + "-[" + a + "-" + c + "]-" + params.Z1; % string for filename
    coords = a:-b:c;
elseif ques_profile == 'Z'
    a = min(round(params.Z1,1));
    b = mode(params.Zi);
    if b == 0
        b = mode(diff(params.Z2));
    else
    end
    c = max(params.Z2);
    coord = params.X(1) + "-" + params.Y1(1) + "-[" + a + "-" + c + "]"; % string for filename
    coords = a:b:c;
    if c ~= coords(end)
        idx = find(abs(params.Z2- params.Z1) - params.Zi >= params.Zi); % find test case where measurement range > Zi
        inc = (params.Z1(idx):params.Zi(idx):params.Z2(idx))'; % incremental array to be injected into coordinates array
        list = [params.Z1;inc;params.Z2]';
        coords = unique(sort(list));
    else
    end
elseif ques_profile == 'R'
    a_z = min(params.Z1);
    b = mode(diff(params.Z1));
    c_z = max(params.Z2);
    a_y = mode(params.Y1);
    b_y = mode(params.Yi);
    c_y = mode(params.Z2);
    % Check dy and dz are equal
    if b ~= b_y
        error('Need to alter the code to accommodate different traverse increments (mm) for Y and Z.')
    else
    end
    coord = params.X(1) + "-[" + c_y + "-" + a_y + "]" + "-[" + a_z + "-" + c_z + "]"; % string for filename
    coords_y = round([c_y:b_y:a_y],3);
    coords_z = a_z:b:c_z;
    clear a_z c_z a_y b_y c_y y_pos
end

%% ESTIMATE PERIOD BETWEEN SIGNAL PEAKS

% Time window of interest (secs)
if params.ID(1) == 13
    p = 265;
elseif b < 0.002 % i.e. boundary layer measurement
    p = 30;
else
    p = 300;
end

% Options for traverse speed and distance between measurements
traverse.vel_trav = mode(params.Trav_V); % traverse speed (m/s)
traverse.acc_trav = 0.001; % accel and decel to target velocity
if p <= 30
    traverse.t_int_nom = 2;
else
    traverse.t_sample = 320; % duration of sampling (secs) at one measurement point
    traverse.s_nom = 0.1; % nominal distance between measurement points. Note that this might not always work if dz and dy in a transect are different
    if ques_profile == 'R'
        traverse.s_nom = b;
    else
    end
    % Nominal period to traverse between measurement points
    traverse.t_acc = traverse.vel_trav/traverse.acc_trav; 
    traverse.s_accel = 0.5*traverse.acc_trav*(traverse.t_acc^2); 
    traverse.s_vmax_nom = traverse.s_nom - (2*traverse.s_accel); % distance travelled at max traverse velocity
    traverse.t_vmax_nom = traverse.s_vmax_nom/traverse.vel_trav; % time elapsed at max traverse velocity
    traverse.t_int_nom = traverse.t_sample+traverse.t_vmax_nom+(2*traverse.t_acc); 
end

%% DATA WINDOWING (1) - PRE-ROTATION

% Parameters
win = 5; % period of interest (secs) for time averaging for 'movmean' function
if ques_profile == 'R'
    meas_points = length(coords_z)*length(coords_y);
elseif params.ID(1) == 583
    meas_points = length(coords)-3; % no. measurement points along the traverse axis
else
    meas_points = length(coords); % no. measurement points along the traverse axis
end

% Preallocation
if ques_profile == 'X' | ques_profile == 'Y' | ques_profile == 'Z'
    x = zeros(1,meas_points); % preallocation for 'x' measurement locations. Represents centre of window
    x_win_a = zeros(1,meas_points); % values for start of data windows
    x_win_b = zeros(1,meas_points); % values for end of data windows
    data_win_pre_rot = cell(1,meas_points); % data from each window to be stored here
elseif ques_profile == 'R'
    x = zeros(1,length(coords_y));
    x_win_a = zeros(1,length(coords_y));
    x_win_b = zeros(1,length(coords_y));
    data_win_pre_rot = cell(length(coords_z),length(coords_y));
end

% Time series with windows
set(0,'DefaultFigureVisible','on')
if ques_profile == 'X'
    
    for i = 1:length(x)
        if i == 1
            x_win_a(i) = 1;
            x_win_b(i) = length(t_ref{i});
        elseif i > 1
            x_win_a(i) = x_win_a(i-1) + length(t_ref{i-1});
            x_win_b(i) = x_win_b(i-1) + length(t_ref{i});
        else
        end
        data_win_pre_rot{i} = u_i(x_win_a(i):x_win_b(i),:);
    end
    clear x
elseif ques_profile == 'Y'
    
    for i = 1:length(x)
        x(i) = 1 + (p/2) + (traverse.t_int_nom*(i-1));
        x_win_a(i) = x(i)-(0.5*p);
        x_win_b(i) = x(i)+(0.5*p);
        if fs*x_win_b(i) > length(u_i)
            delta = fs*x_win_b(i) - length(u_i);
            x_win_b(i) = floor(length(u_i)/fs);
        else
        end
        data_win_pre_rot{i} = u_i((fs*x_win_a(i)) : (fs*x_win_b(i)),:);
    end
elseif ques_profile == 'Z'

    if p >= 300
        for i = 1:width(t_ref)
            if i == 1
                x(2*i - 1) = 1 + (p/2); % p secs = period of interest for averaging
                if width(t_ref) == 1
                    for j = 2:meas_points
                        x(j) = x(j-1) + traverse.t_int_nom;
                    end
                elseif width(t_ref) == 4 & range(params.Z1) >= 0.8
                    x(i+1) = x(2*i - 1) + traverse.t_int_nom;
                else
                    x(i+1) = (height(t_ref{1,1})/fs) + (p/2);   
                    % x(i+1) = x(2*i - 1) + traverse.t_sample; % (note: no allowance for traverse movement at Z = -0.4 and Z = -0.3 as these are separate datasets anyway
                end
            elseif i == 2
                if width(t_ref) == 4 & range(params.Z1) >= 0.8
                    x(2*i - 1) = (height(t_ref{1,1})/fs) + (p/2);
                elseif width(t_ref) > 4 || all(params.Z1 > 0.399)
                    x(2*i - 1) = (height(t_ref{1,1})/fs) + (height(t_ref{1,2})/fs) + (p/2);
                else
                end
                x(2*i) = x(2*i - 1) + traverse.t_int_nom;
                if width(t_ref) > 5 || all(params.Z1 > 0.399)
                    x(2*i) = (height(t_ref{1,1})/fs) + (height(t_ref{1,2})/fs) + (height(t_ref{1,3})/fs) + (p/2);
                else
                end
            elseif i == 3
                if width(t_ref) == 4 & range(params.Z1) >= 0.8
                    x(2*i - 1) = (height(t_ref{1,1})/fs) + (height(t_ref{1,2})/fs) + (p/2);
                elseif width(t_ref) == 5
                    x(2*i - 1) = (height(t_ref{1,1})/fs) + (height(t_ref{1,2})/fs) + (height(t_ref{1,3})/fs) + (p/2);
                elseif width(t_ref) > 5
                    x(2*i - 1) = (height(t_ref{1,1})/fs) + (height(t_ref{1,2})/fs) + (height(t_ref{1,3})/fs) + (height(t_ref{1,4})/fs) + (p/2);
                end
                if width(t_ref) > 5
                    x(2*i) = (height(t_ref{1,1})/fs) + (height(t_ref{1,2})/fs) + (height(t_ref{1,3})/fs) + (height(t_ref{1,4})/fs) + (height(t_ref{1,5})/fs) + (p/2);
                elseif all(params.Z1 > 0.399)
                    % No action
                else
                    x(2*i) = x(2*i - 1) + traverse.t_int_nom;
                    x(2*i + 1) = x(2*i) + traverse.t_int_nom;
                end
            elseif i == 4
                if width(t_ref) == 4 & range(params.Z1) >= 0.8
                    x(2*i) = (height(t_ref{1,1})/fs) + (height(t_ref{1,2})/fs) + (height(t_ref{1,3})/fs) + (p/2);
                    x(2*i + 1) = x(2*i) + traverse.t_int_nom;
                elseif width(t_ref) == 5
                    x(2*i) = (height(t_ref{1,1})/fs) + (height(t_ref{1,2})/fs) + (height(t_ref{1,3})/fs) + (height(t_ref{1,4})/fs) + (p/2);
                    x(2*i + 1) = x(2*i) + traverse.t_int_nom;
                elseif width(t_ref) > 5
                    x(2*i - 1) = (height(t_ref{1,1})/fs) + (height(t_ref{1,2})/fs) + (height(t_ref{1,3})/fs) + (height(t_ref{1,4})/fs) + (height(t_ref{1,5})/fs) + (height(t_ref{1,6})/fs)+ (p/2);
                end
            end
        end
    elseif p == 30
        for i = 1:width(t_ref)
            xx = 31;
            % 1 s at the start to account for manual operation
            if i == 1
                for j = 1:32
                    if j == 1
                        x(j) = 1 + (p/2);
                    else
                        x(j) = x(j-1) + p + traverse.t_int_nom;
                    end
                end
            elseif i == 2
                for j = 1:xx
                    x(j+((i-1)*32)) = x((j)+((i-1)*xx)) + p + traverse.t_int_nom;
                end
            elseif i == 3
                for j = 2:xx
                    x(j+((i-1)*xx)) = x((j-1)+((i-1)*xx)) + p + traverse.t_int_nom;
                end
            end
        end
    end
    x_win_a = round(x-(p/2));
    x_win_b = round(x+(p/2));
    yl = ylim;
    for i = 1:length(x_win_a)
    end
    xline(x(:),':')
    for i = 1:length(x)
        data_win_pre_rot{i} = u_i((fs*x_win_a(i)) : (fs*x_win_b(i)),:);
    end
elseif ques_profile == 'R'
    win = traverse.t_vmax_nom; % period of interest (secs) for time averaging
    
    win_end = cellfun('length',t_ref); % end of each data window before a change in height
    win_end = [0,win_end];
    win_end_cum = cumsum(win_end); % cumulative
    % Establish centre of data window at each value of 'z' and define the window limits
    x = 1 + p/2 + traverse.t_int_nom*[0:length(coords_y)-1];
    x_win_a = round(x-(p/2));
    x_win_b = round(x+(p/2));
    % Data windowing using the centre of the window as a datum
    for i = 1:height(data_win_pre_rot)
        for j = 1:width(data_win_pre_rot)
            % Top to bottom = earliest to latest data acquisition
            data_win_pre_rot{i,j} = u_i( ((fs*x_win_a(j) : (fs*x_win_b(j))) + win_end_cum(i)) , :);
        end
    end
end

clear xx yl

%% PITCH AND YAW - PRE-ROTATION

% Preallocation
if ques_profile ~= 'R'
    p_bar_pre_rot = zeros(1,meas_points); % mean pitch for each data window
    y_bar_pre_rot = zeros(1,meas_points); % mean yaw for each data window
else
    p_bar_pre_rot = zeros(length(coords_z),length(coords_y));
    y_bar_pre_rot = zeros(length(coords_z),length(coords_y));
end

if ques_profile ~= 'R'
    for i = 1:length(data_win_pre_rot)
        data_win_py = data_win_pre_rot{i};
        [y_bar_pre_rot(i),p_bar_pre_rot(i)] = pitch_yaw(data_win_py);
    end
else
    for i = 1:height(data_win_pre_rot)
        for j = 1:width(data_win_pre_rot)
            data_win_py = data_win_pre_rot{i,j};
            [y_bar_pre_rot(i,j),p_bar_pre_rot(i,j)] = pitch_yaw(data_win_py);
        end
    end
end

% Select pitch and yaw corrections
if ques_profile == 'Z'
    idx = any(coords == 0);
    if idx == logical(1) 
        % Calibration point in the centre of the transect (Z = 0), provided the transect is positioned at Y = 0
        idx = find(coords == 0);
    elseif idx == logical(0)
        % Relax requirement that uses the centre point to calibrate the rotation
        if meas_points > 90
            % Identify the relevant centreline measurement file and load the data
            centreline_file = num2str(params.ID(1)-1);
            if length(centreline_file) < 4
                centreline_file = strcat('0',centreline_file);
            end
            base_path = '...add path here';
            filepath = strcat(base_path,'\',centreline_file,'.mat'); % create filepath
            load(filepath,'u_i_pre_rot')
        else
            % iii) use final measurement - useful for streamwise profiling where decay is significant
            idx = meas_points;
        end
    else
    end
    if meas_points < 90
        y_bar_pre_rot_ref = y_bar_pre_rot(idx);
        p_bar_pre_rot_ref = p_bar_pre_rot(idx);
    else
        % Derive pitch and yaw from reference dataset along the centreline
        [y_bar_pre_rot_ref,p_bar_pre_rot_ref] = pitch_yaw(u_i_pre_rot);
    end
elseif ques_profile == 'X'
    if any(params.Y1 & params.Z1) ~= 0
        error('You need a calibration point along the centreline.')
    else
    end
    % Each measurement shall be rotated of its own accord, i.e. centreline points are all calibration points - this is logical as the ADV may undergo slight rotation as the traverse is moved manually down the flume
    y_bar_pre_rot_ref = y_bar_pre_rot;
    p_bar_pre_rot_ref = p_bar_pre_rot;
elseif ques_profile == 'Y'
    if rem(meas_points, 2) == 1
        idx_1 = round(meas_points/2);
    elseif rem(meas_points, 2) == 0
        idx_1 = meas_points/2;
    else
    end
    idx_2 = idx_1 + 1;
    idx = [idx_1 idx_2]; 
    y_bar_pre_rot_ref = 0.5*(y_bar_pre_rot(idx_1) + y_bar_pre_rot(idx_2)); % in the absence of a centrepoint measurement, take the mean of the two points closest to it
    p_bar_pre_rot_ref = 0.5*(p_bar_pre_rot(idx_1) + p_bar_pre_rot(idx_2));
elseif ques_profile == 'R'
    idx_y = any(coords_y == 0);
    idx_z = any(coords_z == 0);
    if idx_y == 1 & idx_z == 1
        % Calibration point in the centre of the transect (Y, Z = 0)
        idx_y = find(coords_y == 0);
        idx_z = find(coords_z == 0);
        y_bar_pre_rot_ref = y_bar_pre_rot(idx_y,idx_z);
        p_bar_pre_rot_ref = p_bar_pre_rot(idx_y,idx_z);
    else
        error('Rethink your calibration point.')
    end
else
end

clear base_path centreline_file i opts path t_ref

%% MEAN VELOCITIES - PRE-ROTATION

% Mean velocity
if ques_profile ~= 'R'
    u_i_bar_pre_rot = cell2mat(cellfun(@mean, data_win_pre_rot(:), 'UniformOutput', false))'; % mean of each column in cell array
else
    u_i_bar_pre_rot = cell2mat(cellfun(@mean, data_win_pre_rot(:), 'UniformOutput', false))';
    u_i_bar_pre_rot = reshape(u_i_bar_pre_rot,[3,height(data_win_pre_rot),width(data_win_pre_rot)]); % 3 velocity components (rows) x no. Z points (cols) x no. Y points (dim3)
end

%% APPLY ROTATION

% Mean reference velocities from time series using the selected index
if meas_points > 90
    u_i1_bar_pre_rot = mean(u_i_pre_rot(:,[1:3]));
    clear u_i_pre_rot
elseif ques_profile == 'Y' | ques_profile == 'Z'
    u_i1_bar_pre_rot = (u_i_bar_pre_rot(:,idx))';
elseif ques_profile == 'R'
    u_i1_bar_pre_rot = (u_i_bar_pre_rot(:,idx_y,idx_z))';
elseif ques_profile == 'X'
    % Here every centreline measurement constitutes a calibration point
    u_i1_bar_pre_rot = u_i_bar_pre_rot;
end
                    
if ques_profile == 'R' | ques_profile == 'Y' | ques_profile == 'Z'
    % Reference 2D array data (pre-rotation) - for both vertical velocities (w1 and w2)
    u_i1_pre_rot = u_i;
    u_i2_pre_rot = u_i1_pre_rot;
    u_i2_pre_rot(:,3) = [];
    u_i2_pre_rot(:,3) = w2;

    % Define rotation matrix based on centrepoint measurement (Y = 0, Z = 0)
    A = rotation_solver(u_i1_bar_pre_rot,DoR,-p_bar_pre_rot_ref,y_bar_pre_rot_ref);
                    
    % Apply rotation - to both 2D velocity arrays, covering both w1 and w2
    u_i1 = u_i1_pre_rot * (A');
    u_i2 = u_i2_pre_rot * (A');
    u_i_pre_rot = [u_i1_pre_rot u_i2_pre_rot(:,3)]; % array of pre-rotated velocity components u v w1 w2
    u_i = [u_i1 u_i2(:,3)]; % now with rotated velocity components u v w1 w2
                    
elseif ques_profile == 'X'
    % Preallocation
    u_i1 = cell(1,meas_points);
    u_i2_pre_rot = data_win_pre_rot; % reference array for 2D array including w2
    u_i2 = cell(1,meas_points);
                    
    % Treat each coordinate individually
    for i = 1:meas_points
        % Define rotation matrix based on time-averaged velocities (covering w1)
        A = rotation_solver(u_i1_bar_pre_rot(:,i)',DoR,-p_bar_pre_rot_ref(i),y_bar_pre_rot_ref(i));

        % Reference data for w2 (2D array already prepared - data_win_pre_rot)
        u_i2_pre_rot{i}(:,3) = w2(x_win_a(i):x_win_b(i),:);

        % Apply rotation - to both 2D velocity arrays, covering both w1 and w2
        u_i1{i} = data_win_pre_rot{i} * (A');
        u_i2{i} = u_i2_pre_rot{i} * (A');
    end
                    
    % Generate continuous time series via concatentation
    % i) Pre-rotation
    u_i1_pre_rot = cat(1, data_win_pre_rot{:});
    u_i2_pre_rot = cat(1, u_i2_pre_rot{:});
    temp = u_i2_pre_rot; % dummy variable
    temp(:,1:2) = [];
    u_i_pre_rot = [cat(1, data_win_pre_rot{:}) temp];
                    
    % ii) Rotation applied
    u_i2_temp = cat(1, u_i2{:});
    u_i2_temp(:,1:2) = [];
    u_i = [cat(1, u_i1{:}) u_i2_temp]; % contains 4 columns, each representing a velocity component with the appropriate rotational correction applied
    u_i1 = cat(1, u_i1{:}); 
    u_i2 = cat(1, u_i2{:});
                    
    clear p_bar_pre_rot_ref temp u_i2_temp y_bar_pre_rot_ref
end

clear a b c A data_win_py u_i1_pre_rot u_i2_pre_rot u_i1 u_i2 u_i1_bar_pre_rot u_i_pre_rot w2

%% DATA WINDOWING (2) - POST-ROTATION

% Preallocation
if ques_profile ~= 'R'
    data_win_post_rot = cell(1,meas_points); % data from each window to be stored here
else
    data_win_post_rot = cell(length(coords_z),length(coords_y));
end

% Data windowing
if ques_profile == 'X'
    for i = 1:meas_points
        data_win_post_rot{i} = u_i(x_win_a(i):x_win_b(i),:);
    end
elseif ques_profile == 'Y' | ques_profile == 'Z'
    for i = 1:length(x)
        data_win_post_rot{i} = u_i((fs*x_win_a(i)) : (fs*x_win_b(i)),:);
    end
elseif ques_profile == 'R'
    for i = 1:height(data_win_post_rot)
        for j = 1:width(data_win_post_rot)
            data_win_post_rot{i,j} = u_i( ((fs*x_win_a(j) : (fs*x_win_b(j))) + win_end_cum(i)) , :); % 3 velocity components per cell - Z coordinates = rows; Y coordinates = cols
        end
    end
else
end

%% PITCH AND YAW - POST-ROTATION

% Preallocation
if ques_profile ~= 'R'
    p_bar_post_rot = zeros(1,meas_points); % mean pitch for each data window
    y_bar_post_rot = zeros(1,meas_points); % mean yaw for each data window
else
    p_bar_post_rot = zeros(length(coords_z),length(coords_y));
    y_bar_post_rot = zeros(length(coords_z),length(coords_y));
end

% Calculation
if ques_profile ~= 'R'
    for i = 1:length(data_win_post_rot)
        data_win_py = data_win_post_rot{i};
        [y_bar_post_rot(i),p_bar_post_rot(i)] = pitch_yaw(data_win_py);
    end
else
    for i = 1:height(data_win_post_rot)
        for j = 1:width(data_win_post_rot)
            data_win_py = data_win_post_rot{i,j};
            [y_bar_post_rot(i,j),p_bar_post_rot(i,j)] = pitch_yaw(data_win_py);
        end
    end
end

clear data_win_py

%% DEFINE REFERENCE VELOCITY

if ques_profile ~= 'R'

    % First calculate the mean modulus at each measurement point
    if meas_points > 90
        load(filepath,'U0_RP_primary','U0_mod')
    else
        U0_mod_bar = zeros(1,length(data_win_post_rot)); % preallocation
        for i = 1:length(data_win_post_rot)
            U0_mod_bar(i) = mean(sqrt(data_win_post_rot{i}(:,1).^2 + data_win_post_rot{i}(:,2).^2 + data_win_post_rot{i}(:,3).^2));
        end
    end

    % i) Mean modulus
    if meas_points > 90
        % Currently no centreline measurement to extract a modulus from discrete centreline measurement to be taken at some point
    elseif exist('idx','var') == 1
        if length(idx) == 2
            % Condition to extract the modulus from two coordinates as Y profiles do not cross the centre of the flume 
            U0_mod_centre_bar = mean([U0_mod_bar(idx_1) U0_mod_bar(idx_2)]);
        else
            % We go with the centreline velocity as reference
            U0_mod_centre_bar = U0_mod_bar(idx);
        end
    elseif ques_profile == 'X'
        U0_mod_centre_bar = U0_mod_bar; % enter all reference values into the 'find_U0' function... all values are centreline velocities
    else
    end
    
    % ii) Mean scalar
    % We need to extract two values here - one at the rotor plane and another at the plane of interest (if we are just considering that plane in isolation)
    if meas_points < 90
%         [U0_RP_secondary,U0_RP_primary] = find_U0(params,params.Mot_percent(1),U0_mod_centre_bar);
    else
        % Boundary layer measurements need a discrete measurement in the centreline to capture the modulus (see above)
    end

elseif ques_profile == 'R'

    % First calculate the mean modulus at each measurement point
    % TARGET STRUCTURE = rows: depth coordinate; cols: x-stream coordinate
    U0_mod_bar = zeros(width(data_win_post_rot),height(data_win_post_rot));
    for i = 1:height(data_win_post_rot)
        for j = 1:width(data_win_post_rot)
            U0_mod_bar(i,j) = mean(sqrt(data_win_post_rot{i,j}(:,1).^2 + data_win_post_rot{i,j}(:,2).^2 + data_win_post_rot{i,j}(:,3).^2));
        end
    end

    % i) i) Mean modulus along the centreline
    U0_mod_centre_bar = U0_mod_bar(idx_y,idx_z);

    % ii) Mean scalar
    % We need to extract two values here - one at the rotor plane, i.e. 5 m, and another at the plane of interest (if we are just considering that plane in isolation)
    [U0_RP_secondary,U0_RP_primary] = find_U0(params,params.Mot_percent(1),U0_mod_centre_bar);

else
end

clear idx_1 idx_2

%% NOISE FLOOR

if p >= 250
    
    % Create folders for saving spectra across each measurement point
    if ques_profile == 'Y'
        path = fullfile(cd,sprintf('%04.f',params.ID));
    elseif ques_profile == 'Z' | ques_profile == 'R' | ques_profile == 'X'
        path = fullfile(cd,sprintf('%04.f-%04.f',params.ID(1),params.ID(end)));
    else
    end
    mkdir(path)
    
    % Spectra - non-scalar structure
    Spectra = struct([]);
    
    if ques_profile ~= 'R'
        for i = 1:length(coords)
            pos = coords(i);
            [Spectra(i).pos] = noise_floor_analysis(data_win_post_rot{i},U0_mod_bar(i),fs,params,path,pos);
        end
    elseif ques_profile == 'R'
        % Structure should contain Y points in the rows and Z points along the columns
        seq = "z_"+coords_z; % Title of each column
        coords_y = fliplr(coords_y);
        for i = 1:height(data_win_post_rot)
            path_dyn = [path filesep seq{i}];
            mkdir(path_dyn)
            for j = 1:width(data_win_post_rot)
                name = convertCharsToStrings(seq{i});
                name = strrep(name,'.','p'); % p = decimal (p)oint
                name = strrep(name,'-','m'); % m = (m)inus
                pos = coords_y(j);
                % Spectra structure = Y pos (rows) x Z pos (cols)
                [Spectra(j).(name)] = noise_floor_analysis(data_win_post_rot{i,j},U0_mod_bar(i,j),fs,params,path_dyn,pos);
            end
        end
        % Define coordinate system for data analysis... only applies if the no. of spatial increments in both directions are equal
        coords(:,:,1) = coords_y; 
        coords(:,:,2) = coords_z;
    end

else
    fprintf('\nSampling time insufficient to ascertain noise floor via spectral method.\n\n')
end

%% PASS VARIABLES OF INTEREST

% Clear variables
clear coords_y coords_z data_win_pre_rot DoR fs h i idx idx_y idx_z inc j leg list name params path path_dyn pos ques_profile seq t u_i win win_end win_end_cum x x_win_a x_win_b yl

% Save workspace variables into a structure
list = who;
for iVar = 1:length(list)
  Post_process.(list{iVar}) = eval(list{iVar});
end

clearvars -except coords Post_process
close all
end

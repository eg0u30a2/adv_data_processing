
function  [Spectra] = noise_floor_analysis(u_i_data,U0_mod,fs,params,fig_path,variable,ques_anal,quiet)

% PURPOSE: this function calculates the power spectral density (PSD) via FFT, with and without averaging 
% (moving average with equally spaced frequency interval in log-space) as a function of frequency.
% Noise floor calculated as the variance included in the tail of the spectrum. Total noise estimated via
% integration over the whole range of frequencies, assuming white noise.

% Arguments IN
% u_i_data = 2D velocity array (m/s)
% fs = sampling freq. (Hz)
% params = test parameters
% fig_path = path where figures will be saved to
% U0_mod = reference modulus velocity (m/s)
% variable = variable of interest, e.g. speed point, measurement coordinate, etc. (for saving filenames)
% ques_anal = nature of analysis, i.e. a profile or a single point measurement
% quiet = logical value for whether the data has been generated from a still tank

% Arguments OUT
% Spectra = power spectral density of velocity fluctuations (m^2/s)

% Version control:
% 31/05/2022 - Method A for noise floor identification chnages from mean of final half decade to 95-100 Hz. Method B changes from min. gradient (polynomial fit) to the final value, i.e. @ Nyquist freq.
% 06/06/2022 - filter order as input to low-pass filter
% 08/06/2022 - f_min changed to 99 Hz
% 11/06/2022 - sigma revised to sigma_fluc, i.e. std dev of the velocity signal is now RMS of the velocity fluctuations. f_min reverted to 95 Hz to improve robustness against sudden deviations at the tail.
% 19/06/2022 - f_cutoff changed from 100 Hz (f_w_cut(end)) to 99.9 Hz in order to low-pass filter at circa Nyquist freq. rather than half-power freq. f50
% 08/10/2022 - noise floor method B inferred from the raw power spectra rather than the Welch spectra
% 12/11/2022 - noise induced bias corrections referred to integrated area under PSD rather than using the approximation cited in Richard, et al. (2013)
% 13/01/2023 - U0 becomes an input argument to the function, rather than calculated inside the function
% 08/03/2023 - yticks for PSD plot with noise floor are calculated parametrically

%% SETUP

% Traverse in a still tank?
if nargin == 6
    ques_anal = 2;
    quiet = NaN;
elseif nargin == 7
    quiet = NaN;
else
end

% Figures
if ques_anal == 2
    coord = [];
    if variable > 50
        suffix = 'rpm';
    else
        suffix = 'm';
    end
    fprintf('\nIncrement: %.3f %s \n----------------------',variable,suffix)
    % set(0,'DefaultFigureVisible','off')
elseif ques_anal == 1
    clear path
    fig_path = '...add path here';
    coord = params.X + "-" + params.Y1 + "-" + params.Z1; % string for filename
    set(0,'DefaultFigureVisible','on')
end

%% DETRENDING

data_ref = u_i_data; % duplicate the input array for later use

% Detrend the data as the mean value does not play a role in PSD spectrum
u_i_data = detrend(u_i_data(:,:));
if mod(height(u_i_data),2) ~= 0
    u_i_data = u_i_data(1:end-1,:);
end    

%% PSD CALCULATION - FOURIER TRANSFORM

% Calculate spectra via FFT
[len,u_i_psd,f_psd] = psd_fft(u_i_data,fs);

%% IDENTIFY NOISE FLOOR
psd_cut_min = i(end,:);

% METHOD B
% Mean PSD over a freq. band
f_min = 95; % select the lowest freq. (Hz) in the band of the tail (upper = Nyquist fs/2)
u_i_psd_mean_final_portion = mean((u_i_psd((f_psd>f_min),:))); % portion over the final half-freq decade band

%% DETRENDING

% Detrend the data (leaving pure fluctuations) as the mean value does not play a role in PSD spectrum
u_i_fluc = detrend(data_ref(:,:));

%% CALCULATE VARIANCE DUE TO NOISE

% Noise contribution (Methods A and B)
noise_var_B = ones(length(f_psd),width(data_ref)).*u_i_psd_mean_final_portion;

% Integrate areas, which equate to a variance
area_psd = trapz(f_psd,u_i_psd);
area_noise_B = trapz(f_psd,noise_var_B);

%% PASS VARIABLES OF INTEREST

% Clear variables
clear area_noise_A area_psd b coeff d1y vec u_i_fluc data_ref fs f_min f_threshold f_cut fig_path h i I_3D_delta I_i_delta
clear interval leg len low mu n noise_A_percent_area noise_B_percent_area noise_std_A noise_var_A noise_var_B params path plot_len psd_cut psd_cut_min ques_anal quiet
clear sigma_i s S vel_var u_i_PWelch_cut u_i_data upp xl yl X X1 Xf y Y

% Save workspace variables into a structure
list = who;
for iVar = 1:length(list)
  Spectra.(list{iVar}) = eval(list{iVar});
end

clearvars -except Spectra
close all
end


function [t,u_i_raw,SNR,corr,amp] = variables(data)

% PURPOSE - Assign variables and headers
%
% Arguments IN:
% data = 2D velocity array (m/s)
%
% Arguments OUT:
% t = 2D time array (secs)
% u_i_raw = 2D velocity array (m/s)
% SNR = 2D SNR array
% corr = 2D corr array
% amp = 2D amplitude array

% Assign headers
data.Properties.VariableNames = {'X1','Time','Step','X2','Vx','Vy','Vz1','Vz2','Amp1','Amp2','Amp3','Amp4','SNR1','SNR2','SNR3','SNR4','Corr1','Corr2','Corr3','Corr4'};

% Assign time t (secs)
t = data{:,2};

% Assign velocity V
u_i_raw = data{:,5:8};

% Assign SNR
SNR = data{:,13:16};

% Assign correlation
corr = data{:,17:20};

% Assign amplitude
amp = data{:,9:12};

end

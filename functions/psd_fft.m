
function  [len,psd,f] = psd_fft(data,fs)

% This function calculates the discrete Fourier transform (DFT) using FFT.
% Spectrum is sampled at a no. of points equal to a power of 2, allowing
% the matrix algebra to be sped up. The FFT samples the signal energy at discrete frequencies.
%
% Arguments IN
% data = 2D velocity array
% fs = sampling freq. (Hz)
%
% Arguments OUT
% f = freq. vector (Hz)
% len = length of dataset
% psd = PSD estimate using FFT
%
% Version control
% 05/06/2022 - Initial version

% Calculate spectra via FFT
len = length(data);
psd = (abs(fft(data(:,:),len)).^2)/(len*fs); % Fourier transform (2nd argument = vector length, such that the DFT returns n-points). Amplitude squared for PSD estimate, as this equals power
psd = 2.*psd(2:len/2+1,:); % single-sided spectrum. (1:N/2+1) takes only one half of the spectrum; the other half is its reflection
f = fs/2*linspace(0,1,len/2+1).'; % Nyquist frequency of the signal = fs/2
f = f(2:end); % remove zero Hz component

end

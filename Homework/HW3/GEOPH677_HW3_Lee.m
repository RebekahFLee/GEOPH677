%% Homework 3: F-k analysis
%% Part 1
%%
% * 1. Load station coordinates and plot in map view*
clear all
close all
clc

% load stations
load('stationCoordinatesYX_in_km.txt')
Y = stationCoordinatesYX_in_km(:,1);
X = stationCoordinatesYX_in_km(:,2);

% plot stations
scatter(X,Y,'filled');
xlabel('X coordinates [km]')
ylabel('Y coordinates [km]')

% figure properties
fsize=16;
figproperties
grid on;
%%
% * 2. Write a function that computes the array response for a given station
% array *

% function [ARF,s,theta,kx,ky] =arrayResp(c_min,freq,phase,X,Y)
% calculates the array response given a minimum velocity, frequency, phase from data or test source, and
% station coordinates. Plots in both cartesian and polar coordinates. Created for GEOPH 677 HW3 assignment 3/2/17
% ** requires polarPcolor function ** 
% Rebekah Lee
 % USAGE ARF =arrayResp(c_min,freq,phase,X,Y)
 %INPUTS
 % c_min: minimum velocity to search over
 % frequ: freq to test
 % phase: phase of source or data: (k0 dot r)
 % X: x coordinates of stations
 % Y: y coordinates of stations
 % plotflag 1 to plot 
 %---------------------------------------------------------------------
%  %---------------------------------------------------------------------
%  
% % create slowness vector (hardcoded numbers are arbitrary, they seem to work
% % well)
% c_max = c_min*50;
% nc = 2000;
% c = linspace(c_max,c_min,nc);
% s=1./c;
% 
% % create grid of thetas
% ntheta = 360;
% theta = linspace(0,360,ntheta);
% Theta = repmat(theta',1,nc);
% 
% % create grid of ray parameter
% P = repmat(s,ntheta,1);
% 
% % calculate K components
% w = 2*pi*freq;
% kx = w.*P.*cosd(Theta); 
% ky = w.*P.*sind(Theta);
% %plot(kx,ky,'b*')
% 
% % loop and sum through stations
% n = length(X);
% temp = zeros(ntheta,nc);
% for ista = 1:n
% k_dot_r = kx*X(ista) + ky*Y(ista);
% k0_dot_r = phase(ista);
% temp = temp + exp(1i*k_dot_r).*exp(-1i*k0_dot_r);
% end
% 
% % Array Response Function
% ARF = abs((1/n*temp));

%end

%% 
% * 3. Compute array response function Using f= 1 Hz, theta = 45 degrees, cmin = 1km/s and capp = 3km/s

% set source variables
f0= 1; %frequency Hz
w0 = 2*pi*f0; %angular frequency
v0 = 3; %V apparent m/s

% Find K components for source
k0x = w0/v0*cosd(45); 
k0y = w0/v0*sind(45);

%inputs
c_min = 1;
freq = 1; %[hz]
phase = k0x * X + k0y * Y;

%get array response and plot
[ARF,s,theta,kx,ky] = arrayResp(c_min,freq,phase,X,Y);
%Db = 20*log10(ARF);

% plot array response
fsize=16;
% plot in Cartesian and Polar Coordinates
%Cartesian plot
figure
subplot(1,2,1)
pcolor(kx,ky,ARF); 
shading flat
cl = colorbar;
ax = caxis;
ylabel(cl,'Normalized Power')
xlabel('K_x')
ylabel('K_y')
title('Array Response: \theta = 45^{\circ}, V_{app} = 3 km/s ')
figproperties

% Polar plot
subplot(1,2,2)
polarPcolor(s,theta,ARF,'Nspokes',9,'Ncircles',4)
caxis(ax);
text(.3533,.05,'Slowness')
title('Array Response: Polar Coordinates')
figproperties
set(gcf,'Position',[52 353 1383 469])

%% Part 2
%%
% * 1. Load the data
load('array_data_example.mat')
raw = dta;
%%
% Define sampling and nyquist frequencies
fs= 200; % sample frequency
dt = 1/fs; % time sample interval
fNyq = fs/2; % Nyquist sampling frequency
[n,npts]=size(dta);
t = 0:npts-1;
t = t*dt;
%%
% preprocessing: demean and filter data. I just played around until I thought the data looked good. 
temp = repmat(mean(dta,2),1,npts);
dta = dta - temp;

% create butterworth filter
fmin=.5;
fmax=99;
btord = 2;
[z,q]=butter(btord,[fmin/fNyq fmax/fNyq],'bandpass');

% filter the data
dtaFilt=filtfilt(z,q,dta');
raw = raw';
dta = dta';
%plot(raw(1,:))
%hold on;
%plot(test(1,:))
% plot(dtaFilt(:,1))
% %legend('raw','demeaned','filtered')
% axis tight
% hold off
%%
% * 2. Plot time domain and amplitude spectral densities for each station

% find nearest power of two to the lenght of the data to speed up fft
nfft = 2^nextpow2(npts);

% create frequency vector (for positive part)
df = fs/nfft;
fArray = (0:(nfft-1)/2 +1).*df;

fDomain = fft(dtaFilt,nfft);
fDomain = fDomain(1:nfft/2+1,:); % just get the positive part of the spectrum

figure;
for ii= 1:n % loop over stations
    %plot time series
    subplot(8,2,ii*2-1)
    plot(t,dtaFilt(:,ii))
    axis tight
    %plot amplitude spectrum
    subplot(8,2,ii*2)
    plot(fArray,abs(fDomain(:,ii)))
    axis tight
    xlim([.5 4])
end

% set some general properties to the figure
set( findall( gcf, '-property', 'LineWidth' ), 'LineWidth', 2 );
set(gcf,'units','inches')
set(gcf,'position',[1.5 .98 16.2 10])

% Add labels
fsize=24;

% xlabels
xlabel('Frequency [Hz]')
set(get(gca,'XLabel'),'fontsize',fsize);

subplot(8,2,15)
xlabel('Time (s)')
set(get(gca,'XLabel'),'fontsize',fsize);

% Ylabels
subplot(8,2,10)
ylabel('        Amplitude')
set(get(gca,'YLabel'),'fontsize',fsize);

subplot(8,2,9)
ylabel('        Velocity (m/s)')
set(get(gca,'YLabel'),'fontsize',fsize);


% add titles
subplot(8,2,1)
title('Filtered Time Series, Stations 1-8','fontsize',fsize)

subplot(8,2,2)
title('Amplitude Spectral Densities, stations 1-8','fontsize',fsize)


%% 
% * 3. Beamform at 1Hz
Phase = angle(fDomain); 
[~,idx] = min(abs(abs(fDomain)-freq));
for ista = 1:n
phase(ista) = angle(fDomain(idx(ista),ista));
end
[beam,s,theta,kx,ky] = arrayResp(c_min,freq,phase,X,Y);
%%
% plot
figure
[~,c] = polarPcolor(s,theta,beam,'Nspokes',18,'Ncircles',4);
caxis(ax);
text(.3533,.05,'Slowness')
title(['Beamforming: f= ',num2str(freq),'Hz'])
ylabel(c,'Power')
figproperties
%set(gcf,'Position',[52 353 1383 469])
%% 
% * 4. Choose two other frequencies to beamform and plot
freq2 = 1.5; % [Hz]
freq3 = 2;
beam2 = arrayResp(c_min,freq2,phase,X,Y);
beam3 = arrayResp(c_min,freq3,phase,X,Y);

% plot
%Cartesian plot
figure
subplot(1,2,1)
[~,c] = polarPcolor(s,theta,beam2,'Nspokes',18,'Ncircles',4);
caxis(ax); % keep original colorbar axis
text(.3533,.05,'Slowness')
title(['Beamforming: f= ',num2str(freq2),'Hz'])
ylabel(c,'Power')
figproperties
set(gcf,'Position',[52 353 1383 469])


subplot(1,2,2)
[~,c] = polarPcolor(s,theta,beam3,'Nspokes',18,'Ncircles',4);
caxis(ax);
text(.3533,.05,'Slowness')
title(['Beamforming: f= ',num2str(freq3),'Hz'])
ylabel(c,'Power')
figproperties
set(gcf,'Position',[52 353 1383 469])
%% 
% * 5 Discussion
%%
% Generally speaking the three plots show energy between 275 and 296
% degrees as well as between about 84 and 106 degrees. 
% At 1.5 and 2 Hz there is a signal at around 339, 233, 180 and ~140 degrees azimuth that are not apparent (or barely visible) at
% 1Hz. The higher frequencies both plot energy at around these azimuths but
% with different slownesses. This could reflect different phase arrivals
% for the same source. However, some of these signals are artifacts. Perhaps the lower hemisphere signals are the artifacts, since this pattern would be similar to the array response pattern in part 1. 
%%
% All plots show something between 275 and 318 degree
% azimuth at a small slowness (higher velocity of at least 3 km/s).
%%
% I detrended and then filtered the data to get rid of the high amplitude
% in low frequencies that I saw in the fft. Removing this significantly
% altered the beamforming (not plotted). I removed the low frequency in my
% bandpass filter. If this data is for a larger amplitude earthquake that
% may be a mistake. 
% I went up to near the nyquist frequency for my filter.
% Most of the energy is under about 3.5 or 4 Hz. I could try filtering out
% the higher frequencies.
%% Extra Credit
% * F-k at multiple frequencies
start = 1.5; % beginning frequency [Hz]
stop = 3; %last frequency [Hz]

frequencies = linspace(start,stop,20);
nf = length(frequencies);
multiF = zeros(size(beam));
for ii = 1: nf
   temp = arrayResp(c_min,frequencies(ii),phase,X,Y);
   multiF = multiF +temp;
end

figure;
[~,c] = polarPcolor(s,theta,multiF./length(frequencies),'Nspokes',18,'Ncircles',4);
caxis(ax);
text(.3533,.05,'Slowness')
title(['Beamforming: f= ',num2str(frequencies(1)),' - ',num2str(frequencies(end)),' Hz'])
ylabel(c,'Normalized Power')
figproperties

%%
% I used the amplitude spectrum to choose the frequency band. There were a
% lot of peaks in the amplitude between about 1.75 and 2.25 Hz so I used
% these as the limits of frequencies I tested. There is some smearing of
% the image as it sums over the different frequencies/phases. This could be
% because of the small number of stations
%%
% areas I noted above are visible at about 296, 339 and 140 degrees
% azimuth. Either these respresent separate sources with different phases
% or they are artifacts of the beamforming. 
%%
% 
%%

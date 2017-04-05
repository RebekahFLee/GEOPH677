function [ARF,s,theta,kx,ky] =arrayResp(c_min,freq,phase,X,Y)
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
 %---------------------------------------------------------------------
 
% create slowness vector (hardcoded numbers are arbitrary, they seem to work
% well)
c_max = c_min*50;
nc = 2000;
c = linspace(c_max,c_min,nc);
s=1./c;

% create grid of thetas
ntheta = 360;
theta = linspace(0,360,ntheta);
Theta = repmat(theta',1,nc);

% create grid of ray parameter
P = repmat(s,ntheta,1);

% calculate K components
w = 2*pi*freq;
kx = w.*P.*cosd(Theta); 
ky = w.*P.*sind(Theta);
%plot(kx,ky,'b*')

% loop and sum through stations
n = length(X);
temp = zeros(ntheta,nc);
for ista = 1:n
k_dot_r = kx*X(ista) + ky*Y(ista);
k0_dot_r = phase(ista);
temp = temp + exp(1i*k_dot_r).*exp(-1i*k0_dot_r);
end

% Array Response Function
ARF = abs((1/n*temp));


% fsize=16;
% % plot in Cartesian and Polar Coordinates
% %Cartesian plot
% figure
% subplot(1,2,1)
% pcolor(kx,ky,ARF); 
% shading flat
% cl = colorbar;
% %ylabel(cl,'Label')
% xlabel('K_x')
% ylabel('K_y')
% title('Array Response: \theta = 45^{\circ}, V_{app} = 3 km/s ')
% figproperties
% 
% % Polar plot
% subplot(1,2,2)
% polarPcolor(s,theta,ARF,'Nspokes',9,'Ncircles',4)
% text(.3533,.05,'Slowness')
% title('Array Response: Polar Coordinates')
% figproperties
% set(gcf,'Position',[52 353 1383 469])

end

function seismicMoment(M,magoff)
% function to plot moment tensor in 2D and 3D. 
% USAGE: seismicMoment(M, magoff)
%
% INPUTS:
%   M 3rd order tensor (3x3 matrix) for seismic moment
%   magoff - set to 1 to just plot the sign, otherwise it will keep the
%   magnitude
% Rebekah Lee 
% GEOPH 677 Earthquake Seismology, 
% April 5 2017

dphi = 1; %step in azimuth in degrees
di = 1; % step in take off angle in degrees

% Get vectors for azimuth and take-off angles
AZ = .01:dphi:360.01; 
n=length(AZ);
take_off = 0.01:di:90.01;
m = length(take_off);

%preallocate ------------------------------------
lx = zeros(n,m);
ly = zeros(n,m);
lz = zeros(n,m);

px = zeros(n,m);
py = zeros(n,m);
pz = zeros(n,m);

phi_x = zeros(n,m);
phi_y = zeros(n,m);
phi_z = zeros(n,m);

U_p = zeros(n,m);
U_SV = zeros(n,m);
U_SH = zeros(n,m);
%------------------------------------
% Find unit vectors for every azimuth, take-off angle pair

for iAZ = 1:n % loop over azimuth
   for iray = 1:m % loop over take off angles
      % l hat
      lx(iAZ,iray) = sind(take_off(iray))*cosd(AZ(iAZ));
      ly(iAZ,iray) = sind(take_off(iray))*sind(AZ(iAZ));
      lz(iAZ,iray) = cosd(take_off(iray));
      lhat = [lx(iAZ,iray), ly(iAZ,iray), lz(iAZ,iray)];
      
      %p hat
      px(iAZ,iray) = cosd(take_off(iray))*cosd(AZ(iAZ));
      py(iAZ,iray) = cosd(take_off(iray))*sind(AZ(iAZ));
      pz(iAZ,iray) = -sind(take_off(iray));
      phat = [px(iAZ,iray), py(iAZ,iray), pz(iAZ,iray)];
      
      %phi hat
      phi_x(iAZ,iray) = -sind(AZ(iAZ));
      phi_y(iAZ,iray) = cosd(AZ(iAZ));
      phi_hat = [phi_x(iAZ,iray), phi_y(iAZ,iray), phi_z(iAZ,iray)];
      
      % find first motions
      tmp = [dot(lhat,M(1,:)), dot(lhat,M(2,:)), dot(lhat,M(3,:))];
      U_p(iAZ,iray) = dot(tmp,lhat'); % P wave
      
      tmp = [dot(phat,M(1,:)), dot(phat,M(2,:)), dot(phat,M(3,:))];
      U_SV(iAZ,iray) = dot(tmp,lhat'); % SV wave
      
      tmp = [dot(phi_hat,M(1,:)), dot(phi_hat,M(2,:)), dot(phi_hat,M(3,:))];
      U_SH(iAZ,iray) = dot(tmp,lhat'); % SH wave
   end
end

% create grids of spherical coordinates
dist = 1; % create distances (unitless)
elev = deg2rad(90 - take_off); % get elevation angle for sph2cart
azimuth = deg2rad(360-AZ); % convert azimuths to radians and azimuth from X axis counterclockwise
[elev, azimuth, distance] = meshgrid(elev,azimuth,dist); 

% convert coordinates to cartesian
[y, x, z] = sph2cart(azimuth,elev,distance);
z = -z;
if magoff
  U_p = sign(U_p); 
  U_SV = sign(U_SV);
  U_SH = sign(U_SH);
end

% Plot P wave radiation pattern ----------------------------------------
figure 
subplot(3,2,1) % 3D plot
surf(x,y,z,U_p,'EdgeColor','none')

colorbar;
ax = gca;
ax.Color = 'blue';
%view([-38 45])

fsize = 16;
xlabel('x - North','Fontsize',fsize)
ylabel('y - East','Fontsize',fsize)
zlabel('z - Up','Fontsize',fsize)
title('Radiation pattern P','Fontsize',fsize)
set( findall( gcf, '-property', 'TickLabelInterpreter' ), 'TickLabelInterpreter', 'Latex' );

subplot(3,2,2) %2D plot
surf(x,y,z,U_p,'EdgeColor','none')

view([-90 90])
colorbar;
axis ij
ax = gca;
ax.Color = 'blue';

fsize = 16;
xlabel('x - North','Fontsize',fsize)
ylabel('y - East','Fontsize',fsize)
title('Radiation pattern P','Fontsize',fsize)
set( findall( gcf, '-property', 'TickLabelInterpreter' ), 'TickLabelInterpreter', 'Latex' );

% Plot SV wave radiation pattern ----------------------------------------

subplot(3,2,3) % 3D plot
surf(x,y,z,U_SV,'EdgeColor','none')

colorbar;
ax = gca;
ax.Color = 'blue';
% view([-38 45])

fsize = 16;
xlabel('x - North','Fontsize',fsize)
ylabel('y - East','Fontsize',fsize)
zlabel('z - Up','Fontsize',fsize)
title('Radiation pattern SV','Fontsize',fsize)
set( findall( gcf, '-property', 'TickLabelInterpreter' ), 'TickLabelInterpreter', 'Latex' );

subplot(3,2,4) %2D plot
surf(x,y,z,U_SV,'EdgeColor','none')

view([-90 90])
colorbar;
axis ij
ax = gca;
ax.Color = 'blue';

fsize = 16;
xlabel('x - North','Fontsize',fsize)
ylabel('y - East','Fontsize',fsize)
title('Radiation pattern SV','Fontsize',fsize)
set( findall( gcf, '-property', 'TickLabelInterpreter' ), 'TickLabelInterpreter', 'Latex' );


% Plot SH wave radiation pattern ----------------------------------------
subplot(3,2,5) % 3D plot
surf(x,y,z,U_SH,'EdgeColor','none')

colorbar;
caxis([-1 1])
ax = gca;
ax.Color = 'blue';
% view([-38 45])

fsize = 16;
xlabel('x - North','Fontsize',fsize)
ylabel('y - East','Fontsize',fsize)
zlabel('z - Up','Fontsize',fsize)
title('Radiation pattern SH','Fontsize',fsize)
set( findall( gcf, '-property', 'TickLabelInterpreter' ), 'TickLabelInterpreter', 'Latex' );

subplot(3,2,6) %2D plot
surf(x,y,z,U_SH,'EdgeColor','none')

view([-90 90])
colorbar;
caxis([-1 1])
axis ij
ax = gca;
ax.Color = 'blue';

fsize = 16;
xlabel('x - North','Fontsize',fsize)
ylabel('y - East','Fontsize',fsize)
title('Radiation pattern SH','Fontsize',fsize)
set( findall( gcf, '-property', 'TickLabelInterpreter' ), 'TickLabelInterpreter', 'Latex' );

set(gcf,'Position',[192 56 736 749])

if magoff ==1
colormap('gray');
else
colormap('parula')
end
end
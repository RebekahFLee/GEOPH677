%% Part 1: Moment Tensors and Radiation Patterns
%%
dphi = 10; %step in azimuth in degrees
di = 2; % step in take off angle in degrees

% Moment Tensor
M = [1, 0, 0; 0, 0, 0; 0, 0, 0];

% Get vectors for azimuth and take-off angles
AZ = 0:dphi:360-dphi; 
n=length(AZ);
take_off = 0:di:90;
m = length(take_off);

%preallocate
lx = zeros(n,m);
ly = zeros(n,m);
lz = zeros(n,m);

px = zeros(n,m);
py = zeros(n,m);
pz = zeros(n,m);

phi_x = zeros(n,m);
phi_y = zeros(n,m);
phi_z = zeros(n,m);

% Find unit vectors for every azimuth, take-off angle pair
for iAZ = 10%:n % loop over azimuth
   for iray = 10%:m % loop over take off angles
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
      
      Up(iAZ,iray) = lhat * M *;
   end
end
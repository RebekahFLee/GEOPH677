function M = momentTensor(strike,dip,rake)
% momentTensor(strike,dip,rake) 

% convert to radians
strike = deg2rad(strike);
dip = deg2rad(dip);
rake = deg2rad(rake);

M0 = 1;

% Find moment tensor (Box 4.4 Aki and Richards)
M = zeros(3,3);
M(1,1) = -M0*(sin(dip)*cos(rake)*sin(2*strike) + sin(2*dip)*sin(rake)*sin(strike)^2);
M(1,2) = M0*(sin(dip)*cos(rake)*cos(2*strike) + 0.5* sin(2*dip)*sin(rake)*sin(2*strike));
M(1,3) = -M0*(cos(dip)*cos(rake)*cos(strike) + cos(2*dip)*sin(rake)*sin(strike));

M(2,1) = M(1,2);
M(2,2) = M0*(sin(dip)*cos(rake)*sin(2*strike) - sin(2*dip)*sin(rake)*cos(strike)^2);
M(2,3) = -M0*(cos(dip)*cos(rake)*sin(strike) - cos(2*dip)*sin(rake)*cos(strike));

M(3,1) = M(1,3);
M(3,2) = M(2,3);
M(3,3) = M0*sin(2*dip)*sin(rake);
end


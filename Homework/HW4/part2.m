%%
% Part 2, Number 1
%
% #1 and #2 nodal plane azimuths are both 0 and 90 degrees. There is zero displacement on the nodal planes.
% This happens along the arc lines the white and black.
%
% # 3 looks like the line is at about 35 degrees from North and on the
% second arc at about 125 degrees
%
% # 4  Nodal plane at about 75 and 170 degrees. I think this one is
% slightly dipping and the 2D version doesn't look like it's quite at 165.
%
% # 5 We are given an Azimuth of 358 which is what looks like the first arc with a thin shaded area.
% The other plane the other nodal plane is close to 90 degrees. 
%
% # 6 We are given an Azimuth of 294, the other nodal plane is at about 135
% degrees. Not exactly 90 degrees apart since the fault is dipping at 61
% degrees.
%%

M = [-0.5587, -0.7589, -0.1949; -0.7589, 0.6203, -0.1833; -0.1949, -0.1833, -0.0618]; 
m0 = 1/3*trace(M).*eye(3);
Mprime = M - m0;
seismicMoment(Mprime,magoff);

%%
% Part 2, Number 2

magoff = 0;
% 1. Right- lateral strike slip fault
M = [0, 1, 0; 1, 0, 0; 0, 0, 0];
seismicMoment(M,magoff);

%%
% 5.
strike = 358;
dip = 80;
rake = 59;

M = momentTensor(strike,dip,rake);
seismicMoment(M,magoff)
%%
% 6.
strike = 294;
dip = 61;
rake = 152;

M = momentTensor(strike,dip,rake);
seismicMoment(M,magoff)

%%
% For the P wave: The absolute value of the amplitude increases radially in
% all directions (except the nodal planes) for a perfictly vertical fault
% (strike slip faults in 1 and 2). For lesser dips the amplitude reaches a
% maximum and starts to fall again radially before reaching the outside of
% the beachball. For the strike slip there are two lows and two highs on
% the very rim of the beachball. For #6 and #5 (order of increasing
% azimuth) the highs and lows start to shift so that there is only one for
% each close to 350 degrees (# 5). In other words the amplitudes reach an
% absolute maximum in only two of the four map-vie quadrients.
% 
% 
% close all
%% Part 1: Moment Tensors and Radiation Patterns
%%

magoff = 1; % 1 to plot the sign
%% 
% 1. Right- lateral strike slip fault
M = [0, 0, 0; 1, 0, 0; 0, 0, 0];
seismicMoment(M,magoff)
%%
% 2. Left-lateral strike slip fault
M = [0, -1, 0; -1, 0, 0; 0, 0, 0];
seismicMoment(M,magoff)
%%
% 3.
M = [.4581, -0.2460, 0.4192;-0.2460, -.9332, -.2560; 0.4192 -0.256 0.4751];
seismicMoment(M,magoff)
%%
% 4.
M = [-0.5587, -0.7589, -0.1949; -0.7589, 0.6203, -0.1833; -0.1949, -0.1833, -0.0618]; 
seismicMoment(M,magoff)
%%
% 5.
strike = 358;
dip = 80;
rake = 59;

M = momentTensor(strike,dip,rake)
seismicMoment(M,magoff)
%%
% 6.
strike = 294;
dip = 61;
rake = 152;

M = momentTensor(strike,dip,rake)
seismicMoment(M,magoff)
%% online test
M = [0,0,0;0,0,-1;0,-1,0];
seismicMoment(M,magoff)
%%
M = [1,0,0;0,-2,0;0,0,1]*1/sqrt(6);
seismicMoment(M,magoff)
%%
M = [-2,0,0;0,1,0;0,0,1]*1/sqrt(6);
seismicMoment(M,magoff)
%%
M = [1,0,0;0,1,0;0,0,-2]*-1/sqrt(6);
seismicMoment(M,magoff)
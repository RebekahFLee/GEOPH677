function polarview(theta,rho,z,prn,maxii)

%POLARVIEW  Polar image plot.
%   POLARVIEW(THETA, RHO, Z, PRN) makes a geographically oriented plot using 
%   polar coordinates of the angle THETA, in radians, versus the radius RHO.
%
%   Clement.Ogaja@gmail.com
if isstr(theta) | isstr(rho)
    error('Input arguments must be numeric.');
end
if ~isequal(size(theta),size(rho))
    error('THETA and RHO must be the same size.');
end

% get hold state
cax = newplot;
next = lower(get(cax,'NextPlot'));
hold_state = ishold;

% get x-axis text color so grid is in same color
tc = get(cax,'xcolor');
ls = get(cax,'gridlinestyle');

% Hold on to current Text defaults, reset them to the
% Axes' font attributes so tick marks use them.
fAngle  = get(cax, 'DefaultTextFontAngle');
fName   = 'Helvetica';%get(cax, 'DefaultTextFontName');
fSize   = get(cax, 'DefaultTextFontSize');
fSize=16;
fWeight = 'bold';get(cax, 'DefaultTextFontWeight');
fUnits  = get(cax, 'DefaultTextUnits');
% set(cax, 'DefaultTextFontAngle',  get(cax, 'FontAngle'), ...
%     'DefaultTextFontName',   get(cax, 'FontName'), ...
%     'DefaultTextFontSize',   get(cax, 'FontSize'), ...
%     'DefaultTextFontWeight', get(cax, 'FontWeight'), ...
%     'DefaultTextUnits','data')
set(cax, 'DefaultTextFontAngle', fAngle , ...
    'DefaultTextFontName',   fName , ...
    'DefaultTextFontSize',   fSize, ...
    'DefaultTextFontWeight', fWeight, ...
    'DefaultTextUnits',fUnits );


% only do grids if hold is off
if ~hold_state

% make a radial grid
    hold on;
    maxrho = max(abs(rho(:)));
    hhh=plot([-maxrho -maxrho maxrho maxrho],[-maxrho maxrho maxrho -maxrho]);
    set(gca,'dataaspectratio',[1 1 1],'plotboxaspectratiomode','auto')
    v = [get(cax,'xlim') get(cax,'ylim')];
    ticks = sum(get(cax,'ytick')>=0);
    delete(hhh);
% check radial limits and ticks
    rmin = 0; rmax = v(4); rticks = max(ticks-1,2);
    if rticks > 5   % see if we can reduce the number
        if rem(rticks,2) == 0
            rticks = rticks/2;
        elseif rem(rticks,3) == 0
            rticks = rticks/3;
        end
    end

% define a circle
    th = 0:pi/50:2*pi;
    xunit = cos(th);
    yunit = sin(th);
% now really force points on x/y axes to lie on them exactly
    inds = 1:(length(th)-1)/4:length(th);
    xunit(inds(2:2:4)) = zeros(2,1);
    yunit(inds(1:2:5)) = zeros(3,1);
% plot background if necessary
    if ~isstr(get(cax,'color')),
       patch('xdata',xunit*rmax,'ydata',yunit*rmax, ...
             'edgecolor',tc,'facecolor',get(gca,'color'),...
             'handlevisibility','off');
    end

% draw radial circles
    c82 = cos(82*pi/180);
    s82 = sin(82*pi/180);
    rinc = (rmax-rmin)/rticks;
    icount=1;
    for i=(rmin):rinc:rmax
        if icount > 1
        hhh = plot(xunit*i,yunit*i,ls,'color',tc,'linewidth',1,...
                   'handlevisibility','on');
        end
        if icount < 6      
        text((i+(rinc)/20)*c82,(i+rinc/20)*s82, ...
           ['  ' num2str(90-i),'\circ'],'verticalalignment','bottom',...
           'handlevisibility','on')
        end
       icount=icount+1;
    end
    set(hhh,'linestyle','-') % Make outer circle solid

% plot spokes
    th = (1:6)*2*pi/12;
    cst = sin(th); snt = cos(th);
    cs = [-cst; cst];
    sn = [-snt; snt];
    plot(rmax*cs,rmax*sn,ls,'color',tc,'linewidth',1,...
         'handlevisibility','on')

% annotate spokes in degrees
    rt = 1.1*rmax;
    for i = 1:length(th)
        text(rt*cst(i),rt*snt(i),[int2str(i*30),'\circ'],...
             'horizontalalignment','center',...
             'handlevisibility','on');
        if i == length(th)
            loc = int2str(0);
        else
            loc = int2str(180+i*30);
        end
        text(-rt*cst(i),-rt*snt(i),[loc,'\circ'],'horizontalalignment','center',...
             'handlevisibility','on')
    end

% set view to 2-D
    view(2);
% set axis limits
    axis(rmax*[-1 1 -1.15 1.15]);
end

% Reset defaults.
% set(cax, 'DefaultTextFontAngle', fAngle , ...
%     'DefaultTextFontName',   fName , ...
%     'DefaultTextFontSize',   fSize, ...
%     'DefaultTextFontWeight', fWeight, ...
%     'DefaultTextUnits',fUnits );

% transform data to Cartesian coordinates.
xx = rho.*sin(theta);
yy = rho.*cos(theta);

% plot data on top of grid
colormap('jet');
%plotclr(xx,yy,z,prn);
pointplot2(xx,yy,z,'o',10)
hold on
for ii=1:length(prn)
    if maxii(ii)~=0
        if prn(ii)<10
text(xx(maxii(ii),ii),yy(maxii(ii),ii),['0',num2str(prn(ii))],'fontsize',22,'fontweight','bold')
        else
text(xx(maxii(ii),ii),yy(maxii(ii),ii),['',num2str(prn(ii))],'fontsize',22,'fontweight','bold')            
        end    
    end
%set(ht,fontsize,8,'verticalalignment','bottom',...
%           'handlevisibility','on')
end
if ~hold_state
    set(gca,'dataaspectratio',[1 1 1]), axis off; set(cax,'NextPlot',next);
end
set(get(gca,'xlabel'),'visible','on')
set(get(gca,'ylabel'),'visible','on')

cosmetique(12)
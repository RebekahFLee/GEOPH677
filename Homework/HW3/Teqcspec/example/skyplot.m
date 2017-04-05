close all
%clear all
clear z
% event='SAMOA_092909';
% cd('~/SKYPLOTS/SAMOA_090929/')
% if (1)
% out1=teqcplot('mkea2730.ele');
% out2=teqcplot('mkea2730.azi');
% end
% satnum=[1:31];
% %satnum=[12 2 9 21 5 26 15 27 18 24 12 10 29];
% seismepoch=17+48/60;
% ti=0:60/3600:3;

% event='TOK-OKI_030925';
% cd('~/SKYPLOTS/TOK-OKI_030925/')
% if (1)
% out1=teqcplot('0490268s.ele');
% out2=teqcplot('0490268s.azi');
% end    
% satnum=[7 24 27 4 29 13 17 10]
%  seismepoch=19+50/60;
%  %ti=19+50/60+60/3600:60/3600:20+50/60;
% ti=19+30/60+60/3600:60/3600:21;
% event='WENSHUAN_080512';
% cd(['~/SKYPLOTS/',event])
% if (1)
% out1=teqcplot('gmsd1330.ele');
% out2=teqcplot('gmsd1330.azi');
% end    
% satnum=[9 5 9 12 14 15 18 21 22 30];
%  seismepoch=6+28/60;
 %ti=6+28/60+60/3600:60/3600:8+28/60;
%ti=6+60/3600:60/3600:8;
% event='NZ_090715';
% cd(['~/SKYPLOTS/',event])
% if (1)    
% out1=teqcplot('nlsn1960.ele');
% out2=teqcplot('nlsn1960.azi');
% end    
% satnum=[1 16 30 31 22 12 11 14 20];
%  seismepoch=9+22/60+30/3600;
%  ti=seismepoch+60/3600:60/3600:seismepoch+1-30/3600;

% event='HAITI_100112';
% cd(['~/SKYPLOTS/',event])
% if (1)
% out1=teqcplot('scub0120.ele');
% out2=teqcplot('scub0120.azi');
% end    
% satnum=[10 27 26 18 24 29 9 22 15 5 21];
% seismepoch=21+53/60;
% ti=seismepoch+60/3600:60/3600:seismepoch+1-30/3600;

% event='SUMATRA_041226';
% cd(['~/SKYPLOTS/',event])
% if (1)
% out1=teqcplot('sa043610.ele');
% out2=teqcplot('sa043610.azi');
% end    
% satnum=[10 3 23 21 25 16 15 1 20 19 13];
% seismepoch=0+59/60;
% ti=seismepoch+60/3600:60/3600:seismepoch+3+30/3600;



event='CHILE_071114';
cd(['~/SKYPLOTS/',event])
if (1)
out1=teqcplot('ctlr3180.ele');
out2=teqcplot('ctlr3180.azi');
end    
satnum=[1:32]
 seismepoch=15+40/60;
 %ti=19+50/60+60/3600:60/3600:20+50/60;
ti=15+30/60+2*60/3600:2*60/3600:17+30/60;

rho=abs(out1.ele-90);
theta=(out2.azi)*pi/180;

for i=1:32
z(:,i)=ti';
if length(find(isnan(rho(:,i)))) > length(rho)/2
rho(:,i)=NaN;    
end    
if ~isempty(find(rho(:,i)<30*pi/180))
    theta(find(rho(:,i)<30*pi/180))=NaN;
    rho(find(rho(:,i)<30*pi/180))=NaN;
end    
if mean(rho(find(~isnan(rho(:,i))),i)) > 80 
rho(:,i)=NaN;    
end
if isempty(intersect(i,satnum))
rho(:,i)=NaN; 
end    
rho2(:,i)=downsample(rho(:,i),5);
theta2(:,i)=downsample(theta(:,i),5);
z2(:,i)=downsample(z(:,i),5);

if ~isempty(find(~isnan(rho2(:,i))))
maxii(i)=max(find(~isnan(rho2(:,i))));
else
maxii(i)=0;    
end

end



prn=[1:32];
scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(2) scrsz(3) scrsz(4)/1])
%figure
%mmpolar(theta,rho,'Rlimit',[0 80])
cmin=15.5;
cmax=17.5;
polarview(theta2,rho2,z2,prn,maxii);
colormap('jet')
h=colorbar('location','Southoutside')
ctitle(h,'Time (GMT)',18)
%caxis([21 22])
%title('Skyplot 03/09/25 19.5h-21h GMT','fontsize',14)
cosmetique(18)
set(gca,'fontname','Helvetica')
caxis([15.5 17.5])
print(gcf, '-depsc',['~/ARTICLES/2009/JGR/FIGURES/Skyplot',event,'2big.eps'])
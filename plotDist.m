function plotDist(axd)
% plotDist(axd)
% 
% description:
%   Plots a distribution of axons into gca
%
% inputs:
%   axd - structure containing vectors xr,yr,ri,ro, and scalar Lx
% 
% by Kevin Harkins (kevin.harkins@vanderbilt.edu)

rectangle('Position',[0 0 axd.Lx axd.Lx],'facecolor',[1 1 1],'linestyle','none');
hold on
for n=1:numel(axd.xr)
    rectangle('Position',[ [axd.xr(n) axd.yr(n)]-axd.ro(n) [axd.ro(n) axd.ro(n)]*2],...
        'facecolor',[0 0 0],'curvature',[1 1],'linestyle','none')
    rectangle('Position',[ [axd.xr(n) axd.yr(n)]-axd.ri(n) [axd.ri(n) axd.ri(n)]*2],...
        'facecolor',[.75 .75 .75],'curvature',[1 1],'linestyle','none')
end
ylim([0 axd.Lx])
xlim([0 axd.Lx])
axis square off
hold off

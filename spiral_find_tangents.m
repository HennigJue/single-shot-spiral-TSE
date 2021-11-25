function [tanxy,itanxy]=spiral_find_tangents(kt,sp,nc,pflag)
%%
% finds all tangential lines to spiral kt through point sp. sp musst be on
% x-axis.
%sp=[100 0]; 
if (nargin<4),
    pflag=0;
end
%% 
% draw circle through sp,[0 0]
r=abs(sp(1)/2);
cr=linspace(0,2*pi,nc);
circ(:,1)=r*(cos(cr));
circ(:,2)=r*(sin(cr));
circ(:,1)=circ(:,1)+sp(1)/2;
% intersect circle with spiral
[xsc,ysc,ii]=polyxpoly(kt(:,1),kt(:,2),circ(:,1),circ(:,2));
% find points on spiral closest to intersection
tanxy(:,1)=kt(ii(:,1),1);
tanxy(:,2)=kt(ii(:,1),2);
itanxy=ii(:,1);
ntan=floor(length(itanxy)/2);
itanxy_1=sort(itanxy(1:ntan));
itanxy_2=sort(itanxy(ntan+1:end));
itanxy=[itanxy_1; itanxy_2];
if(pflag>0),
figure
plot(kt(:,1),kt(:,2))
axis square
hold on
plot(circ(:,1),circ(:,2),'r-');
plot(xsc,ysc,'go')
plot(kt(ii(:,1),1),kt(ii(:,1),2),'ko')
end
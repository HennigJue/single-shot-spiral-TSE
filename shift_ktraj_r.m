
function [kt_s]=shift_ktraj_r(kt,dt,plotflag)
%% shift k-space trajectory
% consisting of trajectory as real numbered matrix kt=ktraj(:,[x y z])
%close all
if(nargin<3),plotflag=0; end;
skt=size(kt);
t=1:skt(2);
kt_s=0*kt;


for k=1:skt(1),
kt_s(k,:) = interp1(t,kt(k,:),t+dt(k),'spline');
end

if (plotflag==1),
    figure
    hold on
    plot(kt,'.-')
    plot(kt_s,'r.-')
end


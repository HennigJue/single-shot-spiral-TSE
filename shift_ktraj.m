
function [kt_s]=shift_ktraj(kt,dx,dy,plotflag)
%% shift k-space trajectory
close all
if(nargin<4),plotflag=0; end;
skt=size(kt);
t=1:skt(1);
ktx=real(kt);ktx_s=0*ktx;
kty=imag(kt);kty_s=0*kty;

for k=1:skt(2),
ktx_s(:,k) = interp1(t,ktx(:,k),t+dx,'spline');
kty_s(:,k) = interp1(t,kty(:,k),t+dy,'spline');
end
kt_s=ktx_s+i*kty_s;
if (plotflag==1),
    figure
    hold on
    plot(kt,'.-')
    plot(kt_s,'r.-')
end


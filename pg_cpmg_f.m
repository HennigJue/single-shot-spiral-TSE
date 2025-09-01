% function mph = pg_cpmg(m0,kper,rflip,rphas,te,r1,r2);
%
% calculates magnetization at echo times for RARE-sequence with different flip angles
%
%	requires:
%	m0: magnetization vector for initialization (=m after excitation pulse): 
%       m0=[1 0 0 0]' for pure x-magnetization
%	
%	kper: 	number of periodic refocusings (=number of echoes/period length
%	rflip: 	flip angle of refocusing pulses in echo train (array)
%	rphas: 	array containing phase values for ref.pulse. If size(rphas).ne.size(rflip)->rphas=0.
%	te:		list of te (in ms symmetrical around pulse. default: 1.
%	r1,r2:	relaxivities.
%
function [mph,mpe,pg] = pg_cpmg_f(m0,kper,rflip,rphas,te,r1,r2,plotflag);
if (nargin<8),
    plotflag=2;
end;

if (nargin<7),
   r2=0;
end;
if (nargin<7),
   r1=0;
end;
if (nargin<5),
   te=([1]);
end;
if (nargin<4),
   rphas=([0]);
end;

winc=1000;
rflip;
sm=size(m0); sp=size(rflip);
sph=size(rphas);
sw=size(te);
if (sph(2)==1),
   rphas=rphas*ones([sp(1) sp(2)]);
end;
if (sw(2)==1),
   te=te*ones([sp(1) sp(2)]);
end;
sph=size(rphas);
sw=size(te);
if (sph(2)~=sp(2)),
   rphas=zeros([sp(1) sp(2)]);
end;
if (sw(2)~=sp(2)),
   te=ones([sp(1) sp(2)])
end;
rflip; 
rphas;
te;
   nper=sp(2);
nech=kper*nper;
mph=zeros([4 nech+1+sm(2) nech+1+sm(2)]);

mph(1:sm(1),1:sm(2),1)=m0;
mpe=mph;
size(mph);
k=2;


for m=1:kper,
for n=1:nper;
phi=rphas(n);
phi;
rphas(n);
mph(:,:,k-1)=pg_pulse(mpe(:,:,k-1),rflip(n),rphas(n));
k;
size(mph);
mpe(:,:,k)=pg_evo_cpmg(mph(:,:,k-1),te(n),r1,r2);

k=k+1;
end;
end;
%mph(:,:,k-1)=pg_pulse(mph(:,:,k-1),rflip(n),rphas(n));

pg=pg_comb(mph);
%pg(1,1,1)=1;pg(2,1,1)=-0.5;


smp=size(mph);
k;
x=1:(k-1);

if (plotflag==1),
figure 
plot(x,squeeze(abs(mph(2,1,1:k-1))),'o-'),axis([0,k-2,0,1]),xlabel('# of echoes'),ylabel('intensity');
end;

if (plotflag==2),
figure
set(gcf,'Position',[200 0 800 300]);
subplot(1,2,1), plot(x,squeeze(abs(mph(2,1,1:k-1)))),axis([0,k-2,0,1]),xlabel('# of echoes'),ylabel('intensity');
subplot(1,2,2), imagesc(abs(pg(:,1:(k-2),1)),[-0.5 1]),colormap('jet');
end

if (plotflag==3)
figure
set(gcf,'Position',[200 0 1200 300]);
subplot(1,3,1), plot(x,squeeze(abs(mph(2,1,1:k-1)))),axis([0,k-2,0,1]),xlabel('# of echoes'),ylabel('intensity');
subplot(1,3,2), imagesc(abs(pg(:,1:(k-2),1)),[-0.5 1]),colormap('jet');
subplot(1,3,3), imagesc(abs(pg(:,1:(k-2),2)),[-0.5 1]),colormap('jet');

end


if (plotflag==4),
figure
%set(gcf,'Position',[200 0 200 400]);
imagesc(real(pg(:,1:(k-2),1)),[-0.5 1]),colormap(jet);;
% ,colormap('gray');
end;

%subplot(1,2,2), imagesc(real(pg(:,:,1)));
%title(strcat(' flip =',num2str(rflip),' int=',num2str(amp(k-1,1))));



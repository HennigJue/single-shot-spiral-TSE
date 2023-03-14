function [Graster, traster, Ginit, tinit] = spiral_init2sp(kBegin, kEnd, G, k, tsp, dG, system);
%%
tflip=0;
dK=kEnd-kBegin;
m=dK(2)/dK(1);
nxy=2;
if(abs(m)<=1), nxy=1; end;
nxy

cflag='1';
dKG(nxy)=G(nxy)*(tsp-dG/2);
if(dKG(nxy)>=dK(nxy)), cflag='2';
    dKG(nxy)=G(nxy)*dG+G(nxy)*(tsp-2*dG)/2;
    if(dKG(nxy)>=dK(nxy)), cflag='3';
        dKG(nxy)=G(nxy)*dG/2;
    end
end
cflag;
%end

switch cflag
    case '1',
        Gamp=(dK(nxy)-G(nxy)*dG/2)/(tsp-dG);
        Ginit=[0 Gamp Gamp G(nxy)];
        
    case '2',
        Gamp=(dK(nxy)-G(nxy)*dG)/(tsp/2-dG/2);
        Ginit=[0 Gamp G(nxy) G(nxy)];
    case '3',
        Gamp=(dK(nxy)-G(nxy)*dG/2)/(tsp/2-dG/2);
        Ginit=[0 0 Gamp G(nxy)];    
end
tinit=[0 dG tsp-dG tsp];
traster=0:system.gradRasterTime:tsp;
Graster=zeros([length(traster) 2]);
Graster(:,nxy)=interp1(tinit,Ginit,traster);
if (nxy==1), nyx=2; Graster(:,nyx)=Graster(:,nxy)*m; end;
if (nxy==2), nyx=1; Graster(:,nyx)=Graster(:,nxy)/m; end;

%Graster(:,nxy)=Graster(:,nxy).*dK(nxy)/(mean(Graster(:,nxy))*tsp);



% mean(Graster)*tsp;
% figure
% plot(tinit,Ginit,'ro')
% hold on
% plot(traster,Graster,'b.')

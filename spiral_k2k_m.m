function [Graster, traster, Gtran, tG] = spiral_k2k_m(kBegin, kEnd, GBegin, GEnd, tsp, slew, dW);
%%
dK=kEnd-kBegin;
dK1=0*dK; dKr=dK1;
t1=0*kBegin;t2=t1; tamp=t1;Gamp=t1;
k=1;


B=-(GBegin(k)+GEnd(k)+tsp*slew);
C=(GBegin(k).^2+GEnd(k).^2)/2+dK(k)*slew;
Gampq1=-B/2+sqrt((B/2).^2-C);
Gampq2=-B/2-sqrt((B/2).^2-C);

t11temp=((Gampq1-GBegin(k))/slew);
t12temp=((Gampq1-GEnd(k))/slew);
tamp1=tsp-t11temp-t12temp;

t21temp=((Gampq2-GBegin(k))/slew);
t22temp=((Gampq2-GEnd(k))/slew);
tamp2=tsp-t21temp-t22temp;

if(tamp1>0&t11temp>0&t12temp>0)
    t1(k)=abs(t11temp); t2(k)=abs(t12temp); Gamp(k)=Gampq1;
else if(tamp2>=0&t21temp>0&t22temp>0)
    t1(k)=abs(t21temp); t2(k)=abs(t22temp); Gamp(k)=Gampq2;
    else
        if(GBegin(k)<=GEnd(k)),
        Gamp(k)=(dK(k)*slew+1/2*GBegin(k).^2-1/2*GEnd(k).^2)/(slew.*tsp+GBegin(k)-GEnd(k));
        t1(k)=(Gamp(k)-GBegin(k))./slew;t2(k)=(GEnd(k)-Gamp(k))./slew;
        else
        Gamp(k)=(dK(k)*slew-1/2*GBegin(k).^2+1/2*GEnd(k).^2)/(slew.*tsp-GBegin(k)+GEnd(k));
        t1(k)=(-Gamp(k)+GBegin(k))./slew;t2(k)=(-GEnd(k)+Gamp(k))./slew;
        end
            
    end
end
tamp(k)=tsp-t1(k)-t2(k);
dK1(k)=(GBegin(k)+Gamp(k)).*t1(k)/2+Gamp(k).*tamp(k)+(GEnd(k)+Gamp(k)).*t2(k)/2;




k=2;

B=-(GBegin(k)+GEnd(k)+tsp*slew);
C=(GBegin(k).^2+GEnd(k).^2)/2+dK(k)*slew;
Gampq1=-B/2+sqrt((B/2).^2-C);
Gampq2=-B/2-sqrt((B/2).^2-C);

t11temp=((Gampq1-GBegin(k))/slew);
t12temp=((Gampq1-GEnd(k))/slew);
tamp1=tsp-t11temp-t12temp;

t21temp=((Gampq2-GBegin(k))/slew);
t22temp=((Gampq2-GEnd(k))/slew);
tamp2=tsp-t21temp-t22temp;

if(tamp1>0&t11temp>0&t12temp>0)
   % disp('1')
    t1(k)=abs(t11temp); t2(k)=abs(t12temp); Gamp(k)=Gampq1;
else if(tamp2>=0&t21temp>0&t22temp>0)
    % disp('2')
        t1(k)=abs(t21temp); t2(k)=abs(t22temp); Gamp(k)=Gampq2;
    else
        if(GBegin(k)<=GEnd(k)),
        %disp('3')
            Gamp(k)=(dK(k)*slew+1/2*GBegin(k).^2-1/2*GEnd(k).^2)/(slew.*tsp+GBegin(k)-GEnd(k));
        t1(k)=(Gamp(k)-GBegin(k))./slew;t2(k)=(GEnd(k)-Gamp(k))./slew;
        else
        %disp('4')    
        Gamp(k)=(dK(k)*slew-1/2*GBegin(k).^2+1/2*GEnd(k).^2)/(slew.*tsp-GBegin(k)+GEnd(k));
        t1(k)=(-Gamp(k)+GBegin(k))./slew;t2(k)=(-GEnd(k)+Gamp(k))./slew;
        end
            
    end
end
tamp(k)=tsp-t1(k)-t2(k);
dK1(k)=(GBegin(k)+Gamp(k)).*t1(k)/2+Gamp(k).*tamp(k)+(GEnd(k)+Gamp(k)).*t2(k)/2;

% figure
% k=1;
% plot([0 t1(k) t1(k)+tamp(k) t1(k)+tamp(k)+t2(k)],[GBegin(k) Gamp(k) Gamp(k) GEnd(k)],'b-')
% hold on
% k=2;
% plot([0 t1(k) t1(k)+tamp(k) t1(k)+tamp(k)+t2(k)],[GBegin(k) Gamp(k) Gamp(k) GEnd(k)],'r-')

Gtran=[GBegin; Gamp; Gamp; GEnd];
tG=[[0 0];t1;t1+tamp;t1+tamp+t2];
traster=dW/2:dW:tsp;
Graster=zeros([length(traster) 2]);
Gintp=interp1(tG(:,1),Gtran(:,1),traster);
Graster(:,1)=Gintp;
Gintp=interp1(tG(:,2),Gtran(:,2),traster);
Graster(:,2)=Gintp;
dKr=sum(Graster)*dW;





    

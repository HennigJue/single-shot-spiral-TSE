function [GPoi,tPoi,t] = spiral_k2k_min(kBegin, kEnd, GBegin, GEnd, Gmax, slew);
%% [GRastery, tRastery, GPoiy, tPoiy, touty]
% if(nargin<8),
%     tmin=0;
% end

G1=GBegin;
G2=GEnd;
K1=kBegin;
K2=kEnd;
dK=K2-K1;
if(dK<0),Gmax=-Gmax; 
end
t1=abs((Gmax-G1)/slew);
dK1=(Gmax+G1)*t1/2;
t2=abs((Gmax-G2)/slew);
dK2=(Gmax+G2)*t2/2;
tm=(dK-dK1-dK2)/Gmax;
t=tm+t1+t2;

tPoi=[0 t1 t1+tm t1+tm+t2];
GPoi=[GBegin Gmax Gmax GEnd];

if(tm<0),
    t12=abs((G2-G1)/slew);
    if(G2>G1)
        ddK=(G1+G1+t12*slew)/2;
    else
        ddK=(G1+G1-t12*slew)/2;
    end
    if(dK<ddK), slew=-slew; end
    %slew=-slew;
    A=1/4*slew;
    B=1/2*(GEnd+GBegin);
    C=-(dK+1/4*(GEnd-GBegin)^2/slew);
    ttot(1)=(-B+sqrt(B.^2-4*A*C))/(2*A);
    ttot(2)=(-B-sqrt(B.^2-4*A*C))/(2*A);
     
    
%    if(isreal(ttot)),
    ind=find(ttot>0);
    t=ttot(ind);
    %t=t(1);
    t1=1/2*(GEnd-GBegin)/slew+1/2*t;
    t2=t-t1;
    if(length(t1)>1),
        if((t1(1)>=0)&(t2(1)>=0)), kk=1;
        else
            kk=2;
        end
        t1=t1(kk); t2=t2(kk);
    end
    
    Gm=GBegin+slew*t1;

tPoi=[0 t1 t1+t2];
GPoi=[GBegin Gm GEnd];
end

% 
% figure
% plot(tPoi,GPoi); 
%axis([0 25 -abs(Gmax) abs(Gmax)])



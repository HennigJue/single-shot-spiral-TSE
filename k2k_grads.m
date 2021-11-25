function [GRaster, tRaster, GPoi, tPoi, tout] = k2k_grads(kBegin, kEnd, GBegin, GEnd, Gmax, slew,t,tGrast,pflag);
% %%
if (nargin<9),
    pflag=0;
end

dK=kEnd-kBegin;
if((abs(GBegin-GEnd)/slew)>t)
    t=abs(ceil(1/tGrast*(GBegin-GEnd)/slew))*tGrast;
end

if(dK==(GBegin+GEnd)/2*t),
    tPoi=[0 t];
    GPoi=[GBegin GEnd];
    tRaster=tGrast/2:tGrast:t;
    GRaster=interp1(tPoi,GPoi,tRaster);
    tout=t;
    return
end


smode=strings([5,1]);
smode(1:5)='trap';
t11=ceil(1/tGrast*(Gmax-GBegin)/slew)*tGrast;
t21=ceil(1/tGrast*(Gmax-GEnd)/slew)*tGrast;
tm1=t-t11-t21;
t12=ceil(1/tGrast*(Gmax+GBegin)/slew)*tGrast;
t22=ceil(1/tGrast*(Gmax+GEnd)/slew)*tGrast;
tm2=t-t12-t22;
tG=ceil(1/tGrast*abs(GBegin-GEnd)/slew)*tGrast;
if(tm1==0),
 %   disp('hallo')
end
if(tm1>=0),
    K1=(GBegin+Gmax)/2*t11+Gmax*tm1+(GEnd+Gmax)/2*t21;
else
    smode(1)='tri';
    smode(2)='tri';
    Gmaxt=(GBegin+t*slew+GEnd)/2;
    t11=floor(1/tGrast*(Gmaxt-GBegin)/slew)*tGrast;
    Gmaxt=GBegin+t11*slew;
    t21=t-t11; tm1=0;
    K1=(GBegin+Gmaxt)/2*t11+(GEnd+Gmaxt)/2*t21;
end

K2=(GBegin+GEnd)/2*tG+GBegin*(t-tG);
K3=(GBegin+GEnd)/2*tG+GEnd*(t-tG);

if(tm2>=0),
    K4=(GBegin-Gmax)/2*t12-Gmax*tm2+(GEnd-Gmax)/2*t22;
else
    smode(4)='tri';
    smode(5)='tri';
    Gmint=(GBegin-t*slew+GEnd)/2;
    t12=floor(1/tGrast*(GBegin-Gmint)/slew)*tGrast;
    Gmint=GBegin-t12*slew;
    t22=t-t12; tm2=0;
    K4=(GBegin+Gmint)/2*t12+(GEnd+Gmint)/2*t22;
end

Kall=[K1 K2 K3 K4];
smode;

ind=find(dK==Kall);
if(ind),
    gmode=strcat('i',num2str(ind));
else
    Kall=[dK Kall];
    Ksort=sort(Kall,'descend');
    ind=find(Ksort==dK);
    gmode=strcat('n',num2str(ind));
end

if(pflag==1),
    disp(Kall);
    disp(smode);
    disp(gmode);
end


G=0;
switch gmode
    case 'i1',
        G=Gmax;t1=t11; t2=t21; tm=tm1;
        if(smode(ind)=='tri'), G=Gmaxt; end
    case 'i2',
        G=GBegin; t1=tG; t2=tG; tm=t-2*tG;
        if(smode(ind)=='tri'), G=GBegin;t1=t-tG; t2=tG;tm=0; end
    case 'i3',
        G=GEnd; t1=tG; t2=tG; tm=t-2*tG;
    case 'i4',
        G=-Gmax;t1=t12; t2=t22; tm=tm2;
        if(smode(ind)=='tri'), G=Gmint; end
        
    case 'n1',
        t1=t11; t2=t21; tm=tm1;
        ddK=dK-K1;
        if(smode(1)=='trap'),
            tm=tm+ceil(1/tGrast*ddK/Gmax)*tGrast;
            G=(dK-GBegin*t1/2-GEnd*t2/2)/(t1/2+tm+t2/2);
            t=t1+tm+t2;
        end
        if(smode(1)=='tri'),
            t1m=ceil(1/tGrast*(Gmax-GBegin)/slew)*tGrast;
            t2m=ceil(1/tGrast*(Gmax-GEnd)/slew)*tGrast;
            K1m=(GBegin+Gmax)/2*t1m+(GEnd+Gmax)/2*t2m;
            if(dK>K1m),
                t1=t1m; t2=t2m;
                ddK=dK-K1m;
                tm=ceil(1/tGrast*ddK/Gmax)*tGrast;
                G=(dK-GBegin*t1/2-GEnd*t2/2)/(t1/2+tm+t2/2);
                t=t1+tm+t2;
                smode(1)='trap';
            else
                A=slew;
                B=2*GBegin;
                C=(GBegin^2-GEnd^2)/(2*slew)-dK;
                t1=(-B+sqrt(B^2-4*A*C))/(2*A);
                t1=ceil(1/tGrast*t1)*tGrast;
                G=GBegin+t1*slew;
                t2=(G-GEnd)/slew;
                t2=ceil(1/tGrast*t2)*tGrast;
                t=t1+t2;
                G=(2*dK-GBegin*t1-GEnd*t2)/(t1+t2);
            end
        end
    case 'n2',
        t1=t11; t2=t21; tm=tm1;
        G=(dK-GBegin*t1/2-GEnd*t2/2)/(t1/2+tm+t2/2);
        if(smode(2)=='tri'),
            A=slew;
            B=2*GBegin;
            C=(GBegin^2-GEnd^2)/(2*slew)-dK;
            t1=(-B+sqrt(B^2-4*A*C))/(2*A);
            t1=ceil(1/tGrast*t1)*tGrast;
            G=GBegin+t1*slew;
%             t2=(G-GEnd)/slew;
%             t2=ceil(1/tGrast*t2)*tGrast;
            t2=t-t1;
            G=(2*dK-GBegin*t1-GEnd*t2)/(t1+t2);
        end
        
    case 'n3',
        t1=tG; t2=tG; tm=t-2*tG;
        G=(dK-GBegin*tG/2-GEnd*tG/2)/(t-tG);
        if(tm<0),
            t1=tG;t2=t-tG;tm=0;
            G=2*(dK-GBegin/2*t1-GEnd/2*t2)/t;         
        end
    
  
    
      
    
    
    case 'n4',
        t1=t12; t2=t22; tm=tm2;
        G=(dK-GBegin*t1/2-GEnd*t2/2)/(t1/2+tm+t2/2);
        if(smode(4)=='tri'),

 G=2*(dK-GBegin/2*t1-GEnd/2*t2)/t;     
        end
    case 'n5',
        t1=t12; t2=t22; tm=tm2;
        ddK=K4-dK;
        if(smode(5)=='trap'),
            tm=tm+abs(ceil(1/tGrast*ddK/Gmax))*tGrast;
            G=(dK-GBegin*t1/2-GEnd*t2/2)/(t1/2+tm+t2/2);
            t=t1+tm+t2;
        end
        
        if(smode(5)=='tri'),
            t1m=ceil(1/tGrast*(Gmax+GBegin)/slew)*tGrast;
            t2m=ceil(1/tGrast*(Gmax+GEnd)/slew)*tGrast;
            K1m=(GBegin-Gmax)/2*t1m+(GEnd-Gmax)/2*t2m;
            if(dK<K1m),
                t1=t1m; t2=t2m;
                ddK=dK-K1m;
                tm=abs(ceil(1/tGrast*ddK/Gmax)*tGrast);
                G=(dK-GBegin*t1/2-GEnd*t2/2)/(t1/2+tm+t2/2);
                t=t1+tm+t2;
                smode(5)='trap';
            else
                A=slew;
                B=-2*GBegin;
                C=(GBegin^2-GEnd^2)/(2*slew)+dK;
                t1=(-B+sqrt(B^2-4*A*C))/(2*A);
                t1=ceil(1/tGrast*t1)*tGrast;
                G=GBegin-t1*slew;
                t2=(GEnd-G)/slew;
                t2=ceil(1/tGrast*t2)*tGrast;
                t=t1+t2;
                G=(2*dK-GBegin*t1-GEnd*t2)/(t1+t2);
            end
        end
        
        
end

if(tm>0),
    tPoi=[0 t1 t1+tm t1+tm+t2];
    GPoi=[GBegin G G GEnd];
else
    tPoi=[0 t1 t1+t2];
    GPoi=[GBegin G GEnd];
end
tout=t;


tRaster=tGrast/2:tGrast:t;
GRaster=zeros([length(tRaster) 1]);
GRaster=interp1(tPoi,GPoi,tRaster);
dKr=sum(GRaster)*tGrast;

if (pflag==1),
figure
plot(tPoi,GPoi,'r-'),axis([0 t -Gmax Gmax]);
hold on
plot(tRaster,GRaster,'b.'),axis([0 t -Gmax Gmax]);
end


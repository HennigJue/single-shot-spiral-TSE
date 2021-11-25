function [ktacq,indkt,ksel]=spiral_stitch_nom(ktraj_adc,i1,nkseg,nadc,necho,nTE,dninc,accfac);
col=['b.';'go'];
figure
hold on
indkt=zeros([necho 4]);
ks=1+mod([1:necho]-nTE,necho);
ktacq=zeros([2 nadc*necho]);
kti=squeeze(ktraj_adc(1,:,:)+i*ktraj_adc(2,:,:));
ksel=zeros([nadc necho]);
inda=0; inde=0;
dk=9/accfac;

for kn=1:necho,
    k=ks(kn);
    indk=[accfac*(i1(ks(kn))-dninc)+1 accfac*i1(ks(kn))+accfac*(nkseg(ks(kn))-dninc-1)];
    indkt(kn,1:2)=indk(1:2);
    plot(ktraj_adc(1,indk(1):indk(2),kn,1),ktraj_adc(2,indk(1):indk(2),kn,1),col(1,:))
    col=flipud(col);% plot k-st segment
end
% figure
% hold on
indkt(:,4)=indkt(:,2);
indkt(:,3)=indkt(:,1);
for kn=1:necho,
   
    k=ks(kn);
    indk=[accfac*(i1(ks(kn))-dninc)+1 accfac*i1(ks(kn))+accfac*(nkseg(ks(kn))-dninc-1)];
    indkt(kn,1:2)=indk(1:2);
    temp1=[ktraj_adc(1,indk(1):indk(2),kn)];
    temp2=[ktraj_adc(2,indk(1):indk(2),kn)];
    plot(ktraj_adc(1,indk(1):indk(2),kn,1),ktraj_adc(2,indk(1):indk(2),kn,1),col(1,:))  % plot k-st segment
    
    if(kn<necho-1),
        
        indt=find(abs(kti(:,kn)-kti(indkt(kn+1,1),kn+1))<dk);  % find k-space point identical to the last point in the previous segment)
        try, indkt(kn,4)=indt(1);
            plot(kti(indkt(kn,4),kn),'rd')
            plot(kti(indkt(kn+1,1),kn+1),'kd')
        end
    end
    
    if(kn>1),
        indt=find(abs(kti(:,kn)-kti(indkt(kn-1,2),kn-1))<dk);  % find k-space point identical to the last point in the previous segment)
        try, indkt(kn,3)=indt(1);
            plot(kti(indkt(kn,3),kn),'rd')
            plot(kti(indkt(kn-1,2),kn-1),'kd')
        end
    end
    
    ksel(indk(1):indk(2),kn)=1;
    inde=length(temp1);
    ktacq(1,inda+1:inda+inde)=temp1;
    ktacq(2,inda+1:inda+inde)=temp2;
    inda=inda+inde;
    col=flipud(col);
end

ktacq=ktacq(:,1:inda);

%
% ktacq=reshape(ktraj_adc,[3,nadc*necho]);
% ktacq=ktacq(1:2,indk);
% figure
% plot(ktacq(1,:),ktacq(2,:),'b.-')

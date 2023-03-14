%% plot acquisition scheme
nslices_s=16;
ntp=nslices/nslices_s;
indm=reshape(ind,[nslices/nslices_s nslices_s]);
figure
plot(slOrder,'o-')
hold on
for k=1:nslices_s
    for m=2:2:ntp-1,
    plot(indm(m:m+1,k),slOrder(indm(m:m+1,k)),'r.-')
    end
    for m=1:2:ntp-2,
    plot(indm(m:m+1,k),slOrder(indm(m:m+1,k)),'k.-')
    end
end
%% resort images according to TR
imt1_trs=zeros([240 240 9 16]);
for k=1:16
tr=diff(indm(:,k));
[temp,trind]=sort(tr,'descend');
trs(:,k)=temp;
imt1_trs(:,:,:,k)=imall(:,:,(k-1)*ntp+[1 trind(1:end)'+1]);
bigx=im_mosaic(squeeze(imt1_trs(:,:,:,k)),1,5,150);
end





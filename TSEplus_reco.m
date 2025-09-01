% very basic and crude non-cartesian recon using griddata()
%%

%rawDir = '/raid/groupspace/range/rawdata';
%addpath('/rktraid/home/extern/range/code/mapvbvd')
%cd(pnp)
clear fnr
disp('select rawdata')
[fnr, pnr] = uigetfile('*.dat', 'Select rawdata files','MultiSelect', 'on');
if(iscell(fnr)), nir=numel(fnr); else nir=1; fnr={fnr}; end
%%
clear fns
disp('select sequence(s)')
[fns, pns] = uigetfile('*.seq', 'Select protocol files','MultiSelect', 'on');
if(iscell(fns)), nip=numel(fns); else nip=1; fns={fns}; end
%
%pnp=pns;
for k=1:length(fns)
    temp=cell2mat(fns(k));
    fnp(k)={strcat('p',temp(4:end-3),'mat')};
end

%%
for k=1:length(fas)
    temp=cell2mat(fas(k));
    %fnp(k)={strcat('p',temp(1:end-3),'mat')};
    fnp(k)={strcat(temp,'.mat')};
end
%%

disp('select protocol(s)')
[fnp, pnp] = uigetfile('*.mat', 'Select protocol files','MultiSelect', 'on');
if(iscell(fnp)), nip=numel(fnp); else nip=1; fnp={fnp}; end
%cd(pnr)


%% Load datata
%kim=1;
% kipnp='D:\matlab\si(rawworJulk\pubiglsekimq_stash\myTSE_data\sequences_20210630AA';
% pnr='D:\matlab\workb\pulseq_ssizetash\myTSE_dsizeata\Spirals_AERA_1_7_2021';
%imall=zeros([240 240 15 16]);
nuFFTflag=0;
FFTflag=1;
recnorm='y';r1=0;r2=20;
mp=0;
clear imsos imcoil imall
%for kim=1:lengths(fnp)
for kim=14:15
    clear imsos imcoil
    mp=mp+1
    %         pnp=pnp0;
    %         fnp=fnp0;
    cd(pnp)
    seqname=cell2mat([fnp(kim)]);
    %         seqname=cell2mat(protname(2,kim));
    %         seqname=strcat('p',seqname(1:end-4),'_1')
    %       seqname='pTSE_T1var_on_10_20_125_144';

    load(seqname)
    try
        temp=acqP.inrep;
    catch
        acqP.inrep=1;
    end
    %acqP.nexc=12;
    %%
    %kim=29;
    cd(pnr)
    rawname=cell2mat(fnr(:,kim));
    %rawname=cell2mat(protname(1,kim));
    twix_obj = mapVBVD(rawname);
    %rawdata =double(twix_obj.image.unsorted());
    rawdata = double(twix_obj{2}.image.unsorted());
    %% resort rawdata
    raw0=permute(rawdata,[1 3 2]);
    siraw=size(raw0);
    nex=sum(acqP.nexc);
    echocount=acqP.necho;
    if(acqP.navmode(1)>0), echocount=echocount+1; end
    if(acqP.navmode(2)>0), echocount=echocount+1; end
    if(acqP.nDummy==1)
        raw=reshape(raw0,[siraw(1) echocount acqP.NSlices*acqP.inrep nex siraw(3)]);
    else
        raw=reshape(raw0,[siraw(1) echocount acqP.NSlices*acqP.inrep nex+1 siraw(3)]);
    end
    if(acqP.nDummy==0)
        rawp=raw(:,:,:,1,:);
        rawi=raw(:,:,:,2:end,:);
    else
        rawi=raw;
    end
    if(acqP.navmode(1)>0), rawnav1=rawi(:,1,:,:,:); rawi=rawi(:,2:end,:,:,:); end
    if(acqP.navmode(2)>0), rawnav2=rawi(:,end,:,:,:); rawi=rawi(:,1:end-1,:,:,:); end
    if(strcmp(recnorm,'y'))
        [mph,amp]=pg_cpmg_f([1 0 0 0],1,acqP.flip,0,acqP.TE*1000,r1,r2,0);
        signorm=squeeze(abs(mph(2,1,1:acqP.necho)));
        %signorm(1)=0.8*signorm(1);
        % if(strcmp(spmode,'alt'))
        %     signorm(acqP.PEind)=signorm;
        % end;
        disp('normalize echoes')

        for k=1:acqP.necho
            rawi(:,k,:,:,:)=rawi(:,k,:,:,:)./signorm(k);
        end

    end
    %acqP.PEorder=temp;
    PEorder=acqP.PEorder(:);
    PEorder=PEorder-min(PEorder(:))+1;
    PEcount=acqP.PEcount(:);
    pemax=find(PEorder==0);
    rawsl_all=zeros([siraw(1) acqP.necho*acqP.nexc(1) length(acqP.nexc) siraw(3)]);
    rawsl_tot=zeros([siraw(1) acqP.necho*acqP.nexc(1) length(acqP.nexc) acqP.NSlices*acqP.inrep siraw(3)]);
    %%
    for ksl=1:acqP.NSlices*acqP.inrep
        rawsl=reshape(squeeze(rawi(:,:,ksl,:,:)),[siraw(1) acqP.necho*nex siraw(3)]);
        for kPE=1:acqP.necho*nex
            rawsl_all(:,PEorder(kPE),PEcount(kPE),:)=rawsl(:,kPE,:);
            %rawsl_all(:,PEorder(kPE),1,:)=rawsl(:,kPE,:);
        end
        rawsl_tot(:,:,:,ksl,:)=reshape(rawsl_all,[siraw(1) acqP.necho*acqP.nexc(1) length(acqP.nexc) 1 siraw(3)]);
    end

    %% FT-reco
    % mp=0
    % clear imall
    if(length(acqP.wav)==1)
        rawslices=squeeze(rawsl_tot);
    else
        PEOS_reco;
    end

    zf=0;
    rawtemp=rawslices;
    sirs=size(rawslices);
    if(zf==1)

        rawzf=zeros([sirs(1) sirs(1) sirs(3) sirs(4)]);
        rawzf(:,floor((sirs(1)-sirs(2))/2)+1:floor((sirs(1)-sirs(2))/2)+sirs(2),:,:)=rawtemp;
        rawslices=rawzf;
    end
    if(length(sirs)==3), clear temp; temp(:,:,1,:)=rawslices; rawslices=temp; sirs=size(rawslices); end;
    %%
    %PEindex=floor(sort(PEorder)-min(acqP.PEorder(:))-acqP.Ny/2);
    PEindex=floor(sort(acqP.PEorder(:)))+acqP.Ny/2+1;
    if(PEindex(1)==0), PEindex=PEindex+1; end
    rawslicesz=zeros(acqP.Nx*acq.accfac,acqP.Ny,acqP.NSlices,siraw(3));
    rawslicesz(:,PEindex,:,:)=rawslices;
    if(FFTflag==1)
        for ksl=1:acqP.NSlices*acqP.inrep
            for kcoil=1:sirs(4)
                %for kcoil=[1:5 7 9:12]
                imcoil(:,:,kcoil)=squeeze(fftshift(ifft2(fftshift(squeeze(rawslicesz(:,:,ksl,kcoil))))));
            end
            %temp=sos(squeeze(imcoil(:,:,[1:5 7 9:12])));
            temp=sos(squeeze(imcoil(:,:,:)));
            imsos(:,:,ksl)=temp(end:-1:1,:)';
        end
    end
    %%


    [~,slorder]=sort(acqP.OSlices,'ascend');
    if(acqP.inrep>1)
        temp=repmat(slorder,[4 1]);
        slorder=temp(:);
    end
    imsos1=0*imsos;
    imsos1=imsos(:,:,slorder);
    if (acqP.NSlices*acqP.inrep>1)
        %bigx=im_mosaic(squeeze(imsos1),1,acqP.NSlices*acqP.inrep,200);
    else
        figure; im_tight(imsos);colormap(gray);
    end
    %try imall(:,:,1:acqP.NSlices*acqP.inrep,mp)=imsos1(:,:,:); catch warning('size of imall does not fit for kim = %i',kim); end
    imall(:,:,:,mp)=imsos1(:,:,:);

    %end
    %end
    if(nuFFTflag==1)
        myTSE_reco_nufft
        imall_n(:,:,:,mp)=imsos_nuFFT;
    end
end
%%

if(acqP.navmode(1)>0),
    rawnav1=squeeze(rawnav1);
    sipro=size(rawnav1);
    proj_1c=0*rawnav1;
    for ksl=1:acqP.NSlices*acqP.inrep
        for kproj=1:sipro(3)
            for kcoil=1:sirs(4)
                temp=squeeze(rawnav1(:,ksl,kproj,kcoil));
                proj_1c(:,ksl,kproj,kcoil)=fftshift(fft(temp));
            end
            proj_1(:,ksl,kproj)=sos(proj_1c(:,ksl,kproj,:));
        end
    end
    proj_1(:,slorder,:)=proj_1;
end

if(acqP.navmode(2)>0),
    rawnav2=squeeze(rawnav2);
    sipro=size(rawnav2);
    proj_2c=0*rawnav2;
    for ksl=1:acqP.NSlices*acqP.inrep
        for kproj=1:sipro(3)
            for kcoil=1:sirs(4)
                temp=squeeze(rawnav2(:,ksl,kproj,kcoil));
                proj_2c(:,ksl,kproj,kcoil)=fftshift(fft(temp));
            end
            proj_2(:,ksl,kproj)=sos(proj_2c(:,ksl,kproj,:));
        end
    end
    proj_2(:,slorder,:)=proj_2;
end









%%

script nufftreco
x=1

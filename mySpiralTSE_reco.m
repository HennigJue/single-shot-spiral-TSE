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
disp('select protocol(s)')
[fnp, pnp] = uigetfile('*.mat', 'Select protocol files','MultiSelect', 'on');
if(iscell(fnp)), nip=numel(fnp); else nip=1; fnp={fnp}; end
%cd(pnr)


%% Load datata
%kim=1;
% kipnp='D:\matlab\si(rawworJulk\pubiglsekimq_stassizeh\myTSE_data\sequences_20210630AA';
% pnr='D:\matlab\workb\pulseq_stash\myTSE_dsizeata\Spirals_AERA_1_7_2021';
%imall=zeros([240 240 15 16]);

%for kim=1:length(fnp),
mp=0
clear imall
for kim=1:length(fnp)
    mp=mp+1;
    cd(pnp)
    segname=cell2mat([fnp(:,1)]);
    load(segname)
    %%
    cd(pnr)
    rawname=cell2mat(fnr(:,kim))
    twix_obj = mapVBVD(rawname);
    rawdata =double(twix_obj.image.unsorted());
    %rawdata = double(twix_obj{2}.image.unsorted());
    %%
    nex=acqP.nex;    
    nslices=acqP.NSlices;
    necho=acqP.necho;
    ktraj_adc0=ktraj_adc;
    sikt0=size(ktraj_adc0);
    siktn=sikt0(2)/nex;
    %
    recnorm='y';
    spmode='cont';
    r2=0;r1=0;
    sikt0=size(ktraj_adc0);
    siktn=sikt0(2)/nex;
    ktraj_adc=ktraj_adc0(:,1:siktn);
    ktraj=ktraj_adc0(:,1:siktn./nslices);
    %
    %[ktraj_adc]=shift_ktraj_r(ktraj_adc0,[2.01,1.37,0],0);
    [ktraj_adc]=shift_ktraj_r(ktraj_adc,[0.0,0.0,0],0);
    ktraj_adc=reshape(ktraj_adc,[3,acq.nadc,necho,nslices]);
    
    %for kex=1:6,
    % Define FOV and resolution
    [nadc,nc,nacq]=size(rawdata);
    slOrder=acqP.sltemp;
    
    % stitch trajectory
    plotflag='1011';
    Nx=256;
    Ny=Nx;
       
    if(strcmp(spmode,'alt')),
        ktemp=0*ktraj_adc;
        ktemp(:,:,acqP.PEind)=ktraj_adc;
        ktemp(:,end:-1:1,2:2:end)=ktemp(:,:,2:2:end);
        ktraj_adc=ktemp;
    end
    %
    [ktacq,indkt,ksel]=spiral_stitch_nom(ktraj_adc(:,:,:,1),seg.i1,seg.nkseg,acq.nadc,necho,acqP.nTE,acq.dninc,acq.accfac);
    nTE=acqP.nTE;
     
    %
    close
     kselr=repmat(ksel,[1 1 nc]);
    
    ktraj_sp=ktacq;
    
   sikt=size(ktraj_sp);
    % here we expect Nx, Ny, deltak to be set already
    % and rawdata ktraj_adc loaded (and having the same dimensions)
    %
    os=1;
    %recnorm='n';
    % oversampling factor (we oversample both in image and k-space)
    % [acq.nadc,nc,necho]=size(rawdata);
    % nslices=floor(necho/necho);
    
    siraw=size(rawdata);
    nex=prod(siraw)/nc/length(ktraj_adc0);
    rawdata_all=reshape(rawdata,[acq.nadc nc necho nslices nex]);
    %for kex=1:1,
    for kex=1:nex,
        rawdata0=squeeze(rawdata_all(:,:,:,:,kex));
        
        if(strcmp(spmode,'alt')),
            rawtemp=0*rawdata0;
            rawtemp(:,:,acqP.PEind,:)=rawdata0;
            rawtemp(:,:,2:2:end,:)=rawtemp(end:-1:1,:,2:2:necho,:);
            rawdata0=rawtemp;
        end
        %raw=zeros([2*Nx+1 2*Ny+1 nc]);
        %im=zeros([2*Nx+1 2*Ny+1 nc+1]);
        %
        raw_sp_a=zeros([sikt(2) nc nslices]);
        Nx=spiral.Nx;Ny=Nx;
        imsos=zeros([2*Nx 2*Ny nslices]);
        imsos1=zeros([2*Nx 2*Ny nslices]);
        for nsl=1:nslices,
        %% 
            rawdata_sl=double(squeeze(rawdata0(:,:,:,nsl)));
            signorm=ones([necho 1]);
            if(strcmp(recnorm,'y'))
                
                [mph,amp]=pg_cpmg_f([1 0 0 0],1,acqP.rflip,0,acqP.TE*1000,r1,r2,0);
                signorm=squeeze(abs(mph(2,1,1:necho)));
                signorm(1)=0.8*signorm(1);
                if(strcmp(spmode,'alt'))
                    signorm(acqP.PEind)=signorm;
                end;
                disp('normalize echoes')
                for kc=1:nc,
                    for k=1:necho
                        rawdata_sl(:,kc,k)=rawdata_sl(:,kc,k)./signorm(k);
                    end
                end
            end
            
            if(strcmp(spmode,'alt')),
                rawdata_sl(:,:,2:2:necho)=rawdata_sl(end:-1:1,:,2:2:necho);
            end
            %
            
            rawdata_slp=permute(rawdata_sl,[1 3 2]);
            rawtemp=reshape(rawdata_slp,[acq.nadc*necho nc]);
            ktemp=reshape(kselr,[acq.nadc*necho nc]);
            ind=find(ktemp==1);
            raw_sp=rawtemp(ind);
            rs=size(raw_sp);
            raw_sp=reshape(raw_sp,[rs(1)/nc nc]);
            raw_sp_a(:,:,nsl)=raw_sp;
            ktraj_sp=ktacq;
            mat = [2*Nx,2*Ny];
            
            
            %compute density weightening
            ktr=-1i*(ktraj_sp(1,:))'+ktraj_sp(2,:)';
            ktr=ktr/(10*0.24/acqP.fov*Nx);
            nkp=length(ktraj_sp);
            kz=find(abs(ktr)==min(abs(ktr(:))));
            ww=abs(gradient((kt(:,1)+1i*kt(:,2))));
            ww=interp1(ww,linspace(1,length(ww),nkp))';
            w=0*ww;
            w(kz:end)= ww(1:end-kz+1);
            w(1:kz-1)=ww(end-kz+2:end);
            w=w./max(w(:))+0.01;
              %
            % % % perform SVD coil compression
            % % disp('perform coil compression at  5\% tollerance.')
            % % D = reshape(data,size(data,1)*size(data,2),size(data,3));
            % % [U,S,V] = svd(D,'econ');
            % % nCoil = max(find(diag(S)/S(1)>0.05));
            % % data = reshape(D*V(:,1:nCoil),size(data,1),size(data,2),nCoil);
            % % disp(sprintf('Using %d virtual channels',nCoil'));
            % %
            
            % Calibration
            % In this example, we calibrate from a fully sampled acquisition
            % reconstructed with gridding and density compensation.
            % Then undersample it to see how the reconstruction performs.
            
            disp('Starting Calibration');
            disp('Generating NUFFT object for calibration')
            
            tic
            % Gridding and density compensation reconstruction
            GFFT1 = NUFFT(ktr,w, [0,0] , mat);
            im = GFFT1'*(reshape(raw_sp,[length(w) 1 nc]).*repmat(sqrt(w),[1,1,nc]));
            toc
           
            imsos(:,:,nsl)=abs(sos(im));
        
            
        end
        %%
        imsos1=0*imsos;
        [ind,indr]=sort(slOrder);
        for k=1:nslices
             imsos1(:,:,k)=flipud(permute(squeeze(imsos(:,:,indr(k))),[2 1]));
        end
        imsos1=imsos1(:,:,nslices:-1:1);
        %%
         mat(3)=1;
        kcount=sprintf('%02d',kim);
        save(strcat('paramTSE',kcount),'ktraj','ksel','necho','nslices','signorm','w','slOrder','mat');
        imall(:,:,:,mp)=imsos1;
     eval(strcat('imall_',num2str(kim),'=imsos1;'));
     bigx=im_mosaic(imsos1(:,:,:,mp),2,6,150);
       

    end
end
%%
cd(pns)
seq.read(cell2mat(fns(:,kim)));
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();
fig=seq.plot();
t_ktraj(end)
t_ktraj(end)/acqP.NSlices;
%%
% rescale imsos
for k=1:nslices
    temp=imsos(:,:,k);
    temp=(temp-min(temp(:)))./(max(temp(:))-min(temp(:)));
    imsos(:,:,k)=temp;
end
im_all=im_mosaic(imsos,nslices,1,250);
%%
% resort imsos
imsos1=0*imsos;
for k=1:nslices
    temp=imsos(:,:,k);
    temp=(temp-min(temp(:)))./(max(temp(:))-min(temp(:)));
    imsos1(:,:,acqP.Oslices(k)+ceil(nslices/2))=temp;
end
im_all=im_mosaic(imsos1,1,nslices,250);

%%
raw_enc=rawdata_ch1;
% raw_enc(:,1:nsl+1:end)=0;
% raw_enc=reshape(raw_enc,[acq.nadc necho nsl]);
% raw_enc=sum(raw_enc,3);
raw_sp=raw_enc(indk);


%% single echo
raw_sp=rawdata_ch1(seg.i1:seg.i1+seg.nkseg);
ktraj_sp=ktraj_adc(1:2,seg.i1:seg.i1+seg.nkseg);



%% fixed trajectory, full adc
indk=seg.i1(1):seg.i1(1)+seg.nkseg(1);
ktraj_sp=ktraj_adc(1:2,:);
raw_sp=double(rawdata_ch1(:));


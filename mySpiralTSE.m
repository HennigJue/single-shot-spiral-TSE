%% Create a spiral TSE sequence and export for execution
%
% The |Sequence| class provides functionality to create magnetic
% resonance sequences (MRI or NMR) from basic building blocks.
%
% This provides an implementation of the open file format for MR sequences
% described here: http://pulseq.github.io/specification.pdf
%
% This example performs the following steps:
%
% # Create slice selective RF pulse for imaging.
% # Create readout gradient and phase encode strategy.
% # Loop through phase encoding and generate sequence blocks.
% # Write the sequence to an open file format suitable for execution on a
% scanner.
%
% The sequence runs in different modes defined by :
%   scanmode=
%       'init'  : initialization, creates empty k-space corrections dKA and dKE
%       'trim'  : used for experimental trimming,
%                 the first and every 2:decho:acqP.necho echo periods are
%                 readout as a spin echo
%       'run'   : run sequence
%
%   segmode=
%       'fix'   : fixed number of points per segment
%       'tan'   : tangential transition to and from each segment
%       'single': single echo spiral
%
%   spmode=
%       'cont'  : continuous segmentation
%       'alt'   : alternate between spiral-out and spiral-in in odd and
%                 even segments.
%
%   initmode=
%       'no'    : no gradient in first half of first segment
%       'rev'   : first half is reverse of second half
%       'outin' : first half 180° rotated and out-in
%       'sym'   : first half mirrored at y-axis and out-in;
%
%   accmode=    acceleration mode
%       'dd'    : dual density. First seg.n_full*acq.nadc/2 points are fully sampled,
%                 others undersampled by a factor seg.n_acc
%       'vd'    : variable density with Fcoeff2 = - Fcoeff/seg.n_acc
%     
%   fatsat
%       'on'
%       'no'    :  with fatsat gradient, but no rf
%       'off'

%% Instantiation and gradient limits
% The system gradient limits can be specified in various units _mT/m_,
% _Hz/cm_, or _Hz/m_. However the limits will be stored internally in units
% of _Hz/m_ for amplitude and _Hz/m/s_ for slew. Unspecificied hardware
% parameters will be assigned default values.
%acqPclear all
dG=200e-6;
seqname='TSE';
plotflag='00011';
initmode='no';
seq=mr.Sequence(system);


%% Sequence events
% Some sequence parameters are defined using standard MATLAB variables
setpar=1;   % if setpar =0 parameters should be read from file

if(setpar==1),
    system = mr.opts('MaxGrad', 35, 'GradUnit', 'mT/m', ...
        'MaxSlew', 170, 'SlewUnit', 'T/m/s', 'rfRingdownTime', 100e-6, ...
        'rfDeadTime', 100e-6);
    acq.slewfac=0.99;
    acq.gradfac=0.99;
    B0=2.89; 
    
    segmode='fix';
    spmode='cont';
    scanmode='run';
    fatsat='on';
    accmode='vd';
        
    if(strcmp(scanmode,'trim')), seqname=strcat(seqname,'_',segmode,'_',spmode,'_trim'); end
    %    acqP.fov=240e-3;
    spiral.Nx=120;
    spiral.kOffset=50;
    spiral.N=0.7;
    %spiral.kOffset=-spiral.Nx/acqP.fov;
    acq.dninc=2; acq.accfac=4;
    
    acqP.NSlices=15; acqP.sliceGAP=1.2; acqP.nex=1;
    seg.n_acc=1.4; 
    seg.n_full=1;
    
    
    acqP.sltemp=[1:2:acqP.NSlices 2:2:acqP.NSlices];
    acqP.Oslices=acqP.sltemp-ceil(acqP.NSlices/2);
    acqP.dummy=1;
    acqP.flipref=60; acqP.flipflag=2;
    acqP.sliceThickness=3e-3;
    acqP.TE=10e-3;
    acqP.nTE=round(70e-3/acqP.TE);
    %acqP.nTE=1;
    acqP.TR=10e-3;
    acqP.TEeff=acqP.nTE*acqP.TE;
    acqP.necho=10;
    acqP.PEtype='linear';
    
    segP.tEx=2.5e-3;
    segP.tExwd=segP.tEx+system.rfRingdownTime+system.rfDeadTime;
    segP.tRef=2e-3;
    segP.tRefwd=segP.tRef+system.rfRingdownTime+system.rfDeadTime;
    segP.tSp=1.05e-3;
    segP.tSpex=0.5*(acqP.TE-segP.tExwd-segP.tRefwd);
    segP.fspS=0.25;
    
end
%%


n_acc_s=num2str(seg.n_acc);
if(seg.n_acc==floor(seg.n_acc)),
    n_acc_s=num2str(seg.n_acc);
else
    ind=find(n_acc_s=='.'); n_acc_s(ind)='_';
end
% seqname=strcat('TSE',num2str(round(1000*acqP.TE)),'_',num2str(floor(1000*acqP.TEeff)), ...
%     '_',num2str(acqP.flipref),'_',accmode,'_',num2str(acqP.NSlices));
seqname=strcat('TSE',num2str(round(1000*acqP.TE)),'_',num2str(floor(1000*acqP.TEeff)),'_',num2str(acqP.flipref),'_',num2str(10*seg.n_acc),'_',num2str(10*spiral.N),'_',num2str(acq.accfac),'_',num2str(acqP.NSlices));

if(exist(strcat(seqname,'.seq'),'file')>0),
    m=1; testname=seqname;
    while(exist(strcat(testname,'.seq'),'file')>0)
        testname=strcat(seqname,'_',num2str(m));
        m=m+1;
    end
    seqname=testname;
end
%%

acq.readoutTime = acqP.TE-2*segP.tSp-segP.tRefwd-2*20e-6;



%%
%%% Base gradients
%%% Slice selection
% Key concepts in the sequence description are *blocks* and *events*.
% Blocks describe a group of events that are executed simultaneously. This
% hierarchical structure means that one event can be used in multiple
% blocks, a common occurrence in MR sequences, particularly in imaging
% sequences.
%
% First, the slice selective RF pulses (and corresponding slice gradient)
% are generated using the |makeSincPulse| function.
% Gradients are recalculated such that their flattime covers the pulse plus
% the rfdead- and rfringdown- times.
%
rfex_phase=pi/2; % MZ: we need to maintain these as variables because we will overwrtite phase offsets for multiple slice positions
rfref_phase=0;

acqP.flipex=90*pi/180;
[rfex, gz] = mr.makeSincPulse(acqP.flipex,system,'Duration',segP.tEx,...
    'sliceThickness',acqP.sliceThickness,'apodization',0.5,'timeBwProduct',4,'maxSlew',system.maxSlew*4);
GSex = mr.makeTrapezoid('z',system,'amplitude',gz.amplitude,'FlatTime',segP.tExwd,'riseTime',dG);
% plotPulse(rfex,GSex);

% acqP.flipref=pi;
[rfref, gz] = mr.makeSincPulse(pi,system,'Duration',segP.tRef,...
    'sliceThickness',acqP.sliceThickness,'apodization',0.5,'timeBwProduct',4,'use','refocusing','maxSlew',system.maxSlew*4);
GSref = mr.makeTrapezoid('z',system,'amplitude',GSex.amplitude,'FlatTime',segP.tRefwd,'riseTime',dG);
refenvelope=rfref.signal;
%
% Gauss pulse for refocussing
% [rfref, gz] = mr.makeGaussPulse(pi,system,'Duration',segP.tRef,...
%     'sliceThickness',acqP.sliceThickness,'apodization',0.5,'timeBwProduct',1,'use','refocusing','maxSlew',system.maxSlew*4);
% GSref = mr.makeTrapezoid('z',system,'amplitude',GSex.amplitude,'FlatTime',segP.tRefwd,'riseTime',dG);
% refenvelope=rfref.signal;
%  plotPulse(rfref,GSref);
%
AGSex=GSex.area/2;
GSspr = mr.makeTrapezoid('z',system,'area',AGSex*(1+segP.fspS),'duration',segP.tSp,'riseTime',dG);
GSspex = mr.makeTrapezoid('z',system,'area',AGSex*segP.fspS,'duration',segP.tSpex,'riseTime',dG);

% Create fat-sat pulse
% (in Siemens interpreter from January 2019 duration is limited to 8.192 ms, and although product EPI uses 10.24 ms, 8 ms seems to be sufficient)
% 1.5 2.89 3.0

acqP.sat_ppm=-3.45;
acqP.sat_freq=acqP.sat_ppm*1e-6*B0*system.gamma;
%rf_fst=1e-5*floor(1e5*min(8e-3*2.89/B0,24e-3));
rf_fs = mr.makeGaussPulse(110*pi/180,'system',system,'Duration',8e-3,...
    'bandwidth',abs(acqP.sat_freq),'freqOffset',acqP.sat_freq);
gz_fs = mr.makeTrapezoid('z',system,'delay',mr.calcDuration(rf_fs),'Area',1/1e-4); % spoil up to 0.1mm

%%
%%% Readout gradient
% To define the remaining encoding gradients we need to calculate the
% $k$-space sampling. The Fourier relationship
%
% $$\Delta k = \frac{1}{acqP.fov}$$
%
% Therefore the area of the readout gradient is $n\Delta k$.

spiral.deltak=1/acqP.fov;
spiral.kWidth = 2*spiral.Nx*spiral.deltak;

acq.ntot=round((acq.readoutTime+2*20e-6)/system.gradRasterTime);
acq.nadc=acq.accfac*round((acq.readoutTime)/system.gradRasterTime);
adc = mr.makeAdc(acq.nadc,'Duration',acq.readoutTime, 'Delay', 20e-6);%,'Delay',GRacq.riseTime);

% use spiral tools from Hargreaves to calculate spiral gradients

% convert pulseq-parameters to Hargreaves units
%%

gamma = 42.576e6;
smax = acq.slewfac*100*system.maxSlew/gamma;	 % 150 T/m/s
gmax = acq.gradfac*100*system.maxGrad/gamma;	 % G/cm

T = system.gradRasterTime;	 % Seconds

% Interleaves
%Fcoeff = [acqP.fov*100 -acqP.fov*50] ;	% acqP.fov decreases linearly from 24 to 12cm.
Fcoeff = acqP.fov*100 ;
if(strcmp(accmode,'vd')), Fcoeff=[acqP.fov*100 -acqP.fov*100/(seg.n_acc)]; end

res=100*acqP.fov/spiral.Nx;
rmax = 2*spiral.kWidth/100;		% cm^(-1), corresponds to 1mm resolution.
rmax=1/res;

disp('Calculating Gradient');

[k,g,s,time,r,theta] = vds(smax,gmax,T/acq.accfac,spiral.N,Fcoeff,rmax);
if ((spiral.N>1)&&strcmp(accmode,'dd')),
    [k1,g,s,time,r,theta] = vds(smax,gmax,T/acq.accfac,1,Fcoeff,rmax);
    k1=k1(1:floor(seg.n_full*acq.nadc/2));
    k2=k;
    rads=abs(k1(end));
    ind=find(abs(k2)>=rads);
    ang1=angle(k1(end));
    ang2=angle(k2(ind(1)));
    k2r=k2(ind(1)+1:end)*exp(i*(ang1-ang2));
    k=[k1 k2r];
end

k=k(1:acq.accfac:end);

g = 10^4/gamma*([k 0]-[0 k])/T;
g = g(1:length(k));
%plot(k)
%%
disp('Plotting Gradient');
g = [real(g(:)) imag(g(:))];
figure
plotgradinfo(g,T);
%%
if(plotflag(1)=='0'), close; end
np=length(k);
kt=zeros([np 2]);
kt(:,1)=real(k*100);
kt(:,2)=imag(k*100);
figure
plot(kt(:,1),kt(:,2),'b.');
if(plotflag(2)=='0'), close; end
Gsp=g*gamma/100;
%Gsp=Gsp.*system.maxGrad./max(Gsp(:));

acq.ninc=acq.ntot-2*acq.dninc;
spiral.kStart=floor(spiral.kWidth/2)+spiral.kOffset;
%spiral.kStart=spiral.kOffset;
if(strcmp(segmode,'single')),
    acqP.necho=1; seg.i1(1)=floor(acq.ntot/2);
    seg.ikseg=zeros([acqP.necho 2]);
    seg.ikseg(1,1)=1; seg.ikseg(1,2)=floor(acq.ninc/2);
    seg.nkseg=acq.ninc*ones([acqP.necho 1]); seg.nkseg(1)=floor(acq.ninc/2);
end


if(strcmp(segmode,'tan')),
    [seg.kseg, seg.nkseg, seg.ikseg, acqP.necho, seg.i1] = spiral_seg(kt,acq.ninc,spiral.kStart,acq.ntot,spmode);
    axlim=max(spiral.kStart,spiral.kWidth/2);
    axis([-axlim axlim -axlim axlim])
    if(plotflag(3)=='0'), close; end
end

if(strcmp(segmode,'fix')),
    
    acqP.necho=floor(np/acq.ninc);
    seg.nkseg=acq.ninc*ones([acqP.necho 1]); seg.nkseg(1)=floor(acq.ninc/2);
    seg.i1=0*seg.nkseg+acq.dninc; seg.i1(1)=floor(acq.ntot/2);
    seg.ikseg=zeros([acqP.necho 2]);
    seg.ikseg(1,1)=1; seg.ikseg(1,2)=floor(acq.ninc/2);
    
    for k=2:acqP.necho,
        seg.ikseg(k,1)=seg.ikseg(k-1,2);
        seg.ikseg(k,2)=seg.ikseg(k,1)+seg.nkseg(k)-1;
    end
    if(strcmp(spmode,'alt')),
        temp=fliplr(seg.ikseg);
        seg.ikseg(2:2:end,:)=temp(2:2:end,:);
    end
    
    if(plotflag(3)=='1'),
        figure
        axlim=max(spiral.kStart,spiral.kWidth/2);
        axis([-axlim axlim -axlim axlim])
        
        axis square
        hold on
        col=['b-';'r-'];
        if(strcmp(spmode,'cont')),
            for k=1:2:acqP.necho,
                plot(kt(seg.ikseg(k,1):seg.ikseg(k,2),1),(kt(seg.ikseg(k,1):seg.ikseg(k,2),2)),col(1,:));
                plot(kt(seg.ikseg(k,2),1),kt(seg.ikseg(k,2),2),'bo');
            end
            col=flipud(col);
            for k=2:2:acqP.necho,
                plot(kt(seg.ikseg(k,1):seg.ikseg(k,2),1),(kt(seg.ikseg(k,1):seg.ikseg(k,2),2)),col(1,:));
                plot(kt(seg.ikseg(k,1),1),kt(seg.ikseg(k,1),2),'ro');
                
            end
        end
        if(strcmp(spmode,'alt')),
            for k=1:2:acqP.necho,
                plot(kt(seg.ikseg(k,1):seg.ikseg(k,2),1),(kt(seg.ikseg(k,1):seg.ikseg(k,2),2)),col(1,:));
                plot(kt(seg.ikseg(k,2),1),kt(seg.ikseg(k,2),2),'bo');
            end
            col=flipud(col);
            for k=2:2:acqP.necho,
                plot(kt(seg.ikseg(k,1):-1:seg.ikseg(k,2),1),(kt(seg.ikseg(k,1):-1:seg.ikseg(k,2),2)),col(1,:));
                plot(kt(seg.ikseg(k,1),1),kt(seg.ikseg(k,1),2),'ro');
                
            end
        end
    end
    
    
end

%%


dKA=zeros([acqP.necho 2]);
dKE=zeros([acqP.necho 2]);


%% split gradients and recombine into blocks
% lets start with slice selection....
GS1times=[0 GSex.riseTime];
GS1amp=[0 GSex.amplitude];
GS1 = mr.makeExtendedTrapezoid('z','times',GS1times,'amplitudes',GS1amp);

GS2times=[0 GSex.flatTime];
GS2amp=[GSex.amplitude GSex.amplitude];
GS2 = mr.makeExtendedTrapezoid('z','times',GS2times,'amplitudes',GS2amp);

GS3times=[0 GSspex.riseTime GSspex.riseTime+GSspex.flatTime GSspex.riseTime+GSspex.flatTime+GSspex.fallTime];
GS3amp=[GSex.amplitude GSspex.amplitude GSspex.amplitude GSref.amplitude];
GS3 = mr.makeExtendedTrapezoid('z','times',GS3times,'amplitudes',GS3amp);

GS4times=[0 GSref.flatTime];
GS4amp=[GSref.amplitude GSref.amplitude];
GS4 = mr.makeExtendedTrapezoid('z','times',GS4times,'amplitudes',GS4amp);

GS5times=[0 GSspr.riseTime GSspr.riseTime+GSspr.flatTime GSspr.riseTime+GSspr.flatTime+GSspr.fallTime];
GS5amp=[GSref.amplitude GSspr.amplitude GSspr.amplitude 0];
GS5 = mr.makeExtendedTrapezoid('z','times',GS5times,'amplitudes',GS5amp);

GS7times=[0 GSspr.riseTime GSspr.riseTime+GSspr.flatTime GSspr.riseTime+GSspr.flatTime+GSspr.fallTime];
GS7amp=[0 GSspr.amplitude GSspr.amplitude GSref.amplitude];
GS7 = mr.makeExtendedTrapezoid('z','times',GS7times,'amplitudes',GS7amp);

% and now the readout gradient....
GRpre = mr.makeTrapezoid('x',system,'area',spiral.kStart,'duration',segP.tSpex,'riseTime',dG);
GRtrim = mr.makeTrapezoid('x',system,'area',2*spiral.kStart,'duration',acqP.TE-GSref.flatTime-2*segP.tSp,'riseTime',dG);



% and filltimes
segS.tEx=GS1.t(end)+GS2.t(end)+GS3.t(end);
segS.tRef=GS4.t(end)+GS5.t(end)+GS7.t(end)+acq.readoutTime;
tend=GS4.t(end)+GS5.t(end);
tETrain=segS.tEx+acqP.necho*segS.tRef+tend;
%TRfill=floor(1e5*(acqP.TR-acqP.NSlices*tETrain)/acqP.NSlices)*1e-5;
TRfill=acqP.TR;
if TRfill<0, TRfill=1e-3;
    disp(strcat('Warning!!! acqP.TR too short, adapted to include all slices to : ',num2str(1000*acqP.NSlices*(tETrain+TRfill)),' ms'));
else
    disp(strcat('TRfill : ',num2str(1000*TRfill),' ms'));
end
delayTR = mr.makeDelay(TRfill);
% and flip angles
acqP.rflip=acqP.flipref+zeros([1 acqP.necho]);
if(acqP.flipflag==1),  acqP.rflip(1)=90+acqP.flipref/2;end
if(acqP.flipflag==2),
    [rf,pow] = fliptraps(acqP.flipref,100,6,'opt',0,2,0,acqP.flipref,acqP.flipref,[6 5 5 acqP.necho]);
    disp(strcat('rel. power :',num2str(pow)));
    acqP.rflip=rf(1:acqP.necho);
end
if(strcmp(acqP.PEtype,'linear')), acqP.PEind=acqP.necho-(mod(acqP.nTE-1+(acqP.necho-[1:acqP.necho]),acqP.necho)); end
if(strcmp(acqP.PEtype,'centric')),
    acqP.PEind=zeros([1 acqP.necho]);
    acqP.PEind(1:acqP.nTE-1)=[2*acqP.nTE-2:-2:2];
    acqP.PEind(acqP.nTE:2*acqP.nTE-1)=[1:2:2*acqP.nTE-1];
    acqP.PEind(2*acqP.nTE:end)=[2*acqP.nTE:acqP.necho];
end

%% Define sequence blocks
% Next, the blocks are put together to form the sequence
if(~exist('scanmode')), scanmode='run', end
ntr=0;
switch scanmode
    case 'init', dKA=zeros([acqP.necho 2]); dKE=zeros([acqP.necho 2]); % dKA and dKE are trim values to correct for trajectory imperfections
    case 'trim',  decho=4; acqP.sliceGAP=0;  ntrim=[1 2:decho:acqP.necho]; acqP.NSlices=length(ntrim); acqP.dummy=1;
    case 'allTE',  acqP.sliceGAP=0; acqP.NSlices=acqP.necho;
end
nk=acqP.necho;

for kex=acqP.dummy:acqP.nex
    for s=1:acqP.NSlices,
        rfex.freqOffset=acqP.sliceGAP*GSex.amplitude*acqP.sliceThickness*(acqP.Oslices(s));
        rfref.freqOffset=acqP.sliceGAP*GSref.amplitude*acqP.sliceThickness*(acqP.Oslices(s));
        rfex.phaseOffset=rfex_phase-2*pi*rfex.freqOffset*mr.calcRfCenter(rfex); % align the phase for off-center slices
        rfref.phaseOffset=rfref_phase-2*pi*rfref.freqOffset*mr.calcRfCenter(rfref); % dito
        
        
        
        if(strcmp(fatsat,'on')), seq.addBlock(rf_fs,gz_fs); end
        if(strcmp(fatsat,'no')), seq.addBlock(gz_fs); end
        seq.addBlock(GS1);
        seq.addBlock(GS2,rfex);
        seq.addBlock(GS3,GRpre);
        
        if(ntr>0), nk=s; end
        %for m=1:1,
        for m=1:acqP.necho,
            
            if(strcmp(scanmode,'allTE')),acqP.nTE=s;end
            %k=nk-(mod(acqP.nTE-1+(nk-m),nk))
            k=acqP.PEind(m);
            % calculate gradients
            if((k==1)&&strcmp(initmode,'rev')),
                % initial gradient in central echo with reversal of the
                % trajectory in first half
                ind=find(kt(1:seg.nkseg(1),2)>0);
                nkseg1=max(ind);
                segS.tSp=GSspr.riseTime+GSspr.flatTime+GSspr.fallTime+(seg.i1(k)-nkseg1)*system.gradRasterTime;
                kBegin1=[-spiral.kStart+dKA(k,1) 0];
                kEnd1=[kt(nkseg1,1) kt(nkseg1,2)+dKA(k,2)];
                GStart=[0 0];
                GDest=[-Gsp(seg.ikseg(k,2),1) -Gsp(seg.ikseg(k,2),2)];
                if(strcmp(segmode,'tan'))
                    [GBegin, tBegin, Ginit, tinit] = spiral_k2k_m(kBegin1, kEnd1, GStart, GDest, segS.tSp, system.maxSlew, system.gradRasterTime);
                else
                    [GBegin, tBegin, tt] = seqspiral_k2k_opt(kBegin1, kEnd1, GStart, GDest, system.maxGrad, system.maxSlew, segS.tSp, system.gradRasterTime);
                end
            else
                %% initial gradient
                segS.tSp=GSspr.riseTime+GSspr.flatTime+GSspr.fallTime+seg.i1(k)*system.gradRasterTime;
                kBegin1=[-spiral.kStart+dKA(k,1) 0];
                kEnd1=[kt(seg.ikseg(k,1),1) kt(seg.ikseg(k,1),2)+dKA(k,2)];
                GStart=[0 0];
                if(strcmp(spmode,'cont')|(k~=2*floor(k/2))),
                    GDest=[Gsp(seg.ikseg(k,1),1) Gsp(seg.ikseg(k,1),2)];
                end
                if(strcmp(spmode,'alt')&&(k==2*floor(k/2))),
                    GDest=-[Gsp(seg.ikseg(k,1),1) Gsp(seg.ikseg(k,1),2)];
                end
                if(strcmp(segmode,'tan')),
                    [GBegin, tBegin, Ginit, tinit] = spiral_k2k_m(kBegin1, kEnd1, GStart, GDest, segS.tSp, system.maxSlew, system.gradRasterTime);
                else
                    [GBegin, tBegin,Gtran, tG, tt] = spiral_k2k_opt(kBegin1, kEnd1, GStart, GDest, system.maxGrad, system.maxSlew, segS.tSp, system.gradRasterTime);
                end
                
            end
            
            % terminal gradient
            ntk=acq.ntot-seg.i1(k)-seg.nkseg(k);
            segS.tSp=GSspr.riseTime+GSspr.flatTime+GSspr.fallTime+ntk*system.gradRasterTime;
            kBegin2=[kt(seg.ikseg(k,2),1) kt(seg.ikseg(k,2),2)+dKE(k,2)];
            kEnd2=[spiral.kStart+dKE(k,1) 0];
            GDest=[0 0];
            if(strcmp(spmode,'cont')|(k~=2*floor(k/2))),
                GStart=[Gsp(seg.ikseg(k,2),1) Gsp(seg.ikseg(k,2),2)];
            end
            if(strcmp(spmode,'alt')&&(k==2*floor(k/2)))
                GStart=-[Gsp(seg.ikseg(k,2),1) Gsp(seg.ikseg(k,2),2)];
            end
            if(strcmp(segmode,'tan')),
                [GEnd, tEnd, Gfin, tfin] = spiral_k2k_m(kBegin2, kEnd2, GStart, GDest, segS.tSp, system.maxSlew, system.gradRasterTime);
            else
                [GEnd, tEnd, Gfin, tfin] = spiral_k2k_opt(kBegin2, kEnd2, GStart, GDest, system.maxGrad, system.maxSlew, segS.tSp, system.gradRasterTime);
                
            end
            % combine gradients
            
            
            if((k==1)&&strcmp(initmode,'rev')),
                Gtot=[GBegin(1:end,:); -Gsp(nkseg1:-1:seg.ikseg(k,1),:); Gsp(seg.ikseg(k,1)+1:seg.ikseg(k,2),:); GEnd(1:end,:);[0 0]];
            else
                if((k==1)&&strcmp(initmode,'outin')),
                    Gtot(:,1)=[GEnd(end:-1:1,1); Gsp(seg.ikseg(k,2):-1:seg.ikseg(k,1),1); Gsp(seg.ikseg(k,1)+1:seg.ikseg(k,2),1); GEnd(1:end,1);[0]];
                    Gtot(:,2)=[GEnd(end:-1:1,2); Gsp(seg.ikseg(k,2):-1:seg.ikseg(k,1),2); Gsp(seg.ikseg(k,1)+1:seg.ikseg(k,2),2); GEnd(1:end,2);[0]];
                else
                    
                    if((k==1)&&strcmp(initmode,'sym')),
                        Gtot(:,1)=[GEnd(end:-1:1,1); Gsp(seg.ikseg(k,2):-1:seg.ikseg(k,1),1); Gsp(seg.ikseg(k,1)+1:seg.ikseg(k,2),1); GEnd(1:end,1);[0]];
                        Gtot(:,2)=[-GEnd(:,2); -Gsp(seg.ikseg(k,2):-1:seg.ikseg(k,1),2); Gsp(seg.ikseg(k,1)+1:seg.ikseg(k,2),2); GEnd(1:end,2);[0]];
                    else
                        
                        if(strcmp(spmode,'cont')|(k~=2*floor(k/2))),
                            Gtot=[GBegin(1:end,:); Gsp(seg.ikseg(k,1)+1:seg.ikseg(k,2),:); GEnd(1:end,:);[0 0]];
                        end
                        if(strcmp(spmode,'alt')&&(k==2*floor(k/2))),
                            Gtot=[GBegin(1:end,:); -Gsp(seg.ikseg(k,1):-1:seg.ikseg(k,2)+1,:); GEnd(1:end,:);[0 0]];
                        end
                        
                        ttot=[tBegin(1:end) [(1:seg.nkseg(k))*system.gradRasterTime] tEnd(2:end)];
                        
                    end
                end
            end
            
            nGtot=length(Gtot);
            ttot=[1:nGtot]*system.gradRasterTime;
            
            
            
            
            % split Gtot
            nsp=round((GSspr.riseTime+GSspr.flatTime+GSspr.fallTime)/system.gradRasterTime);
            
            Gsp1=Gtot(1:nsp,:);
            Gspiral=Gtot(nsp+1:nsp+acq.ntot,:);
            Gsp2=Gtot(nsp+acq.ntot+1:end,:);
            
            Gsp1_x=mr.makeArbitraryGrad('x','system',system,'waveform',Gsp1(:,1));
            Gsp1_y=mr.makeArbitraryGrad('y','system',system,'waveform',Gsp1(:,2));
            
            if(strcmp(segmode,'single')),
                for kadc=1:1,
                    Gspiral_x(kadc)=mr.makeArbitraryGrad('x','system',system,'maxSlew',4.2576e+11,'waveform',Gspiral(:,1));
                    Gspiral_y(kadc)=mr.makeArbitraryGrad('y','system',system,'maxSlew',4.2576e+11,'waveform',Gspiral(:,2));
                end
                
            else
                Gspiral_x=mr.makeArbitraryGrad('x','system',system,'maxSlew',4.2576e+11,'waveform',Gspiral(:,1));
                Gspiral_y=mr.makeArbitraryGrad('y','system',system,'maxSlew',4.2576e+11,'waveform',Gspiral(:,2));
            end
            
            
            Gsp2_x=mr.makeArbitraryGrad('x','system',system,'waveform',Gsp2(:,1));
            Gsp2_y=mr.makeArbitraryGrad('y','system',system,'waveform',Gsp2(:,2));
            rfref.signal=refenvelope*acqP.rflip(m)/180;
            seq.addBlock(GS4,rfref);
            
            if(strcmp(scanmode,'trim')&&(k==ntrim(s))),
                
                seq.addBlock(GS5);
                seq.addBlock(GRtrim,adc);
                seq.addBlock(GS7);
                
            else
                if(length(Gsp1_x.t)>length(GS5.t))
                    fprintf('\ntrajectory not compatible with sequence parameters. \nChange parameters and/or buy faster/stronger gradients\n');
                    return
                end
                seq.addBlock(Gsp1_x,Gsp1_y,GS5);
                if (kex>0),
                    seq.addBlock(Gspiral_x,Gspiral_y,adc);
                else
                    seq.addBlock(Gspiral_x,Gspiral_y);
                end
                if(length(Gsp2_x.t)>length(GS7.t))
                    fprintf('\ntrajectory not compatible with sequence parameters. \nChange parameters and/or buy faster/stronger gradients\n');
                    return
                end
                seq.addBlock(Gsp2_x,Gsp2_y,GS7);
            end
            
        end
        seq.addBlock(GS4);
        seq.addBlock(GS5);
        seq.addBlock(delayTR);
        
    end
    
    
    
end




%% new single-function call for trajectory calculation
[ktraj_adc, ktraj, t_excitation, t_refocusing] = seq.calculateKspace();

%% plot k-spaces

% figure; plot(ktraj'); % plot the entire k-space trajectory
if(plotflag(4)=='1')
    figure; plot(ktraj(1,:),ktraj(2,:),'b',...
        ktraj_adc(1,:),ktraj_adc(2,:),'r'); % a 2D plot
    set(gcf,'Position',[100 100 550 350]);
    axis equal
end

%% Write to file

% The sequence is written to file in compressed form according to the file
% format specification using the |write| method.
%pn=which('myspiralTSE');
seq.write(strcat(seqname,'.seq'))
save(strcat('p',seqname),'segmode','spmode','scanmode','accmode','fatsat','acq','acqP','seg','segP','spiral','kt','ktraj_adc','system','B0');
%mySpiralTSE_reco_prep;
%save(strcat('./seqs/p',seqname),'segmode','spmode','scanmode','accmode','fatsat','acq','acqP','seg','spiral','kt','ktraj_adc','system');
%%
% Display the first few lines of the output file
s=fileread(strcat(seqname,'.seq'));
disp(s(1:300))
%% plot gradients etc
if(plotflag(5)=='1'),
    fig=seq.plot('TimeRange',[0 (acqP.necho+3)*acqP.TE+acq.readoutTime],'timeDisp','ms');
end

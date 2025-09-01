%   sequence set for 7T (B0 and TimeBandwidth Product of RF-pulses)
%
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
%       'trim'  : used for experimental trimming,
%                 the first and every 2:decho:acqP.necho echo periods are
%                 readout as a spin echo
%       'run'   : run sequence
%
%
%   accmode=    acceleration mode
%       'dd'    : dual density. First seg.n_full*acq.nadc/2 points are fully sampled,
%                 others undersampled by a factor spiral.Nout
%       'vd'    : variable density with Fcoeff2 = - Fcoeff/spiral.Nout
%       (other) : undersampled by a factor spiral.Nout
%
%   fatsat      : 'on' or 'off'
%   spiralmode  : 'fly_through' fast passage through center


%% Instantiation and gradient limits
% The system gradient limits can be specified in various units _mT/m_,
% _Hz/cm_, or _Hz/m_. However the limits will be stored internally in units
% of _Hz/m_ for amplitude and _Hz/m/s_ for slew. Unspecificied hardware
% parameters will be assigned default values.
%clear all
dG=150e-6;
system = mr.opts('MaxGrad', 35, 'GradUnit', 'mT/m', ...
    'MaxSlew', 150, 'SlewUnit', 'T/m/s', 'rfRingdownTime', 100e-6, ...
    'rfDeadTime', 100e-6);
B0=7;
plotflag='00011';

seqname='SoSTSE';
spiralmode='stop';
gradmode='all';

%%
% A new sequence object is created by calling the class constructor.
seq=mr.Sequence(system);
%warning('OFF', 'mr:restoreShape')
warning('OFF', 'all')

%% Set parameters
setpar=1;   % if setpar =0 parameters should be read from file

if(setpar==1),
% mode
spiralmode='stop';
scanmode='3D';
fatsat='on';
accmode='vd';

% acq
acq.dninc=2; acq.accfac=4;

% spiral
spiral.Nx=30;
spiral.Nout=1.2;
spiral.kOffset=120;

% acqP
acqP.fov=200e-3;
acqP.NSlices=1; acqP.sliceGAP=1.5; acqP.nex=1;
acqP.nPartitions=1;
acqP.nPartEx=1;
acqP.sltemp=ceil(acqP.NSlices/2):2:acqP.NSlices;acqP.sltemp=[acqP.sltemp 2:2:floor(acqP.NSlices/2)];
acqP.sltemp=[acqP.sltemp ceil(acqP.NSlices/2)+1:2:acqP.NSlices];acqP.sltemp=[acqP.sltemp 1:2:floor(acqP.NSlices/2)];
acqP.Oslices=acqP.sltemp-ceil(acqP.NSlices/2);
acqP.dummy=1;   % = 0 for dummy scan, otherwise = 1
acqP.flipref=60; acqP.flipflag=2;
acqP.sliceThickness=160e-3;
acqP.oversampling3D=1.10;
acqP.fov3D=acqP.sliceThickness*acqP.oversampling3D;
%acqP.TE=10e-3; 
acqP.nTE=1; acqP.TR=5000e-3;
%acqP.TEeff=acqP.nTE*acqP.TE;
acqP.necho=21;
acqP.PEtype='centric';
acqP.PEmode='quadratic';
acqP.PEfac=1.4;

% seg
seg.n_full=1;

% segP
segP.tEx=2.5e-3;
segP.tExwd=segP.tEx+system.rfRingdownTime+system.rfDeadTime;
segP.tRef=2e-3;
segP.tRefwd=segP.tRef+system.rfRingdownTime+system.rfDeadTime;
segP.tSp=1e-3;
%segP.tSpex=0.5*(acqP.TE-segP.tExwd-segP.tRefwd);
if(strcmp(scanmode,'3D')),
    segP.fspS=-1;
else
    segP.fspS=0.0;
end
end

%% Sequence events
% Some sequence parameters are defined using standard MATLAB variables

decho=4;    %only used for Trim 
delta0=0;

% slice order for multislice (acqP.NSlices should be odd)

if(strcmp(scanmode,'3D')), acqP.nPartitions=acqP.necho.*acqP.nPartEx; end

n_acc_s=num2str(spiral.Nout);
if(spiral.Nout==floor(spiral.Nout)),
    n_acc_s=num2str(spiral.Nout);
else
    ind=find(n_acc_s=='.'); n_acc_s(ind)='_';
end

%%
spiral.deltak=1/acqP.fov;
spiral.kWidth = 2*spiral.Nx*spiral.deltak;




%%

gamma = 42.576e6;
smax = 100*system.maxSlew/gamma;	 % 150 T/m/s
gmax = 100*system.maxGrad/gamma;	 % G/cm

T = system.gradRasterTime;	 % Seconds
%N = spiral.Nout;		 % Interleaves
%Fcoeff = [acqP.fov*100 -acqP.fov*50] ;	% FOV decreases linearly from 24 to 12cm.
Fcoeff = acqP.fov*100 ;
if(strcmp(accmode,'vd')), Fcoeff=[acqP.fov*100 -acqP.fov*100/spiral.Nout]; end
%Fcoeff=[acqP.fov*100 -acqP.fov*100/(spiral.Nout)];
res=100*acqP.fov/spiral.Nx;
rmax = 2*spiral.kWidth/100;		% cm^(-1), corresponds to 1mm resolution.
rmax=1/res;

disp('Calculating Gradient');
%%ac
[k_out,g,s,time,r,theta] = vds(smax,gmax,T/acq.accfac,2,Fcoeff,rmax);
k_out=k_out(1:acq.accfac:end);
acqP.n_int=1;
% k_out(1)=0;
%k_out=k_out(1:acq.accfac:acq.accfac*(round(acq.ninc/2)+1));
k_in=-k_out(end:-1:1);
k_int=[k_in k_out];
%k_int=[k_in(2:end-1) 0*k_in(end) 0*k_in(end) k_out(2:end-1)];
%%
acq.ntot=length(k_int);
acq.ninc=acq.ntot-2*acq.dninc;
acq.readoutTime = length(k_int)*system.gradRasterTime-2*20e-6;
acqP.TE = acq.readoutTime +2*segP.tSp+segP.tRefwd+2*20e-6;
acqP.TEeff=acqP.nTE*acqP.TE;
segP.tSpex=0.5*(acqP.TE-segP.tExwd-segP.tRefwd);
acq.nadc=acq.accfac*round((acq.readoutTime)/system.gradRasterTime);
adc = mr.makeAdc(acq.nadc,'Duration',acq.readoutTime, 'Delay', segP.tSp+20e-6);


%%
%%
spiral.kStart=floor(spiral.kWidth/2)+spiral.kOffset;
ang1=acos(100*abs(k_int(1))/abs(spiral.kStart));
ang2=pi+angle(k_int(1));
k_rot=k_int*exp(i*(-ang1-ang2));

%
g = 10^4/gamma*([k_rot 0]-[0 k_rot])/T;
g = g(1:length(k_rot));
g(1)=g(2)-(g(3)-g(2));

disp('Plotting Gradient');
g = [real(g(:)) imag(g(:))];
figure
plotgradinfo(g,T);

if(plotflag(1)=='0'), close; end
np=length(k_rot);
kt=zeros([np 2]);
kt(:,1)=real(k_rot*100);
kt(:,2)=imag(k_rot*100);
figure
plot(kt(:,1),kt(:,2),'b-');
axis equal
%%
if(plotflag(2)=='0'), close; end
Gsp=g*gamma/100;
spiral.kStart=floor(spiral.kWidth/2)+spiral.kOffset;
%spiral.kStart=spiral.kOffset;
seg.i1=acq.dninc;seg.nkseg=acq.ninc;seg.ikseg=[1 acq.ninc];





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

%%
acq.flipex=90*pi/180;
[rfex, gz] = mr.makeSincPulse(acq.flipex,system,'Duration',segP.tEx,'use','excitation',...
    'sliceThickness',acqP.sliceThickness,'apodization',0.5,'timeBwProduct',2,'maxSlew',system.maxSlew*4,'PhaseOffset',rfex_phase);
GSex = mr.makeTrapezoid('z',system,'amplitude',gz.amplitude,'FlatTime',segP.tExwd,'riseTime',dG);
 %plotPulse(rfex,GSex);
%%

[rfref, gz] = mr.makeSincPulse(pi,system,'Duration',segP.tRef,...
    'sliceThickness',acqP.sliceThickness,'apodization',0.5,'timeBwProduct',2,'PhaseOffset',rfref_phase,'use','refocusing','maxSlew',system.maxSlew*4);
GSref = mr.makeTrapezoid('z',system,'amplitude',GSex.amplitude,'FlatTime',segP.tRefwd,'riseTime',dG);
refenvelope=rfref.signal;
% plotPulse(rfref,GSref);

AGSex=GSex.area/2;
GSspr = mr.makeTrapezoid('z',system,'area',AGSex*(1+segP.fspS),'duration',segP.tSp,'riseTime',dG);
GSspex = mr.makeTrapezoid('z',system,'area',AGSex*segP.fspS,'duration',segP.tSpex,'riseTime',dG);

% Create fat-sat pulse
% (in Siemens interpreter from January 2019 duration is limited to 8.192 ms, and although product EPI uses 10.24 ms, 8 ms seems to be sufficient)
acqP.sat_ppm=-3.45;
acqP.sat_freq=acqP.sat_ppm*1e-6*B0*system.gamma;
rf_fst=8e-3;
if (B0<2), rf_fst=1e-5*floor(1e5*10e-3*1.5/B0); end
rf_fs = mr.makeGaussPulse(180*pi/180,'system',system,'Duration',rf_fst,'use','other',...
    'bandwidth',abs(acqP.sat_freq),'freqOffset',acqP.sat_freq);
gz_fs = mr.makeTrapezoid('z',system,'delay',mr.calcDuration(rf_fs),'Area',1/1e-4); % spoil up to 0.1mm

%%
%%% Readout gradient
% To define the remaining encoding gradients we need to calculate the
% $k$-space sampling. The Fourier relationship
%
% $$\Delta k = \frac{1}{FOV}$$
%
% Therefore the area of the readout gradient is $n\Delta k$.

spiral.deltak=1/acqP.fov;
spiral.kWidth = 2*spiral.Nx*spiral.deltak;

acq.ntot=round((acq.readoutTime+2*20e-6)/system.gradRasterTime);
acq.ninc=acq.ntot-2*acq.dninc;
acq.nadc=acq.accfac*round((acq.readoutTime)/system.gradRasterTime);
%acq.readoutTime=7.8e-3;
adc = mr.makeAdc(acq.nadc,'Duration',acq.readoutTime, 'Delay', segP.tSp+20e-6);%,'Delay',GRacq.riseTime);
%adc = mr.makeAdc(acq.nadc,'Duration',acq.readoutTime, 'Delay', 20e-6);%,'Delay',GRacq.riseTime);
% use spiral tools from Hargreaves to calculate spiral gradients

% convert pulseq-parameters to Hargreaves units
%% calculate spiral gradient

% initial gradient
%segP.tSp=GSspr.riseTime+GSspr.flatTime+GSspr.fallTime+seg.i1*system.gradRasterTime;
segP.tSp=GSspr.riseTime+GSspr.flatTime+GSspr.fallTime;
kBegin1=[-spiral.kStart 0];
%kEnd1=[kt(1,1) kt(1,2)];
kEnd1=[-Gsp(1,1)*system.gradRasterTime+kt(1,1) -Gsp(1,2)*system.gradRasterTime+kt(1,2)];

GStart1=[0 0];
GDest1=[Gsp(1,1) Gsp(1,2)];
k=1;
[GBegin, tBegin, Gfin, tfin] = spiral_k2k_opt(kBegin1, kEnd1, GStart1, GDest1, system.maxGrad, system.maxSlew, segP.tSp, system.gradRasterTime);

%%
% terminal gradient
%ntk=acq.ntot-seg.i1-seg.nkseg;
%segP.tSp=GSspr.riseTime+GSspr.flatTime+GSspr.fallTime+ntk*system.gradRasterTime;
kBegin2=[kt(end,1) kt(end,2)];
kEnd2=[spiral.kStart 0];
GStart2=[Gsp(end,1) Gsp(end,2)];
GDest2=[0 0];
[GEnd, tEnd, Gfin, tfin] = spiral_k2k_opt(kBegin2, kEnd2, GStart2, GDest2, system.maxGrad, system.maxSlew, segP.tSp, system.gradRasterTime);

% combine gradients

Gtot=[GBegin; Gsp; GEnd];
%%

nGtot=length(Gtot);
ttot=[1:nGtot]*system.gradRasterTime;

% split Gtot
nsp=(GSspr.riseTime+GSspr.flatTime+GSspr.fallTime)/system.gradRasterTime;

Gsp1=Gtot(1:nsp,:);
Gspiral=Gtot(nsp+1:nsp+acq.ntot,:);
Gsp2=Gtot(nsp+acq.ntot+1:end,:);

%          plot((Gsp2(:,1)));
Gsp1_x=mr.makeArbitraryGrad('x',system,'waveform',Gsp1(:,1),'first',Gsp1(1,1),'last',Gsp1(end,1));
Gsp1_y=mr.makeArbitraryGrad('y',system,'waveform',Gsp1(:,2),'first',Gsp1(1,2),'last',Gsp1(end,2));

Gspiral_x=mr.makeArbitraryGrad('x',system,'maxSlew',4.2576e+11,'waveform',Gspiral(:,1),'first',Gspiral(1,1),'last',Gspiral(end,1));
Gspiral_y=mr.makeArbitraryGrad('y',system,'maxSlew',4.2576e+11,'waveform',Gspiral(:,2),'first',Gspiral(1,2),'last',Gspiral(end,2));

Gsp2_x=mr.makeArbitraryGrad('x',system,'waveform',Gsp2(:,1),'first',Gsp2(1,1),'last',Gsp2(end,1));
Gsp2_y=mr.makeArbitraryGrad('y',system,'waveform',Gsp2(:,2),'first',Gsp2(1,2),'last',Gsp2(end,2));

%%
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
%%
GS_alltimes=[GS5times GS5times(end)+acq.ntot*system.gradRasterTime GS5times(end)+acq.ntot*system.gradRasterTime+GS7times(2:end)];
GS_allamp=[GS5amp GS7amp];
GS_allamp0=GS_allamp;
GS_all = mr.makeExtendedTrapezoid('z','times',GS_alltimes,'amplitudes',GS_allamp);
%%

% and now the readout gradient....
GRpre = mr.makeTrapezoid('x',system,'area',spiral.kStart,'duration',segP.tSpex,'riseTime',dG);
GRtrim = mr.makeTrapezoid('x',system,'area',2*spiral.kStart,'duration',acqP.TE-GSref.flatTime-2*segP.tSp,'riseTime',dG);
%%
%%% Phase encoding
% To move the $k$-space trajectory away from 0 prior to the readout a
% prephasing gradient must be used. Furthermore rephasing of the slice
% select gradient is required.
if(strcmp(scanmode,'3D')),
    [PEorder,nPE] = myTSE_PEorder(acqP.necho,acqP.necho,acqP.nTE,acqP.PEtype)
    
%     nex=1;
%     if((acqP.nPartitions/2-floor(acqP.nPartitions/2))==0),
%     pe_steps=(1:acqP.nPartitions)-0.5*acqP.nPartitions-1;
%     else
%     pe_steps=(1:acqP.nPartitions)-0.5*(acqP.nPartitions+1);  
%     end
%     
%     PEorder=reshape(pe_steps,[nex,acqP.necho])';
    phaseAreas = PEorder/acqP.fov3D;
    if(strcmp(acqP.PEmode,'quadratic')),
        temp=PEorder;
        phaseAreas=(sign(temp).*abs(temp).^(acqP.PEfac))/acqP.fov3D;
    end
end
%%
% and filltimes
segP.tExf=GS1.tt(end)+GS2.tt(end)+GS3.tt(end);
segP.tReff=GS4.tt(end)+GS5.tt(end)+GS7.tt(end)+acq.readoutTime;
tend=GS4.tt(end)+GS5.tt(end);
tETrain=segP.tExf+acqP.necho*segP.tReff+tend;
TRfill=floor(1e5*(acqP.TR-acqP.NSlices*tETrain)/acqP.NSlices)*1e-5;
if TRfill<0, TRfill=1e-3;
    disp(strcat('Warning!!! acqP.TR too short, adapted to include all slices to : ',num2str(1000*acqP.NSlices*(tETrain+TRfill)),' ms'));
else
    disp(strcat('TRfill : ',num2str(1000*TRfill),' ms'));
end
delayTR = mr.makeDelay(TRfill);
% and flip angles
acqP.rflip=acqP.flipref+zeros([1 acqP.necho]);
if(acqP.flipflag==1),  acqP.rflip(1)=90+acqP.flipref(1)/2;end
if(acqP.flipflag==2),
    [rf,pow] = fliptraps(acqP.flipref(1),100,6,'opt',0,2,0,acqP.flipref(1),acqP.flipref(1),[6 5 5 acqP.necho]);
    disp(strcat('rel. power :',num2str(pow)));
    acqP.rflip=rf(1:acqP.necho);
end

%% Define sequence blocks
% Next, the blocks are put together to form the sequence
if(~exist('scanmode')), scanmode='2D', end
ntr=0;
switch scanmode
    
    case 'trim',  acqP.sliceGAP=0;  ntrim=[1 2:decho:acqP.necho]; acqP.NSlices=length(ntrim);
    case 'allTE',  acqP.sliceGAP=0; acqP.NSlices=1;
end
nk=acqP.necho;
%acqP.n_int=2*64;
nseg=acqP.n_int;


delta=pi/nseg;
if(strcmp(gradmode,'xy')), nseg=8; delta=2*pi/nseg; end
    
%nseg=2*nseg;

testg=zeros([nseg 7]);
for kex=acqP.dummy:acqP.nex,
    for kseg=1:nseg,
        %for kseg=3:3,
        if(strcmp(fatsat,'on')), seq.addBlock(rf_fs,gz_fs); end
        
        for s=1:acqP.NSlices,
            %         rfex.freqOffset=acqP.sliceGAP*GSex.amplitude*acqP.sliceThickness*(s-1-(acqP.NSlices-1)/2);
            %         rfref.freqOffset=acqP.sliceGAP*GSref.amplitude*acqP.sliceThickness*(s-1-(acqP.NSlices-1)/2);
            rfex.freqOffset=acqP.sliceGAP*GSex.amplitude*acqP.sliceThickness*acqP.Oslices(s);
            rfref.freqOffset=acqP.sliceGAP*GSref.amplitude*acqP.sliceThickness*acqP.Oslices(s);
            rfex.phaseOffset=rfex_phase-2*pi*rfex.freqOffset*mr.calcRfCenter(rfex); % align the phase for off-center slices
            rfref.phaseOffset=rfref_phase-2*pi*rfref.freqOffset*mr.calcRfCenter(rfref); % dito
    
            
            
             phi=delta0+delta*(kseg-1);
            seq.addBlock(GS1);
            seq.addBlock(GS2,rfex);
            area_x=spiral.kStart*cos(phi);
            area_y=spiral.kStart*sin(phi);
            
            GRpre_x = mr.makeTrapezoid('x',system,'area',area_x,'duration',segP.tSpex,'riseTime',dG);
            GRpre_y = mr.makeTrapezoid('y',system,'area',area_y,'duration',segP.tSpex,'riseTime',dG);
            testg(kseg,1)=area_x;testg(kseg,2)=area_y;testg(kseg,3)=abs(area_x+i*area_y);
            testg(kseg,4)=GRpre_x.amplitude;testg(kseg,5)=GRpre_y.amplitude;testg(kseg,6)=abs(testg(kseg,4)+i*testg(kseg,5));
            testg(kseg,7)=phi*180/pi;
            
            %seq.addBlock(mr.rotate('z',phi,GS3,GRpre));
            
            
            
            
            if(strcmp(gradmode,'y'));
                seq.addBlock(GS3);
            else
            seq.addBlock(mr.rotate('z',phi,GS3,GRpre));
            end
           
            if(ntr>0), nk=s; end
            
            for m=1:acqP.necho,
                %for m=1:1,
                if(strcmp(scanmode,'allTE')),acqP.nTE=s;end
                k=nk-(mod(acqP.nTE-1+(nk-m),nk));
                
                
                rfref.signal=refenvelope*acqP.rflip(m)/180;
           
                seq.addBlock(GS4,rfref);
                
                %             seq.addBlock(mr.rotate('z',phi,Gsp1_x,Gsp1_y,GS5));
                %             seq.addBlock(mr.rotate('z',phi,Gspiral_x,Gspiral_y,adc));
                %             seq.addBlock(mr.rotate('z',phi,Gsp2_x,Gsp2_y,GS7));
                if(strcmp(scanmode,'3D')),
                    if (kex>0)
                        phaseArea=phaseAreas(m,kex);
                    else
                        phaseArea=0;
                    end
                    GPpre = mr.makeTrapezoid('y',system,'Area',phaseArea,'Duration',segP.tSp,'riseTime',dG);
                    GPrew = mr.makeTrapezoid('y',system,'Area',-phaseArea,'Duration',segP.tSp,'riseTime',dG);
                    GS_alltimes=[GS5times GS5times(end)+acq.ntot*system.gradRasterTime GS5times(end)+acq.ntot*system.gradRasterTime+GS7times(2:end)];
                    GS_allamp(2:3)=GS_allamp0(2:3)+GPpre.amplitude;
                    GS_allamp(6:7)=GS_allamp0(6:7)+GPrew.amplitude;
                    GS_all = mr.makeExtendedTrapezoid('z','times',GS_alltimes,'amplitudes',GS_allamp);
                end
                Gtemp_x=mr.makeArbitraryGrad('x',system,'maxSlew',4.2576e+11,'waveform',Gtot(:,1),'first',Gtot(1,1),'last',Gtot(end,1));
                Gtemp_y=mr.makeArbitraryGrad('y',system,'maxSlew',4.2576e+11,'waveform',Gtot(:,2),'first',Gtot(1,2),'last',Gtot(end,2));
                
                if(strcmp(gradmode(1),'x')),Gtemp_y.waveform=0*Gtemp_y.waveform; end
                if(strcmp(gradmode,'y')),Gtemp_x.waveform=0*Gtemp_x.waveform; end
                               
                if (kex>0),
                    seq.addBlock(mr.rotate('z',phi,Gtemp_x,Gtemp_y,GS_all,adc));
                else
                    seq.addBlock(mr.rotate('z',phi,Gtemp_x,Gtemp_y,GS_all));
                end
                
                %             seq.addBlock(Gsp1_x,Gsp1_y,GS5);
                %             seq.addBlock(Gspiral_x,Gspiral_y,adc);
                %             seq.addBlock(Gsp2_x,Gsp2_y,GS7);
                %
            end
            seq.addBlock(GS4);
            seq.addBlock(GS5);
            seq.addBlock(delayTR);
        end
    end
end




%% new single-function call for trajectory calculation
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();

%% plot k-spaces

% figure; plot(ktraj'); % plot the entire k-space trajectory
if(plotflag(4)=='1'),
    figure; plot(ktraj(1,:),ktraj(2,:),'b',...
        ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % a 2D plot
    set(gcf,'Position',[100 100 550 350]);
    axis equal
    XL=xlim;
    YL=ylim;
    figure; plot(ktraj(1,:),ktraj(3,:),'b',...
        ktraj_adc(1,:),ktraj_adc(3,:),'r.'); % a 2D plot
    set(gcf,'Position',[700 100 550 350]);
    axis([XL YL]);
    %axis equal
end
%% 
%seqname='SoS_TSE';
seqname=strcat(seqname,'_',scanmode,'_',accmode,'_',num2str(round(1000*acqP.TE)),'_',num2str(round(1000*acqP.TEeff)), ...
    '_',num2str(acqP.necho),'_',num2str(1000*acqP.sliceThickness));
if(exist(strcat(seqname,'.seq'),'file')>0),
    m=1; testname=seqname;
    while(exist(strcat(testname,'.seq'),'file')>0)
    testname=strcat(seqname,'_',num2str(m));
    m=m+1;
    end
    seqname=testname;
end


%% Write to file

% The sequence is written to file in compressed form according to the file
% format specification using the |write| method.
%pn=which('myspiralTSE');
seq.write(strcat(seqname,'.seq'))
save(strcat('p',seqname),'spiralmode','scanmode','accmode','fatsat','acq','acqP','seg','segP','spiral','kt','ktraj_adc','system');
%%
% Display the first few lines of the output file
s=fileread(strcat(seqname,'.seq'));
disp(s(1:300))
disp(strcat('echo spacing:',num2str(acqP.TE*1000),' ms'))
%% plot gradients etc
if(plotflag(5)=='1'),
    fig=seq.plot('TimeRange',[0 (acqP.necho+1)*acqP.TE+4*acq.readoutTime],'timeDisp','ms');
end
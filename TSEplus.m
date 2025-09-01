%% Create a TSE sequence and export for execution
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
%   Juergen Hennig <juergen.hennig@uniklinik-freiburg.de>
%   Maxim Zaitsev  <maxim.zaitsev@uniklinik-freiburg.de>
%
%   Modes:
%   acqP.wavmode: weighted averaging mode
%       = 'cont': PE-steps for averaging are sorted into continuos array.
%       parameter: acqP.wav:
%                   array to specify weighted averaging, eg [1 0.7 0.3]
%
%   acqP.navmode: include navigator echos in first and/or last echo
%                 echotrain will be prolonged by number of navigators
%       values: x 0  navigator in first echo
%               0 x  navigator in last echo
%               x x  navigator in both echoes
%               x=1: readout direction, x=2: PE direction
%   PE-navigator in first echo is not compatible with acqP.T2prep='on'
%
%   acqP.T2prep: First echo interval with different echo spacing

%% Instantiation and gradient limits
% The system gradient limits can be specified in various units _mT/m_,
% _Hz/cm_, or _Hz/m_. However the limits will be stored internally in units
% of _Hz/m_ for amplitude and _Hz/m for slew. Unspecificied hardware
% parameters will be assigned default values.
clear count* acqP acq
plotflag='11111111';
dG=100e-6;
pulseflag_ex=0;
pulseflag_ref=0;
%%
system = mr.opts('MaxGrad', 45, 'GradUnit', 'mT/m', ...
    'MaxSlew', 130, 'SlewUnit', 'T/m/s', 'rfRingdownTime', 100e-6, ...
    'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6);
B0=3;
%%
% A new sequence object is created by calling the class constructor.
nseq=20;
for kseq=8
    seq=mr.Sequence(system);
    acqP.T2prep='off';
    acqP.fatsat='off';
    acqP.flipback='off';
    acqP.sat_ppm=-3.45;
    acqP.sat_freq=acqP.sat_ppm*1e-6*B0*system.gamma;
    %% Sequence events
    % Some sequence parameters are defined using standard MATLAB variables
    acqP.fov=220e-3;
    acqP.Nx=320;        acqP.Ny=320;        acqP.PEfac=acqP.Nx/acqP.Ny; %=1:same fov; =acqP.Nx/acqP.Ny:same res
    acqP.necho=20;
    acqP.samplingTime= 3.84e-3;                                           %should be a multiple of Nx in microseconds
    acqP.NSlices=1;     acqP.sliceGAP=1.2;  acqP.sliceThickness=2e-3;
    acqP.flipref=120;   acqP.flipflag=4;
    acqP.TE=10e-3;      acqP.TEeff=80e-3;   acqP.TEprep=40e-3;
    acqP.TR=4.66;       acqP.TI=0e-3;      acqP.TImod=1;
    acqP.HF_fac=0;      acqP.PEref=20;      acqP.GR_fac=1;
    acqP.PEtype='linear';
    acqP.nDummy=0;  % =0: dummy scan will be performed
    % acqP-parameters related to special acquisition modes.
    acqP.GRramp=0;
    acqP.nrep=1;
    acqP.navmode=[0 0];         % navigator echoes at beginning and end of echotrain
    acqP.wav=[1];    acqP.wavmode='cont';    acqP.PEinc=0;
    acqP.GDs=30;
    acqP.GDx=0;
    acqP.GDy=0;
    acqP.tD=0e-3;
    
    %%  acq-Parameters
    acq.tEx=2e-3;
    acq.tBwPex=2;
    acq.accfac=2;
    tExwd=acq.tEx+system.rfRingdownTime+system.rfDeadTime;
    acq.tRef=2.56e-3;
    acq.tBwPref=2.56;
    tRefwd=acq.tRef+system.rfRingdownTime+system.rfDeadTime;
    acq.tRefSE=5.12e-3;
    acq.tBwPrefSE=5.12;
    tRefwdSE=acq.tRefSE+system.rfRingdownTime+system.rfDeadTime;
    acq.rfTIdur=5e-3;
    acq.rfTItBwP=5;
    acq.GSexFac=1;
    acq.GSrefFac=1;
    acq.GSSEfac=1;
    acq.GSIRFac=1;
    acq.fspR=0.55;
    acq.fspS=0.55;
    rfex_phase=pi/2;
    rfref_phase=0;
    slofac=0;               %slice order factor
    if(nseq>1), myTSE_para; end
    %% Calculate dependant parameters
    if(strcmp((acqP.PEtype),'paired')&&(acqP.Ny/acqP.necho>1.5))
        disp('paired PE-mode only valid for single-shot acquisition, sequence construction terminated')
        return
    end
    readoutTime = acqP.samplingTime + 2*system.adcDeadTime;
    if(strcmp(acqP.T2prep,'on')), k0=round((acqP.TEeff-acqP.TEprep+acqP.TE)/acqP.TE);
    else
        k0=round(acqP.TEeff/acqP.TE);
    end
    if(strcmp(acqP.T2prep,'on')&&(acqP.navmode(1)==2))
        disp('PE-navigator in first echo is not compatible with acqP.T2prep=on, navigator will be in read direction')
        acqP.navmode(1)=0;
    end
    tExwd=acq.tEx+system.rfRingdownTime+system.rfDeadTime;
    tRefwd=acq.tRef+system.rfRingdownTime+system.rfDeadTime;
    tRefwdSE=acq.tRefSE+system.rfRingdownTime+system.rfDeadTime;
    tSpex=0.5*(acqP.TE-tExwd-tRefwd);
    tSp=0.5*(acqP.TE-readoutTime-tRefwd);
    acqP.OSlices=-[1:2:acqP.NSlices 2:2:acqP.NSlices]+ceil(acqP.NSlices/2);  
    if(slofac==1)
        acqP.OSlices=acqP.OSlices(end:-1:1);
    end
    if(acqP.tD==0),     acqP.GDs=0;    acqP.GDx=0;    acqP.GDy=0; acqP.tD=2*dG; end

    %%
    %%% Base gradients
    %%% Slice selection
    %
    % First, the slice selective RF pulses (and corresponding slice gradient)
    % are generated using the |makeSincPulse| function.
    % Gradients are recalculated such that their flattime covers the pulse plus
    % the rfdead- and rfringdown- times.
    % %
    flipex=90*pi/180;

    [rfex, gzex] = mr.makeSincPulse(flipex,system,'Duration',acq.tEx,...
        'sliceThickness',acqP.sliceThickness,'apodization',0.5,'timeBwProduct',acq.tBwPex,'PhaseOffset',rfex_phase,'use','excitation');
    rfex.delay=system.rfRingdownTime;
    GSex = mr.makeTrapezoid('z',system,'amplitude',gzex.amplitude,'FlatTime',tExwd,'riseTime',dG);
    GSex.amplitude=acq.GSexFac*GSex.amplitude;
    GSex.area=acq.GSexFac*GSex.area;
    %if(plotflag(1)=='1'), plotPulse(rfex,GSex); end
    GSexfirst=GSex.amplitude;
    acq.GSex=GSex.amplitude;

    [rfref, gzref] = mr.makeSincPulse(pi,system,'Duration',acq.tRef,...
        'sliceThickness',acqP.sliceThickness,'apodization',0.5,'timeBwProduct',acq.tBwPref,'PhaseOffset',rfref_phase,'use','refocusing');
    rfref.delay=system.rfRingdownTime;
    GSref = mr.makeTrapezoid('z',system,'amplitude',gzref.amplitude,'FlatTime',tRefwd,'riseTime',dG);
    GSref.amplitude=acq.GSrefFac*GSref.amplitude;
    GSref.area=acq.GSrefFac*GSref.area;
    acq.GSref=GSref.amplitude;
    refenvelope=rfref.signal;
    %%
    [rfrefSE, gzrefSE] = mr.makeSincPulse(pi,system,'Duration',acq.tRefSE,...
        'sliceThickness',acqP.sliceThickness,'apodization',0.5,'timeBwProduct',acq.tBwPrefSE,'PhaseOffset',rfref_phase,'use','refocusing');
    GSrefSE = mr.makeTrapezoid('z',system,'amplitude',gzref.amplitude,'FlatTime',tRefwdSE,'riseTime',dG);
    GSrefSE.amplitude=acq.GSSEfac*GSrefSE.amplitude;
    GSrefSE.area=acq.GSSEfac*GSrefSE.area;
    rfrefSE.delay=system.rfDeadTime;
    %if(plotflag(1)=='1'), plotPulse(rfref,GSref); end
    GSrefirst=GSref.amplitude;
    GSSErefirst=GSrefSE.amplitude;
    AGSex=GSex.area/2;
    AGSref=GSref.area/2;
    GSspr = mr.makeTrapezoid('z',system,'area',AGSref*(1+acq.fspS),'duration',tSp,'riseTime',dG);
    GSspex = mr.makeTrapezoid('z',system,'area',AGSref*(1+acq.fspS)-AGSex,'duration',tSpex,'riseTime',dG);
    GSspr_end = mr.makeTrapezoid('z',system,'amplitude',1e6,'duration',4e-3,'riseTime',dG);
    GSspr_end.delay=dG;
    %% Fatsat
    if(strcmp(acqP.fatsat,'on'))
        rf_fst=8e-3;
        if (B0<2), rf_fst=1e-5*floor(1e5*10e-3*1.5/B0); end
        rf_fs = mr.makeGaussPulse(110*pi/180,'system',system,'Duration',4e-3,...
            'bandwidth',abs(acqP.sat_freq),'freqOffset',acqP.sat_freq,'use','other');
        gz_fs = mr.makeTrapezoid('z',system,'delay',mr.calcDuration(rf_fs),'Area',1/1e-4); % spoil up to 0.1mm
        fs_dur=gz_fs.delay+gz_fs.riseTime+gz_fs.flatTime+gz_fs.fallTime;
    else
        fs_dur=0;
    end

    %% Inversion recovery
    if(acqP.TI>0)
        
        [rfIR, gzIR] = mr.makeSincPulse(pi,system,'Duration',acq.rfTIdur,...
            'sliceThickness',acqP.sliceThickness,'apodization',0.5,'timeBwProduct',acq.rfTItBwP,'PhaseOffset',0,'delay',system.rfDeadTime,'use','other');
        rfIR.delay=rfIR.delay+dG;
        GSIR = mr.makeTrapezoid('z',system,'amplitude',gzIR.amplitude,'FlatTime',rfIR.shape_dur+system.rfRingdownTime+system.rfDeadTime,'riseTime',dG);
        GSIR.amplitude=acq.GSIRFac*GSIR.amplitude;
        GSIRsp=mr.makeTrapezoid('z',system,'amplitude',0.5*system.maxGrad,'duration',4e-3,'riseTime',dG,'delay',mr.calcDuration(GSIR)-dG);
        GSIRtot=mr.addGradients({GSIR,GSIRsp},'system',system);
        GSIRtot_times=[0; GSIRtot.tt+dG; GSIRtot.shape_dur+2*dG];
        GSIRtot_gradients=[0; GSIRtot.waveform; 0];
        GSIRtot=mr.makeExtendedTrapezoid('z','times',GSIRtot_times,'amplitudes',GSIRtot_gradients);
        TI_dur=GSIRtot.shape_dur; TI_t2P=rfIR.delay+rfIR.shape_dur/2;
        acqP.nDummy=0;
    else
        TI_dur=0;

    end
    
    

    %%
    %%% Readout gradient
    % To define the remaining encoding gradients we need to calculate the
    % $k$-space sampling. The Fourier relationship
    %
    % $$\Delta k = \frac{1}{acqP.fov}$$
    %
    % Therefore the area of the readout gradient is $n\Delta k$.
    deltak=1/acqP.fov;
    kWidth = acqP.Nx*deltak;

    GRacq = mr.makeTrapezoid('x',system,'FlatArea',kWidth,'FlatTime',readoutTime,'riseTime',dG);
    adc = mr.makeAdc(acq.accfac*acqP.Nx,'Duration',acqP.samplingTime, 'Delay', system.adcDeadTime);%,'Delay',GRacq.riseTime);
    GRspr = mr.makeTrapezoid('x',system,'area',GRacq.area*acq.fspR,'duration',tSp,'riseTime',dG);
    AGRspr=GRspr.area;%GRacq.area/2*acq.fspR;
    
    


    %%
    %%% Phase encoding
    % To move the $k$-space trajectory away from 0 prior to the readout a
    % prephasing gradient must be used. Furthermore rephasing of the slice
    % select gradient is required.

    nex=floor(acqP.Ny/acqP.necho);
  
    Ny0=acqP.Ny;
    acqP.PEorder=[];acqP.PEcount=[];acqP.PEindex=[];
    clear nexc;
    
    nav=length(acqP.wav);
    for kav=1:nav
      %  for kav=1
            [PEorder,PEindex,NS] = PE_order(acqP.PEtype,floor(acqP.Ny*acqP.wav(kav)),acqP.necho,k0,[],acqP.HF_fac*acqP.wav(kav),acqP.PEref,acqP.GR_fac);
            if(NS==1), acqP.necho=length(PEindex); end
            PEcount=0*PEorder+kav;
            acqP.PEorder=[acqP.PEorder PEorder]+acqP.PEinc*kav-1;
            acqP.PEcount=[acqP.PEcount PEcount];
            acqP.PEindex=[acqP.PEindex acqP.PEindex];
            nexc(kav)=NS;
    end

    nex=sum(nexc);
    if(nav>1)
    if(strcmp(acqP.wavmode,'cont'))
        [PEtemp,ind]=sort(acqP.PEorder(:));
        count=acqP.PEcount(:);
        count=count(ind);
        [PEorder,PEindex,NS] = PE_order(acqP.PEtype,floor(acqP.Ny*acqP.wav(kav)),acqP.necho,k0,PEtemp',acqP.HF_fac*acqP.wav(kav),acqP.PEref,acqP.GR_fac);
        acqP.PEorder=PEorder;
        acqP.PEindex=PEindex;
    end
    end
    if(plotflag(1)=='1')
        figure
        hold on
        plot(acqP.PEorder,'.-')
        temp=ylim;
        plot([k0 k0],[temp(1) 0.8*temp(2)],'k-')
        text(k0,0.9*temp(2),'k0')
        title(strcat('PE-order (k vs. Necho) ,',acqP.PEtype))
    end
    acqP.nexc=nexc;
    phaseAreas= acqP.PEorder*deltak*acqP.PEfac;
    acq.phaseAreas=phaseAreas;

    %% split gradients and recombine into blocks
    % lets start with slice selection....
    GS1times=[0 GSex.riseTime];
    GS1amp=[0 GSexfirst];
    GS1 = mr.makeExtendedTrapezoid('z','times',GS1times,'amplitudes',GS1amp);
    GS1rew_amp=[GSexfirst 0];
    GS1rew = mr.makeExtendedTrapezoid('z','times',GS1times,'amplitudes',GS1rew_amp);

    GS2times=[0 GSex.flatTime];
    GS2amp=[GSexfirst GSexfirst];
    GS2 = mr.makeExtendedTrapezoid('z','times',GS2times,'amplitudes',GS2amp);

    GS3times=[0 GSspex.riseTime GSspex.riseTime+GSspex.flatTime GSspex.riseTime+GSspex.flatTime+GSspex.fallTime];
    GS3amp=[GSexfirst GSspex.amplitude GSspex.amplitude GSrefirst];
    GS3 = mr.makeExtendedTrapezoid('z','times',GS3times,'amplitudes',GS3amp);

    GS3rewamp=[GSrefirst GSspex.amplitude GSspex.amplitude GSexfirst];
    GS3rew = mr.makeExtendedTrapezoid('z','times',GS3times,'amplitudes',GS3rewamp);

    GS4times=[0 GSref.flatTime];
    GS4amp=[GSrefirst GSrefirst];
    GS4 = mr.makeExtendedTrapezoid('z','times',GS4times,'amplitudes',GS4amp);

    GS4SEtimes=[0 GSrefSE.flatTime];
    GS4SEamp=[GSSErefirst GSSErefirst];
    GS4SE = mr.makeExtendedTrapezoid('z','times',GS4SEtimes,'amplitudes',GS4SEamp);

    GS5times=[0 GSspr.riseTime GSspr.riseTime+GSspr.flatTime GSspr.riseTime+GSspr.flatTime+GSspr.fallTime];
    GS5amp=[GSrefirst GSspr.amplitude GSspr.amplitude 0];
    GS5 = mr.makeExtendedTrapezoid('z','times',GS5times,'amplitudes',GS5amp);
    GS5SEamp=[GSSErefirst GSspr.amplitude GSspr.amplitude 0];
    GS5SEtimes=[0 GSspr.riseTime GSspr.riseTime+GSspr.flatTime GSspr.riseTime+GSspr.flatTime+GSspr.fallTime];
    GS5SE = mr.makeExtendedTrapezoid('z','times',GS5SEtimes,'amplitudes',GS5SEamp);

    GS7times=[0 GSspr.riseTime GSspr.riseTime+GSspr.flatTime GSspr.riseTime+GSspr.flatTime+GSspr.fallTime];
    GS7amp=[0 GSspr.amplitude GSspr.amplitude GSrefirst];
    GS7 = mr.makeExtendedTrapezoid('z','times',GS7times,'amplitudes',GS7amp);
    GS7SEamp=[0 GSspr.amplitude GSspr.amplitude GSSErefirst];
    GSarea=calcArea(GS4)+calcArea(GS5)+calcArea(GS7);
    GS7SEtimes=[0 GSspr.riseTime GSspr.riseTime+GSspr.flatTime GSspr.riseTime+GSspr.flatTime+GSspr.fallTime];
    %GS7SEtimes=[0 GSspr.riseTime GSspr.riseTime+GS5SEflatt GSspr.riseTime+GS5SEflatt+GSspr.fallTime];
    GS7SE = mr.makeExtendedTrapezoid('z','times',GS7SEtimes,'amplitudes',GS7SEamp);

    if((acqP.navmode(1)==3)||(acqP.navmode(2)==3))
        GSpredur=system.gradRasterTime*floor(readoutTime/4/system.gradRasterTime);
        GSread = mr.makeTrapezoid('z',system,'Amplitude',5e5,'Duration',readoutTime-2*GSpredur,'riseTime',dG);
        GSprew=mr.makeTrapezoid('z',system,'Area',-GSread.area/2,'Duration',GSpredur,'riseTime',dG);
        GSacqtimes=[0 dG GSpredur-dG GSpredur GSpredur+dG GSpredur+dG+GSread.flatTime readoutTime-GSpredur readoutTime-GSpredur+dG readoutTime-dG readoutTime]
        GSacqamp=[0 GSprew.amplitude GSprew.amplitude 0 GSread.amplitude GSread.amplitude 0 GSprew.amplitude GSprew.amplitude 0 ]
        %plot(GSacqtimes,GSacqamp)
        GSacq = mr.makeExtendedTrapezoid('z','times',GSacqtimes,'amplitudes',GSacqamp);
    end

    if(strcmp(acqP.T2prep,'on'))
        GS3p1=mr.makeExtendedTrapezoid('z','times',[0 GSex.fallTime],'amplitudes',[GSex.amplitude 0]);
        GS3p2=mr.makeTrapezoid('z',system,'area',-AGSex+GSarea/2,'duration',tSpex,'riseTime',dG);
        % and now the readout gradient....
        GS5p1 = mr.makeTrapezoid('z',system,'Area',GSarea/2,'Duration',tSp);
    end

    GR5times=[0 GRspr.riseTime GRspr.riseTime+GRspr.flatTime GRspr.riseTime+GRspr.flatTime+GRspr.fallTime];
    GR5amp=[0 GRspr.amplitude GRspr.amplitude GRacq.amplitude*(1+acqP.GRramp)];
    GR5 = mr.makeExtendedTrapezoid('x','times',GR5times,'amplitudes',GR5amp);

    GR6times=[0 readoutTime/2 readoutTime];
    GR6amp=[GRacq.amplitude*(1+acqP.GRramp) GRacq.amplitude*(1-acqP.GRramp) GRacq.amplitude*(1+acqP.GRramp)];

    GR6 = mr.makeExtendedTrapezoid('x','times',GR6times,'amplitudes',GR6amp);

    GR7times=[0 GRspr.riseTime GRspr.riseTime+GRspr.flatTime GRspr.riseTime+GRspr.flatTime+GRspr.fallTime];
    GR7amp=[GRacq.amplitude*(1+acqP.GRramp) GRspr.amplitude GRspr.amplitude 0];
    GR7 = mr.makeExtendedTrapezoid('x','times',GR7times,'amplitudes',GR7amp);

    GRarea=calcArea(GR5)+calcArea(GR6)+calcArea(GR7);
    GRpreph = mr.makeTrapezoid('x',system,'Area',GRarea/2,'duration',tSpex,'riseTime',dG);
    GR3=GRpreph;
    GRrew = mr.makeTrapezoid('x',system,'area',GRarea/2,'duration',tSp,'riseTime',dG);

    GRpre_s=mr.makeTrapezoid('x',system,'Area',GRarea/2,'duration',tSpex,'riseTime',dG);

    % and now the PE-navigator
    GPacq = mr.makeTrapezoid('y',system,'FlatArea',kWidth,'FlatTime',readoutTime,'riseTime',dG);
    GPRpre = mr.makeTrapezoid('y',system,'Area',-GPacq.area/2,'Duration',tSp-dG,'riseTime',dG);
    GPR5times=[0 dG dG+GPRpre.flatTime 2*dG+GPRpre.flatTime 3*dG+GPRpre.flatTime];
    GPR5amp=[0 GPRpre.amplitude GPRpre.amplitude 0 GPacq.amplitude];
    GPR5 = mr.makeExtendedTrapezoid('y','times',GPR5times,'amplitudes',GPR5amp);

    GPRacqtimes=[0 GPacq.flatTime];
    GPRacqamp=[GPacq.amplitude GPacq.amplitude];
    GPRacq = mr.makeExtendedTrapezoid('y','times',GPRacqtimes,'amplitudes',GPRacqamp);

    GPR7times=[0 dG 2*dG 2*dG+GPRpre.flatTime 3*dG+GPRpre.flatTime];
    GPR7amp=[GPacq.amplitude 0 GPRpre.amplitude GPRpre.amplitude 0 ];
    GPR7 = mr.makeExtendedTrapezoid('y','times',GPR7times,'amplitudes',GPR7amp);
    nPEecho=acqP.necho;
    if(acqP.navmode(1)>0), acqP.necho=acqP.necho+1; end;
    PEshift=acqP.necho-nPEecho;
    if(acqP.navmode(2)>0), acqP.necho=acqP.necho+1; end;
    temp=zeros([acqP.necho nex]);
    temp(PEshift+1:PEshift+nPEecho,:)=acq.phaseAreas;
    acq.phaseAreas=temp;

    %% Calculate filltimes
    acq.Extime=mr.calcDuration(GS1)+mr.calcDuration(GS2)+mr.calcDuration(GS3);
    acq.Reftime=mr.calcDuration(GS4)+mr.calcDuration(GS5)+mr.calcDuration(GS7)+readoutTime;
    tEx2P=rfex.delay+rfex.shape_dur/2;
    if(acqP.flipback==0)
        tend=mr.calcDuration(GS4)+mr.calcDuration(GS5)+mr.calcDuration(GSspr_end);
    else
        tend=mr.calcDuration(GS2)+mr.calcDuration(GS3)+mr.calcDuration(GS4)+mr.calcDuration(GS1rew)+mr.calcDuration(GSspr_end);
    end
    tETrain=acq.Extime+acqP.necho*acq.Reftime+tend+fs_dur+TI_dur;
    if(strcmp(acqP.T2prep,'on'))
        tETrain=tETrain+acqP.TEeff;
    end
    if(acqP.NSlices*tETrain>(acqP.TR)),
        acqP.TR=acqP.NSlices*tETrain;
        disp(strcat('Warning!!! acqP.TR too short, adapted to include all slices to : ',num2str(acqP.TR*1000),' ms'));
    end
    TR_tot=(acqP.TR-acqP.NSlices*tETrain)/acqP.NSlices;
    TIfill=0;
    IRSlices=acqP.OSlices;
    if(acqP.TI>0)
        nTI=floor(acqP.TI/(tETrain+TR_tot));
        IRSlices=circshift(acqP.OSlices,-2);
        if(nTI==0), acqP.TImod=1; end
        if(acqP.TImod==0)
            acqP.TI=nTI*tETrain;
            disp(strcat('acqP.TI adapted t0 : ',num2str(acqP.TI*1000),' ms'));
        else
            TIfill=acqP.TI-nTI*tETrain;
            if(TIfill>TR_tot), TR_tot=(acqP.TI-nTI*tETrain)/(nTI+1);
                TIfill=TR_tot-tEx2P-TI_dur+TI_t2P;                  % TIfill is corrected for pulse timing
                acqP.TR=acqP.NSlices*(tETrain+TR_tot);
                disp(strcat('acqP.TR adapted t0 : ',num2str(acqP.TR*1000),' ms'));
            end
        end
    end
    TRfill=TR_tot-TIfill;

    % round to gradient raster
    TIfill=system.gradRasterTime * round(TIfill / system.gradRasterTime);
    TRfill=system.gradRasterTime * round(TRfill / system.gradRasterTime);
    delayTR = mr.makeDelay(TRfill);
    delayTI = mr.makeDelay(TIfill);
    %% Fill times for SE-TSE
    if(strcmp(acqP.T2prep,'on'))
        TEfill1=acqP.TEprep/2-GS2.shape_dur/2-GS3p1.shape_dur-mr.calcDuration(GS3p2)-mr.calcDuration(GRspr)-GS4SE.shape_dur/2;
        TEfill2=acqP.TEprep/2-GS5SE.shape_dur-GS4SE.shape_dur/2+GS4.shape_dur/2-acqP.TE/2;
        
        if(TEfill1<=0),
            fprintf(2,'\n\nTEprep too short, sequence creation terminated\n')
            return
        end
        if(TEfill1<acqP.tD), acqP.tD=TEfill1; 
        disp(strcat('Warninng, acqP.tD has been reduced to:',num2str(acqP.tD*1000),' ms'))
        end
       
        delayTE1 = mr.makeDelay(TEfill1);
        delayTE2 = mr.makeDelay(TEfill2);
    end
        %% Diffusion gradients
    if(strcmp(acqP.T2prep,'on'))
    GDs=mr.makeTrapezoid('z',system,'amplitude',acqP.GDs/1000*system.gamma,'duration',acqP.tD,'riseTime',dG,'delay',dG);
    GDx=mr.makeTrapezoid('x',system,'amplitude',acqP.GDx/1000*system.gamma,'duration',acqP.tD,'riseTime',4*dG,'delay',dG);
    GDy=mr.makeTrapezoid('y',system,'amplitude',acqP.GDy/1000*system.gamma,'duration',acqP.tD,'riseTime',dG,'delay',dG);
    GDxr=GDx; GDxr.delay=TEfill2-acqP.tD-dG;
    GDyr=GDy; GDyr.delay=TEfill2-acqP.tD-dG;
    GDsr=GDs; GDsr.delay=TEfill2-acqP.tD-dG;
    end

    %% calculate flip angles
    if(acqP.flipflag<=1)
        rf=acqP.flipref+zeros([1 acqP.nrep.*acqP.necho]); end

    if(acqP.flipflag==1),  rf(1)=90+acqP.flipref/2;end
    
    if(acqP.flipflag==2), 
        kz=k0;
        if(kz<5),kz=5; end;
        flipend=acqP.flipref;
        [rf,~] = TRAPS_flip(acqP.flipref,(acqP.nrep).*acqP.necho,kz,'opt',2,acqP.flipref,flipend,[kz kz kz acqP.necho]);
    end

    if(acqP.flipflag==3), flipend=60;
        kz=k0;
        if(kz<5),kz=5; end
        [rf,~] = TRAPS_flip(acqP.flipref,(acqP.nrep).*acqP.necho,k0,'none',0,acqP.flipref,flipend,[kz kz kz acqP.necho]);
        rf(1)=(180+acqP.flipref)/2;
    end

    if(acqP.flipflag==4), flipend=60;
        kz=k0;
        if(kz<5),kz=5; end
        [rf,~] = TRAPS_flip(acqP.flipref,(acqP.nrep)*acqP.necho,k0,'opt',2,acqP.flipref,flipend,[6 6 kz acqP.necho]);
    end

    if(acqP.flipflag==5), 
        ka=k0;
        if(ka==1); ka=2; end
        flipend=60;
        [rf,~] = TRAPS_flip(acqP.flipref,(acqP.nrep)*acqP.necho,k0,'opt',2,180,60,[ceil(ka/2) ka ka acqP.necho]);
        
    end

    rf=rf(1:acqP.nrep*acqP.necho);
    acqP.flip=rf(1:acqP.nrep*acqP.necho);
    if(plotflag(2)=='1')
        figure
        plot(acqP.flip,'.-'), axis([xlim 0 180]);
        title('flip angle along echotrain')
    end
    % the following lines can be usedto calculate the total RF-power of the
    % sequence in a.u. For comparison of different sequences the calculated
    % values can be referenced to a reference sequence (e.g. a sequence using
    % 180

    [RFP,RFPtr,RFPrel] = calcRFP(rfex,rfref,acqP);
    acqP.pow=RFPrel;
    disp(strcat('rel. power compared to 180° refocussing :',num2str(acqP.pow)));
    rfpower=acqP.pow;
    %% Define sequence blocks
    % Next, the blocks are put together to form the sequence
    for kex=acqP.nDummy:nex % MZ: we start at 0 to have one dummy
        for s=1:acqP.NSlices
            rfex.freqOffset=GSex.amplitude*acqP.sliceGAP*acqP.sliceThickness*acqP.OSlices(s);
            rfref.freqOffset=GSref.amplitude*acqP.sliceGAP*acqP.sliceThickness*acqP.OSlices(s);
            rfrefSE.freqOffset=GSrefSE.amplitude*acqP.sliceGAP*acqP.sliceThickness*acqP.OSlices(s);
            rfex.phaseOffset=rfex_phase-2*pi*rfex.freqOffset*mr.calcRfCenter(rfex); % align the phase for off-center slices
            rfref.phaseOffset=rfref_phase-2*pi*rfref.freqOffset*mr.calcRfCenter(rfref); % dito
            rfrefSE.phaseOffset=rfref_phase-2*pi*rfrefSE.freqOffset*mr.calcRfCenter(rfrefSE); % dito
            seq.addBlock(delayTR);
            if(acqP.TI>0),
                rfIR.freqOffset=GSIR.amplitude*acqP.sliceGAP*acqP.sliceThickness*IRSlices(s);
                rfIR.phaseOffset=-2*pi*rfref.freqOffset*mr.calcRfCenter(rfIR);
                seq.addBlock(rfIR,GSIRtot);
                seq.addBlock(delayTI);
            end
            if(strcmp(acqP.fatsat,'on')), seq.addBlock(rf_fs,gz_fs); end
            if(strcmp(acqP.T2prep,'on'))
                seq.addBlock(GS1);
                seq.addBlock(GS2,rfex);
                seq.addBlock(GS3p1);
                seq.addBlock(GS3p2,GRpre_s);
                %seq.addBlock(GDx,GDy,GDs);
                seq.addBlock(delayTE1,GDx,GDy,GDs);
                seq.addBlock(GS7SE,GRspr);
                seq.addBlock(GS4SE,rfrefSE);
                seq.addBlock(GS5SE,GRspr);
                seq.addBlock(delayTE2,GDxr,GDyr,GDsr);
               % seq.addBlock(GDxr,GDyr,GDsr);
            else
                seq.addBlock(GS1);
                seq.addBlock(GS2,rfex);
                seq.addBlock(GS3,GR3);

            end
            phaseArea=0;
            for kech=1:acqP.necho
                if (kex>0)
                    phaseArea=acq.phaseAreas(kech,kex);
                end
                GPpre = mr.makeTrapezoid('y',system,'Area',phaseArea,'Duration',tSp);
                GPrew = mr.makeTrapezoid('y',system,'Area',-phaseArea,'Duration',tSp);
                rfref.signal=refenvelope*acqP.flip(kech)/180;
                if((kech==1)&&strcmp(acqP.T2prep,'on'))
                    seq.addBlock(GR5,GPpre,GS5p1);
                    seq.addBlock(GR6,adc);
                    seq.addBlock(GS7,GR7,GPrew);
                else
                    seq.addBlock(GS4,rfref);
                    if((kech==1)&&acqP.navmode(1)==2)
                        seq.addBlock(GS5,GRrew,GPR5);
                        seq.addBlock(GPRacq,adc);
                        seq.addBlock(GS7,GRrew,GPR7);
                    elseif((kech==acqP.necho)&&acqP.navmode(2)==2)
                        seq.addBlock(GS5,GRrew,GPR5);
                        seq.addBlock(GPRacq,adc);
                        seq.addBlock(GS7,GRrew,GPR7);
                    elseif((kech==acqP.necho)&&acqP.navmode(2)==3)
                        seq.addBlock(GS5,GRrew);
                        seq.addBlock(GSacq,adc);
                        seq.addBlock(GS7,GRrew);
                    else
                        seq.addBlock(GS5,GR5,GPpre);
                        seq.addBlock(GR6,adc);
                        seq.addBlock(GS7,GR7,GPrew);
                    end

                end

            end
            if(strcmp(acqP.flipback,'on'))
                rfex.phaseOffset=pi+rfex_phase-2*pi*rfex.freqOffset*mr.calcRfCenter(rfex);
                seq.addBlock(GS4,rfref);
                seq.addBlock(GS3rew);
                seq.addBlock(GS2,rfex);
                seq.addBlock(GS1rew);
            else
                seq.addBlock(GS4);
                seq.addBlock(GS5);
            end
            seq.addBlock(GSspr_end);



        end
        
    end
    acqP.necho=nPEecho; %reset number of echoes

    %% check whether the timing of the sequence is correct
    [ok, error_report]=seq.checkTiming;

    if (ok)
        fprintf('Timing check passed successfully\n');
    else
        fprintf('Timing check failed! Error listing follows:\n');
        fprintf([error_report{:}]);
        fprintf('\n');
    end

    %% k-space trajectory calculation
    [ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();
    %%
    % plot k-spaces
    if(plotflag(3)=='1')
        figure; plot(t_ktraj,ktraj'); axis([TRfill TRfill+acqP.TR/acqP.NSlices ylim]); title('k-space components as functions of time'); end% plot the entire k-space trajectory
    if(plotflag(4)=='1')
        figure; plot(ktraj(1,:),ktraj(2,:),'b',...
            ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % a 2D plot
        axis('equal'); % enforce aspect ratio for the correct trajectory display
        title('2D k-space (x-y)');
    end
    if(plotflag(5)=='1')
        figure; plot(ktraj(1,:),ktraj(3,:),'b',...
            ktraj_adc(1,:),ktraj_adc(3,:),'r.'); % a 2D plot
        axis('equal'); % enforce aspect ratio for the correct trajectory display
        title('2D k-space(x-z)');
    end
 
    %% Write to file

    % The sequence is written to file in compressed form according to the file
    % format specification using the |write| method.
    if(~exist('seqno')), seqno=1; end
    seqname=strcat('TSE_',num2str(seqno));
    seq.write(strcat(seqname,'.seq'));
    seqno=seqno+1;
    save(seqname,'system','pulseflag_ex','pulseflag_ref','acqP','acq');
    %%

    if(plotflag(6)=='1')
        fig=seq.plot('TimeRange',[TRfill TRfill+0.5],'timeDisp','ms','stacked',1);
    end
    % seq.install('siemens');
   
end
[~, pns_n, pns_c, tpns]=seq.calcPNS('MP_GPA_K2309_2250V_951A_AS82.asc'); % prisma
% [pns_ok, pns_n, pns_c, tpns]=seq.calcPNS('MP_GradSys_P034_X60.asc'); % Cima.X PNS

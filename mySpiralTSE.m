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
% # Create slice selective RF pulse for imaging
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
%       'outin' : first half 180ï¿½ rotated and out-in
%       'sym'   : first half mirrored at y-axis and out-in;
%
%   accmode=    acceleration mode
%       'dd'    : dual density. First seg.n_full*acq.nadc/2 points are fully sampled,
%                 others undersampled by a factor spiral.Nout
%       'vd'    : variable density with Fcoeff2 = - Fcoeff/spiral.Nout
%
%   fatsat
%       'on'
%       'no'    :  with fatsat gradient, but no rf
%       'off'
%   T2prep
%       'on'    : First echo is different from CPMG-train. Usefúll for fMRI
%                and quantitative T2-measurements
%
%   seq_var_mod     switch for specific parameters. Some of the switches require to set kseq>1
%                   Parameters for the different modes are defined in mySpiralTSE_par.m
%
%       seqvar_mod='none';      normal execution
%       seqvar_mod='T1var';     generates sequence with multiple
%                               acquisition of individual slices for T1-calculation
%                               acquisition scheme is determined by
%                               acqP.TRfac, which is an array representing
%                               the slice interleaving scheme.
%       seqvar_mod='MTvar';     used for assessment of MT-effect. acqP.MTfac determines the position of central slice
%       seqvar_mod='TErep';     number of images acquired on a single echotrain.              Used to measure long T2.
%       seqvar_mod='TEvar';     Generates sequence with TE-variation by concatenation
%
%   plotflag:   sequence of 6 binary digits to define, which plots are shown
%       1   :   plots from vds-package
%       2   :   trajectory as cloud of dots to check completeness
%       3   :   plot of individual segments
%       4   :   x-y-trajectory during readout(red) and total trajectory (blue)
%       5   :   x-z-trajectory during readout(red) and total trajectory (blue)
%       6   :   sequence
%% Instantiation and gradient limits
% The system gradient limits can be specified in various units _mT/m_,
% _Hz/cm_, or _Hz/m_. However the limits will be stored internally in units
% of _Hz/m_ for amplitude and _Hz/m/s_ for slew. Unspecificied hardware
% parameters will be assigned default values.


dG=100e-6;
seqname='TSE';
plotflag='000111';
initmode='no';

%count0=20;
count=10;
seqvar_mod='none';
% seqvar_mod='Diff';
% seqvar_mod='T1var';
% seqvar_mod='MTvar';
% seqvar_mod='TErep';
% seqvar_mod='TEvar';
% seqvar_mod='TEall';
% seqvar_mod='accvar';
% seqvar_mod='fsvar';
system = mr.opts('MaxGrad',35, 'GradUnit', 'mT/m', ...
    'MaxSlew', 170, 'SlewUnit', 'T/m/s', 'rfRingdownTime', 100e-6, ...
    'rfDeadTime', 100e-6);
B0=3;
try seq=mr.Sequence(system); end
warning('OFF', 'mr:restoreShape')
nseq=1;
for kseq=1:nseq

    tic
    try seq=mr.Sequence(system); end
    kconc=1;
    %for kconc=1:1,
    count=count+1;
    %% Sequence events
    % Some sequence parameters are defined using standard MATLAB variables
    setpar=1;   % if setpar =0 parameters should be read from file

    if(setpar==1)
        acq.slewfac=1;
        acq.gradfac=0.99;

        acqP.sat_ppm=-3.45;
        acqP.sat_freq=acqP.sat_ppm*1e-6*B0*system.gamma;

        segmode='fix';
        spmode='cont';
        scanmode='run';
        fatsat='on';
        accmode='vd';
        T2prep='on';

        temp=sprintf('%s %s','fatsat:',fatsat);
        disp(temp)
        temp=sprintf('%s %s','T2prep:',T2prep);
        disp(temp)

        if(strcmp(scanmode,'trim')), seqname=strcat(seqname,'_',segmode,'_',spmode,'_trim'); end
        acqP.fov=240e-3;
        spiral.Ninx=120;
        spiral.kOffset=50;
        spiral.Nin=1;
        spiral.Nout=1.4;
        acq.dninc=2;
        acq.accfac=4;
        %%
        acqP.NSlices=1;
        acqP.TRfac=1;
        acqP.MTfac=1;
        SLfac=1;
        acqP.sliceThickness=4e-3*SLfac; acqP.sliceGAP=1.5; acqP.nex=1;

        seg.n_full=1;

        acqP.sltemp=[1:2:acqP.NSlices 2:2:acqP.NSlices];
        if(strcmp(seqvar_mod,'T1var'))
            [temp] = T1var_index(acqP.NSlices,acqP.TRfac);
            acqP.sltemp=temp;
            acqP.Oslices=acqP.sltemp-ceil(acqP.NSlices/2);
            acqP.NSlices=length(temp);
        else
            if(strcmp(seqvar_mod,'MTvar'))
                [temp] = MTvar_index(acqP.sltemp,acqP.MTfac);
                acqP.sltemp=temp;
                acqP.Oslices=acqP.sltemp-ceil(acqP.NSlices/2);
                acqP.NSlices=length(temp);
            else
                acqP.Oslices=acqP.sltemp-ceil(acqP.NSlices/2);
            end
        end
        %%
        acqP.dummy=1;
        acqP.flipref=60; acqP.flipflag=2;
        acqP.TR=500e-3;
        %%
        acqP.TE=10e-3;
        acqP.TEprep=40e-3;
        acqP.nTE=1;
        if(strcmp(T2prep,'off'))
            acqP.TEeff=acqP.nTE*acqP.TE;
        else
            acqP.TEeff=acqP.TEprep; acqP.nTE=1;
        end
        acqP.nrep=1;
        %%
        acqP.PEtype='linear';
        segP.tEx=2.5e-3;
        segP.tExwd=segP.tEx+system.rfRingdownTime+system.rfDeadTime;
        segP.tRef=2e-3;
        segP.tRefwd=segP.tRef+system.rfRingdownTime+system.rfDeadTime;
        segP.tSp=1.05e-3;
        segP.tSpex=0.5*(acqP.TE-segP.tExwd-segP.tRefwd);
        segP.t1=1e-3;
        segP.t2=0.5e-3;
        segP.fspS=0.5;
        segP.TEprep=0.9e-3;
        segP.GSfac=1;
        segP.GXfac=1;
        segP.GYfac=1;

    end
    %%
    if(kconc==2), seqtype=strcat(sprintf('%03s',num2str(count)),'TSC'); end

    if(kseq>0)
        quitflag=0;
        mySpiralTSE_par
        if(quitflag==1), return; end
    end
    if(exist(strcat(seqname,'.seq'),'file')>0)
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
    %%

    rf_fst=8e-3;
    if (B0<2), rf_fst=1e-5*floor(1e5*10e-3*1.5/B0); end
    rf_fs = mr.makeGaussPulse(110*pi/180,'system',system,'Duration',2.5e-3,...
        'bandwidth',abs(acqP.sat_freq),'freqOffset',acqP.sat_freq);
    gz_fs = mr.makeTrapezoid('z',system,'delay',mr.calcDuration(rf_fs),'Area',1/1e-4); % spoil up to 0.1mm
    %figure; plot(rf_fs.t,rf_fs.signal)
    %%
    %%% Readout gradient
    % To define the remaining encoding gradients we need to calculate the
    % $k$-space sampling. The Fourier relationship
    %
    % $$\Delta k = \frac{1}{acqP.fov}$$
    %
    % Therefore the area of the readout gradient is $n\Delta k$.

    spiral.deltak=1/acqP.fov;
    spiral.kWidth = 2*spiral.Ninx*spiral.deltak;

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
    if(strcmp(accmode,'vd')), Fcoeff=[acqP.fov*100 -acqP.fov*100/(spiral.Nout)]; end

    res=100*acqP.fov/spiral.Ninx;
    rmax = 2*spiral.kWidth/100;		% cm^(-1), corresponds to 1mm resolution.
    rmax=1/res;

    disp('Calculating Gradient');

    [k,~,~,time,r,theta] = vds(smax,gmax,T/acq.accfac,spiral.Nin,Fcoeff,rmax);
    if ((spiral.Nin>1)&&strcmp(accmode,'dd'))
        [k1,~,~,time,r,theta] = vds(smax,gmax,T/acq.accfac,1,Fcoeff,rmax);
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
    if(strcmp(segmode,'single'))
        acqP.necho=1; seg.i1(1)=floor(acq.ntot/2);
        seg.ikseg=zeros([acqP.necho 2]);
        seg.ikseg(1,1)=1; seg.ikseg(1,2)=floor(acq.ninc/2);
        seg.nkseg=acq.ninc*ones([acqP.necho 1]); seg.nkseg(1)=floor(acq.ninc/2);
    end


    if(strcmp(segmode,'tan'))
        [seg.kseg, seg.nkseg, seg.ikseg, acqP.necho, seg.i1] = spiral_seg(kt,acq.ninc,spiral.kStart,acq.ntot,spmode);
        axlim=max(spiral.kStart,spiral.kWidth/2);
        axis([-axlim axlim -axlim axlim])
        if(plotflag(3)=='0'), close; end
    end

    if(strcmp(segmode,'fix'))

        acqP.necho=floor(np/acq.ninc);
        seg.nkseg=acq.ninc*ones([acqP.necho 1]); seg.nkseg(1)=floor(acq.ninc/2);
        seg.i1=0*seg.nkseg+acq.dninc; seg.i1(1)=floor(acq.ntot/2);
        seg.ikseg=zeros([acqP.necho 2]);
        seg.ikseg(1,1)=1; seg.ikseg(1,2)=floor(acq.ninc/2);

        for k=2:acqP.necho
            seg.ikseg(k,1)=seg.ikseg(k-1,2);
            seg.ikseg(k,2)=seg.ikseg(k,1)+seg.nkseg(k)-1;
        end
        if(strcmp(spmode,'alt')),
            temp=fliplr(seg.ikseg);
            seg.ikseg(2:2:end,:)=temp(2:2:end,:);
        end

        if(plotflag(3)=='1')
            figure
            axlim=max(spiral.kStart,spiral.kWidth/2);
            axis([-axlim axlim -axlim axlim])

            axis square
            hold on
            col=['b-';'r-'];
            if(strcmp(spmode,'cont'))
                for k=1:2:acqP.necho
                    plot(kt(seg.ikseg(k,1):seg.ikseg(k,2),1),(kt(seg.ikseg(k,1):seg.ikseg(k,2),2)),col(1,:));
                    plot(kt(seg.ikseg(k,2),1),kt(seg.ikseg(k,2),2),'bo');
                end
                col=flipud(col);
                for k=2:2:acqP.necho
                    plot(kt(seg.ikseg(k,1):seg.ikseg(k,2),1),(kt(seg.ikseg(k,1):seg.ikseg(k,2),2)),col(1,:));
                    plot(kt(seg.ikseg(k,1),1),kt(seg.ikseg(k,1),2),'ro');

                end
            end
            if(strcmp(spmode,'alt'))
                for k=1:2:acqP.necho
                    plot(kt(seg.ikseg(k,1):seg.ikseg(k,2),1),(kt(seg.ikseg(k,1):seg.ikseg(k,2),2)),col(1,:));
                    plot(kt(seg.ikseg(k,2),1),kt(seg.ikseg(k,2),2),'bo');
                end
                col=flipud(col);
                for k=2:2:acqP.necho
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
    GS3p1=mr.makeExtendedTrapezoid('z','times',[0 GSex.fallTime],'amplitudes',[GSex.amplitude 0]);


    %GS3p3=mr.makeExtendedTrapezoid('z','times',[0 GSref.fallTime],'amplitudes',[GSref.amplitude 0]);
    GS4times=[0 GSref.flatTime];
    GS4amp=[GSref.amplitude GSref.amplitude];
    GS4 = mr.makeExtendedTrapezoid('z','times',GS4times,'amplitudes',GS4amp);
    %segP.TEprep= GSspr.flatTime;
    GS5times=[0 GSspr.riseTime GSspr.riseTime+GSspr.flatTime GSspr.riseTime+GSspr.flatTime+GSspr.fallTime];
    GS5amp=[GSref.amplitude GSspr.amplitude GSspr.amplitude 0];
    GS5 = mr.makeExtendedTrapezoid('z','times',GS5times,'amplitudes',GS5amp);
    %     GS5p1times=[0 GSspr.riseTime GSspr.riseTime+segP.TEprep GSspr.riseTime+segP.TEprep+GSspr.fallTime];
    %     GS5p1amp=[GSref.amplitude segP.GSfac*GSspr.amplitude segP.GSfac*GSspr.amplitude 0];

    GS5p1times=[0 dG];
    GS5p1amp=[GSref.amplitude 0];
    GS5p1 = mr.makeExtendedTrapezoid('z','times',GS5p1times,'amplitudes',GS5p1amp);

    %%
    GS7times=[0 GSspr.riseTime GSspr.riseTime+GSspr.flatTime GSspr.riseTime+GSspr.flatTime+GSspr.fallTime];
    GS7amp=[0 GSspr.amplitude GSspr.amplitude GSref.amplitude];
    GS7 = mr.makeExtendedTrapezoid('z','times',GS7times,'amplitudes',GS7amp);


    GS7p1times=[0 GSspr.riseTime GSspr.riseTime+segP.TEprep GSspr.riseTime+segP.TEprep+GSspr.fallTime];
    GS7p1amp=[0 segP.GSfac/2*GSspr.amplitude segP.GSfac/2*GSspr.amplitude GSref.amplitude];
    GS7p1= mr.makeExtendedTrapezoid('z','times',GS7p1times,'amplitudes',GS7p1amp);

    GSarea=calcArea(GS4)/2+calcArea(GS7);
    GSprep_area1=calcArea(GS2)/2+calcArea(GS3p1)+calcArea(GS7p1)+calcArea(GS4)/2;
    GSprep_area2=calcArea(GS4)/2+calcArea(GS5p1);
    GS3p2=mr.makeTrapezoid('z',system,'area',segP.GSfac*GSarea-GSprep_area1,'duration',segP.t1,'riseTime',dG);
    GS5p2=mr.makeTrapezoid('z',system,'area',segP.GSfac*GSarea-GSprep_area2,'duration',(acqP.TE-2*segP.tSp-segP.tRefwd)/2,'riseTime',dG);

    % and now the readout gradient....
    GRpre_s = mr.makeTrapezoid('x',system,'area',spiral.kStart,'duration',segP.t1,'riseTime',dG);
    GRpre = mr.makeTrapezoid('x',system,'area',spiral.kStart,'duration',segP.tSpex,'riseTime',dG);

    GRref = mr.makeTrapezoid('x',system,'area',spiral.kStart,'duration',acq.readoutTime,'riseTime',dG);
    GRtrim = mr.makeTrapezoid('x',system,'area',2*spiral.kStart,'duration',acqP.TE-GSref.flatTime-2*segP.tSp,'riseTime',dG);

    GRspoi_x=mr.makeTrapezoid('x',system,'area',(segP.GXfac-1)*spiral.kStart,'duration',GSspr.riseTime+segP.TEprep,'riseTime',dG);
    GRspoi_xs=mr.makeTrapezoid('x',system,'area',(segP.GXfac-1)*spiral.kStart,'duration',GSspr.riseTime+segP.TEprep,'riseTime',dG);
    GRspoi_y=mr.makeTrapezoid('y',system,'area',segP.GYfac*spiral.kStart,'duration',GSspr.riseTime+segP.TEprep,'riseTime',dG);

    % and filltimes
    segS.tEx=GS1.shape_dur+GS2.shape_dur+GS3.shape_dur;
    segS.tRef=GS4.shape_dur+GS5.shape_dur+GS7.shape_dur+acq.readoutTime;
    tend=GS4.shape_dur+GS5.shape_dur;
    tETrain=segS.tEx+acqP.necho*segS.tRef+tend;
    %TRfill=floor(1e5*(acqP.TR-acqP.NSlices*tETrain)/acqP.NSlices)*1e-5;
    TRfill=acqP.TR;
    if TRfill<0, TRfill=1e-3;
        disp(strcat('Warning!!! acqP.TR too short, adapted to include all slices to : ',num2str(1000*acqP.NSlices*(tETrain+TRfill)),' ms'));
    else
        disp(strcat('TRfill : ',num2str(1000*TRfill),' ms'));
    end
    delayTR = mr.makeDelay(TRfill);
    if(strcmp(T2prep,'on'))
        TEfill1=acqP.TEprep/2-GS2.shape_dur/2-GS3p1.shape_dur-segP.t1-GS7p1.shape_dur-GS4.shape_dur/2;
        TEfill2=segP.t2+acqP.TEprep/2+acqP.TE/2-(segP.tRef+GS5p1.shape_dur+segP.tSp+acqP.TE-2*segP.tSp-segP.tRefwd+GS7.shape_dur);
        delayTE1 = mr.makeDelay(TEfill1);
        delayTE2 = mr.makeDelay(TEfill2);
    end
    % and flip angles
    acqP.rflip=acqP.flipref+zeros([1 acqP.necho]);
    if(acqP.flipflag==1),  acqP.rflip(1)=90+acqP.flipref/2;end
    if(acqP.flipflag==2)
        [rf,~] = fliptraps(acqP.flipref,(acqP.nrep+1)*acqP.necho,6,'opt',0,2,0,acqP.flipref,acqP.flipref,[6 5 5 acqP.necho]);
        %if(strcmp(T2prep,'off')), rf=[180 rf]; end
        %rf=[180 rf];
        rf=rf(1:acqP.nrep*acqP.necho); pow=sum(rf.^2)/sum((0*rf+180).^2);
        disp(strcat('rel. power :',num2str(pow)));
        acqP.rflip=rf(1:acqP.nrep*acqP.necho);
    end
    if(strcmp(acqP.PEtype,'linear')), acqP.PEind=acqP.necho-(mod(acqP.nTE-1+(acqP.necho-[1:acqP.necho]),acqP.necho)); end
    if(strcmp(acqP.PEtype,'centric'))
        acqP.PEind=zeros([1 acqP.necho]);
        acqP.PEind(1:acqP.nTE-1)=[2*acqP.nTE-2:-2:2];
        acqP.PEind(acqP.nTE:2*acqP.nTE-1)=[1:2:2*acqP.nTE-1];
        acqP.PEind(2*acqP.nTE:end)=[2*acqP.nTE:acqP.necho];
    end

    %% Define sequence blocks
    % Next, the blocks are put together to form the sequence
    if(~exist('scanmode')), scanmode='run'; end
    ntr=0;
    switch scanmode
        case 'init', dKA=zeros([acqP.necho 2]); dKE=zeros([acqP.necho 2]); % dKA and dKE are trim values to correct for trajectory imperfections
        case 'trim',  decho=4; acqP.sliceGAP=0;  ntrim=[1 2:decho:acqP.necho]; acqP.NSlices=length(ntrim); acqP.dummy=1;
        case 'allTE',  acqP.sliceGAP=0; acqP.NSlices=acqP.necho;
    end
    nk=acqP.necho;

    for kex=acqP.dummy:acqP.nex
        for s=1:acqP.NSlices
            dz=acqP.sliceGAP*acqP.sliceThickness*(acqP.Oslices(s));
            rfex.freqOffset=acqP.sliceGAP*GSex.amplitude*acqP.sliceThickness*(acqP.Oslices(s));
            rfref.freqOffset=acqP.sliceGAP*GSref.amplitude*acqP.sliceThickness*(acqP.Oslices(s));
            rfex.phaseOffset=rfex_phase-2*pi*rfex.freqOffset*mr.calcRfCenter(rfex); % align the phase for off-center slices
            rfref.phaseOffset=rfref_phase-2*pi*rfref.freqOffset*mr.calcRfCenter(rfref); % dito
            rfref.signal=refenvelope;
            % phaseOffset is moved into echo loop


            if(strcmp(fatsat,'on')), seq.addBlock(rf_fs,gz_fs); end
            if(strcmp(fatsat,'no')), seq.addBlock(gz_fs); end

            if(strcmp(T2prep,'on'))
                seq.addBlock(GS1);
                seq.addBlock(GS2,rfex);
                seq.addBlock(GS3p1);
                seq.addBlock(GS3p2,GRpre_s);
                seq.addBlock(delayTE1);
                seq.addBlock(GS7p1,GRspoi_x,GRspoi_y);
                seq.addBlock(GS4,rfref);
                %seq.addBlock(GS5p1,GRspoi_x,GRspoi_y);
                %seq.addBlock(GS5p1,GRspoi_y);
                seq.addBlock(GS5p1);
                seq.addBlock(delayTE2);
                %seq.addBlock(GRref);
                %              seq.addBlock(GS7);
            else
                seq.addBlock(GS1);
                seq.addBlock(GS2,rfex);
                seq.addBlock(GS3,GRpre);
            end
            if(ntr>0), nk=s; end

            %%
            %for m=1:1
            for m=1:acqP.nrep*acqP.necho

                if(strcmp(scanmode,'allTE')),acqP.nTE=s;end
                %k=nk-(mod(acqP.nTE-1+(nk-m),nk))
                indk=mod(m,acqP.necho);
                if(indk==0), indk=acqP.necho; end
                k=acqP.PEind(indk);
                % calculate gradients
                if((k==1)&&strcmp(initmode,'rev'))
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
                        [GBegin, tBegin, tt] = spiral_k2k_opt(kBegin1, kEnd1, GStart, GDest, system.maxGrad, system.maxSlew, segS.tSp, system.gradRasterTime);
                    end
                else
                    %% initial gradient
                    segS.tSp=GSspr.riseTime+GSspr.flatTime+GSspr.fallTime+seg.i1(k)*system.gradRasterTime;
                    if((m==1)&&strcmp(T2prep,'on')), kBegin1=[-segP.GXfac*spiral.kStart+dKA(k,1) -segP.GYfac*spiral.kStart];
                    else
                        kBegin1=[-spiral.kStart+dKA(k,1) 0];
                    end
                    kEnd1=[kt(seg.ikseg(k,1),1) kt(seg.ikseg(k,1),2)+dKA(k,2)];
                    GStart=[0 0];
                    if(strcmp(spmode,'cont')||(k~=2*floor(k/2))),
                        GDest=[Gsp(seg.ikseg(k,1),1) Gsp(seg.ikseg(k,1),2)];
                    end
                    if(strcmp(spmode,'alt')&&(k==2*floor(k/2)))
                        GDest=-[Gsp(seg.ikseg(k,1),1) Gsp(seg.ikseg(k,1),2)];
                    end
                    if(strcmp(segmode,'tan'))
                        [GBegin, tBegin, Ginit, tinit] = spiral_k2k_m(kBegin1, kEnd1, GStart, GDest, segS.tSp, system.maxSlew, system.gradRasterTime);
                    else
                        if((m==1)&&strcmp(T2prep,'on')),[GBegin, tBegin,Gtran, tG, tt] = spiral_k2k_opt(kBegin1, kEnd1, GStart, GDest, system.maxGrad, system.maxSlew, segS.tSp-segP.t2, system.gradRasterTime);
                        else
                            [GBegin, tBegin,Gtran, tG, tt] = spiral_k2k_opt(kBegin1, kEnd1, GStart, GDest, system.maxGrad, system.maxSlew, segS.tSp, system.gradRasterTime);
                        end
                    end

                end
                %%
                % terminal gradient
                ntk=acq.ntot-seg.i1(k)-seg.nkseg(k);
                segS.tSp=GSspr.riseTime+GSspr.flatTime+GSspr.fallTime+ntk*system.gradRasterTime;
                kBegin2=[kt(seg.ikseg(k,2),1) kt(seg.ikseg(k,2),2)+dKE(k,2)];
                kEnd2=[spiral.kStart+dKE(k,1) 0];
                GDest=[0 0];
                if(strcmp(spmode,'cont')||(k~=2*floor(k/2))),
                    GStart=[Gsp(seg.ikseg(k,2),1) Gsp(seg.ikseg(k,2),2)];
                end
                if(strcmp(spmode,'alt')&&(k==2*floor(k/2)))
                    GStart=-[Gsp(seg.ikseg(k,2),1) Gsp(seg.ikseg(k,2),2)];
                end
                if(strcmp(segmode,'tan'))
                    [GEnd, tEnd, Gfin, tfin] = spiral_k2k_m(kBegin2, kEnd2, GStart, GDest, segS.tSp, system.maxSlew, system.gradRasterTime);
                else
                    [GEnd, tEnd, Gfin, tfin] = spiral_k2k_opt(kBegin2, kEnd2, GStart, GDest, system.maxGrad, system.maxSlew, segS.tSp, system.gradRasterTime);

                end
                % combine gradients


                if((k==1)&&strcmp(initmode,'rev'))
                    Gtot=[GBegin(1:end,:); -Gsp(nkseg1:-1:seg.ikseg(k,1),:); Gsp(seg.ikseg(k,1)+1:seg.ikseg(k,2),:); GEnd(1:end,:);[0 0]];
                else
                    if((k==1)&&strcmp(initmode,'outin'))
                        Gtot(:,1)=[GEnd(end:-1:1,1); Gsp(seg.ikseg(k,2):-1:seg.ikseg(k,1),1); Gsp(seg.ikseg(k,1)+1:seg.ikseg(k,2),1); GEnd(1:end,1);[0]];
                        Gtot(:,2)=[GEnd(end:-1:1,2); Gsp(seg.ikseg(k,2):-1:seg.ikseg(k,1),2); Gsp(seg.ikseg(k,1)+1:seg.ikseg(k,2),2); GEnd(1:end,2);[0]];
                    else

                        if((k==1)&&strcmp(initmode,'sym'))
                            Gtot(:,1)=[GEnd(end:-1:1,1); Gsp(seg.ikseg(k,2):-1:seg.ikseg(k,1),1); Gsp(seg.ikseg(k,1)+1:seg.ikseg(k,2),1); GEnd(1:end,1);0];
                            Gtot(:,2)=[-GEnd(:,2); -Gsp(seg.ikseg(k,2):-1:seg.ikseg(k,1),2); Gsp(seg.ikseg(k,1)+1:seg.ikseg(k,2),2); GEnd(1:end,2);[0]];
                        else

                            if(strcmp(spmode,'cont')||(k~=2*floor(k/2)))
                                Gtot=[GBegin(1:end,:); Gsp(seg.ikseg(k,1)+1:seg.ikseg(k,2),:); GEnd(1:end,:);[0 0]];
                            end
                            if(strcmp(spmode,'alt')&&(k==2*floor(k/2)))
                                Gtot=[GBegin(1:end,:); -Gsp(seg.ikseg(k,1):-1:seg.ikseg(k,2)+1,:); GEnd(1:end,:);[0 0]];
                            end

                            ttot=[tBegin(1:end) [(1:seg.nkseg(k))*system.gradRasterTime] tEnd(2:end)];

                        end
                    end
                end

                nGtot=length(Gtot);
                ttot=[1:nGtot]*system.gradRasterTime;




                % split Gtot
                if((m==1)&&strcmp(T2prep,'on')),
                    nsp=round((GSspr.riseTime+GSspr.flatTime+GSspr.fallTime-segP.t2)/system.gradRasterTime);
                else
                    nsp=round((GSspr.riseTime+GSspr.flatTime+GSspr.fallTime)/system.gradRasterTime);
                end

                Gsp1=Gtot(1:nsp,:);
                Gspiral=Gtot(nsp+1:nsp+acq.ntot,:);
                Gsp2=Gtot(nsp+acq.ntot+1:end,:);

                Gsp1_x=mr.makeArbitraryGrad('x','system',system,'waveform',Gsp1(:,1),'first',0,'last',mean(Gtot(nsp:nsp+1,1)));
                Gsp1_y=mr.makeArbitraryGrad('y','system',system,'waveform',Gsp1(:,2),'first',0,'last',mean(Gtot(nsp:nsp+1,2)));

                if(strcmp(segmode,'single'))
                    for kadc=1:1
                        Gspiral_x(kadc)=mr.makeArbitraryGrad('x','system',system,'maxSlew',4.2576e+11,'waveform',Gspiral(:,1));
                        Gspiral_y(kadc)=mr.makeArbitraryGrad('y','system',system,'maxSlew',4.2576e+11,'waveform',Gspiral(:,2));
                    end

                else
                    Gspiral_x=mr.makeArbitraryGrad('x','system',system,'maxSlew',4.2576e+11,'waveform',Gspiral(:,1),'first',Gsp1_x.last,'last',mean(Gtot(nsp+acq.ntot:nsp+acq.ntot+1,1)));
                    Gspiral_y=mr.makeArbitraryGrad('y','system',system,'maxSlew',4.2576e+11,'waveform',Gspiral(:,2),'first',Gsp1_y.last,'last',mean(Gtot(nsp+acq.ntot:nsp+acq.ntot+1,2)));
                end

                Gsp2_x=mr.makeArbitraryGrad('x','system',system,'waveform',Gsp2(:,1),'first',Gspiral_x.last,'last',0);
                Gsp2_y=mr.makeArbitraryGrad('y','system',system,'waveform',Gsp2(:,2),'first',Gspiral_y.last,'last',0);
                rfref.signal=refenvelope*acqP.rflip(m)/180;
                dOffset(indk)=0;
                if(kconc==2),   dOffset(indk)=-acq.concB(indk).*dz^2; end
                rfref.phaseOffset=dOffset(indk) + rfref_phase-2*pi*rfref.freqOffset*mr.calcRfCenter(rfref);
                if (m>1)
                    seq.addBlock(GS4,rfref);
                else
                    if(strcmp(T2prep,'off'))
                        seq.addBlock(GS4,rfref);
                    end
                end
                if(strcmp(scanmode,'trim')&&(k==ntrim(s)))
                    seq.addBlock(GS5);
                    seq.addBlock(GRtrim,adc);
                    seq.addBlock(GS7);
                else
                    if (Gsp1_x.shape_dur>GS5.shape_dur)
                        fprintf('\ntrajectory not compatible with sequence parameters. \nChange parameters and/or buy faster/stronger gradients\n');
                        return
                    end
                    Gsp1_x.first=0;
                    if(m>1)
                        seq.addBlock(Gsp1_x,Gsp1_y,GS5);
                    else
                        if(strcmp(T2prep,'off'))
                            seq.addBlock(Gsp1_x,Gsp1_y,GS5);
                        else
                            seq.addBlock(Gsp1_x,Gsp1_y);
                        end
                    end


                    if (kex>0)
                        if ((m==1)&& strcmp(T2prep,'on'))
                            seq.addBlock(Gspiral_x,Gspiral_y,GS5p2,adc);
                        else
                            seq.addBlock(Gspiral_x,Gspiral_y,adc);
                        end
                    else
                        seq.addBlock(Gspiral_x,Gspiral_y);
                    end
                    if(Gsp2_x.shape_dur>GS7.shape_dur)
                        fprintf('\ntrajectory not compatible with sequence parameters. \nChange parameters and/or buy faster/stronger gradients\n');
                        return
                    end
                    seq.addBlock(Gsp2_x,Gsp2_y,GS7);
                end


            end
            %%
            seq.addBlock(GS4);
            %[duration, numBlocks, eventCount]=seq.duration();
            seq.addBlock(GS5);
            seq.addBlock(delayTR);

        end
    end
    %end



    %% calculate concomitant fields
    [ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();
    %gw=seq.gradient_waveforms();
    %  wave_data=seq.waveforms_and_times();
    %         %...
    %         GRT=system.gradRasterTime;
    %         n_excitation=round(t_excitation/GRT);
    %         n_refocusing=round(t_refocusing/GRT);
    %         k=1;
    %         gx=gw(1,n_refocusing(1)-acqP.TE/GRT/2+1:n_refocusing(1))/gamma;
    %         gy=gw(2,n_refocusing(1)-acqP.TE/GRT/2+1:n_refocusing(1))/gamma;
    %         concB=(gx.^2+gy.^2)./(2*B0);
    %         concBz(k)=sum(concB).*gamma*GRT*2*pi;
    %         figure
    %         hold on
    %         plot(gw(1,1:20000));
    %         plot(n_refocusing(1)-acqP.TE/GRT/2+1:n_refocusing(1),gx,'r.')
    %         %
    %         for k=2:acqP.necho,
    %             gx=gw(1,n_refocusing(k-1):n_refocusing(k))/gamma;
    %             gy=gw(2,n_refocusing(k-1):n_refocusing(k))/gamma;
    %             concB=(gx.^2+gy.^2)./(2*B0);
    %             concBz(k)=sum(concB).*gamma*system.gradRasterTime*2*pi;
    %             if(floor(k/2)==k/2)
    %                 plot(n_refocusing(k-1):n_refocusing(k),gx,'g.');
    %             else
    %             plot(n_refocusing(k-1):n_refocusing(k),gx,'r.');
    %             end
    %
    %
    %         end,
    %         if(plotflag(6)=='0'), close; end;
    %         acq.concB=concBz(1:acqP.necho);




    %% plot k-spaces

    % figure; plot(ktraj'); % plot the entire k-space trajectory
    if(plotflag(4)=='1')
        figure; plot(ktraj(1,:),ktraj(2,:),'b',...
            ktraj_adc(1,:),ktraj_adc(2,:),'r'); % a 2D plot
        set(gcf,'Position',[100 100 550 350]);
        axis([-1.2*spiral.kStart 1.2*spiral.kStart -1.2*spiral.kStart 1.2*spiral.kStart]);
        axis equal
    end
    if(plotflag(5)=='1')
        figure; plot(ktraj(1,:),ktraj(3,:),'b',...
            ktraj_adc(1,:),ktraj_adc(3,:),'r'); % a 2D plot
        set(gcf,'Position',[500 100 550 350]);
        %axis equal
        axis([-1.2*spiral.kStart 1.2*spiral.kStart -1.2*spiral.kStart 1.2*spiral.kStart]);
    end

    %% Write to file

    % The sequence is written to file in compressed form according to the file
    % format specification using the |write| method.
    %pn=which('myspiralTSE');
    if(~exist('seqno')), seqno=1; end
    seqname=strcat('ssTSE_',num2str(seqno));
    %seqname=allname;
    seq.write(strcat(seqname,'.seq'))
    save(strcat('p',seqname),'segmode','spmode','scanmode','accmode','fatsat','T2prep','acq','acqP','seg','segP','spiral','kt','ktraj_adc','system','B0');
    seqno=seqno+1;
    %%
    % Display the first few lines of the output file
    s=fileread(strcat(seqname,'.seq'));
    disp(s(1:300))
    %% plot gradients etc
    if(plotflag(6)=='1')
        fig=seq.plot('TimeRange',[0 0.25],'timeDisp','ms');
        %fig=seq.plot('TimeRange',[0 acqP.TEprep+acqP.TE+acq.readoutTime+20e-3],'timeDisp','ms');
    end
    %
    
end
toc
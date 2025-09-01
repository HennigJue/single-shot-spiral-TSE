if(kseq==1), acqP.T2prep='on'; acqP.navmode=[1 1]; acqP.PEtype='centric'; acqP.TEeff=40e-3; acqP.TEprep=40e-3; acq.GSSEfac=0.2; end
if(kseq==2), acqP.T2prep='on'; acqP.navmode=[1 1]; acqP.PEtype='centric'; acqP.TEeff=40e-3; acqP.TEprep=40e-3; acqP.tD=10e-3; acqP.GDs=30; acq.GSSEfac=0.2;    end
if(kseq==3), acqP.T2prep='on'; acqP.navmode=[1 1]; acqP.PEtype='centric'; acqP.TEeff=40e-3; acqP.TEprep=40e-3; acqP.tD=10e-3; acqP.GDx=30; acq.GSSEfac=0.2;   end

if(kseq==4), acqP.Nx=200; acqP.Ny=100; acqP.necho=100;acqP.samplingTime= 2e-3; acqP.TE=4e-3; acq.fspR=0.2;acq.fspS=0.2;
    acqP.NSlices=1; acqP.TR=10;acqP.T2prep='off'; acqP.navmode=[1 1]; acqP.PEtype='linear'; 
    acqP.TEeff=40e-3; acqP.flipref=40; acqP.flipflag=2; acqP.HF_fac=0.875; acqP.PEref=16; 
    acq.tEx=0.8e-3; acq.tBwPex=1;acq.tRef=0.8e-3;acq.tBwPref=1; 
end
if(kseq==5), acqP.Nx=200; acqP.Ny=100; acqP.necho=100;acqP.samplingTime= 2e-3;  acqP.TE=4.8e-3; acq.fspR=0.5;acq.fspS=0.5;
    acqP.NSlices=1; acqP.TR=10;acqP.T2prep='on'; acqP.navmode=[1 1]; acqP.PEtype='centric'; 
    acqP.TEeff=40e-3; acqP.tD=0e-3; acqP.GDs=0; acqP.flipref=40; acqP.flipflag=2; acqP.HF_fac=0; acqP.PEref=16; 
    acq.tEx=0.8e-3; acq.tBwPex=1;acq.tRef=0.8e-3;acq.tBwPref=1;acq.tRefSE=1.2e-3;acq.tBwPrefSE=1;acq.GSSEfac=0.2;
end
if(kseq==6), acqP.Nx=200; acqP.Ny=100; acqP.necho=100;acqP.samplingTime= 2e-3; acqP.TE=4.8e-3; acq.fspR=0.5;acq.fspS=0.5;
    acqP.NSlices=1; acqP.TR=10;acqP.T2prep='on'; acqP.navmode=[1 1]; acqP.PEtype='centric'; 
    acqP.TEeff=40e-3; acqP.tD=10e-3; acqP.GDs=20; acqP.flipref=40; acqP.flipflag=2; acqP.HF_fac=0; acqP.PEref=16; 
    acq.tEx=0.8e-3; acq.tBwPex=1;acq.tRef=0.8e-3;acq.tBwPref=1;acq.tRefSE=1.2e-3;acq.tBwPrefSE=1;acq.GSSEfac=0.2;
end
if(kseq==7), acqP.Nx=200; acqP.Ny=100; acqP.necho=100;acqP.samplingTime= 2e-3; acqP.TE=4.8e-3; acq.fspR=0.5;acq.fspS=0.5;
    acqP.NSlices=1; acqP.TR=10;acqP.T2prep='on'; acqP.navmode=[1 1]; acqP.PEtype='centric'; 
    acqP.TEeff=40e-3; acqP.tD=10e-3; acqP.GDx=20; acqP.flipref=40; acqP.flipflag=2; acqP.HF_fac=0; acqP.PEref=16; 
    acq.tEx=0.8e-3; acq.tBwPex=1;acq.tRef=0.8e-3;acq.tBwPref=1;acq.tRefSE=1.2e-3;acq.tBwPrefSE=1;acq.GSSEfac=0.2;
end

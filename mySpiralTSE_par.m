
switch seqvar_mod
    case 'none'
        seqname=strcat('TSE_',T2prep,'_',num2str(round(1000*acqP.TE)),'_',num2str(floor(1000*acqP.TEeff)),'_',num2str(100*acqP.sliceGAP),'_',num2str(acqP.NSlices));
    case 'Diff'
        Diffvals=[4 10 20 30 40 4 10 20 30 40];
        if(kseq)>10, quitflag=1; return; end
        segP.GXfac=Diffvals(kseq);
        segP.GYfac=0;
        if(kseq>5),
            segP.GYfac=Diffvals(kseq);
            segP.GXfac=0;
        end
        seqname=strcat('TSE_',T2prep,'_',num2str(round(1000*acqP.TE)),'_',num2str(floor(1000*acqP.TEeff)),'_',num2str(100*acqP.sliceGAP),'_',num2str(acqP.NSlices),'_Diff');
        allname=strcat('TSEall_',num2str(round(1000*acqP.TE)),'_',num2str(100*acqP.sliceGAP),'_',num2str(acqP.NSlices),'_Diff');
    case 'TErep'
        seqname=strcat('TSE_',T2prep,'_',num2str(round(1000*acqP.TE)),'_',num2str(floor(1000*acqP.TEeff)),'_',num2str(100*acqP.sliceGAP),'_',num2str(acqP.NSlices),'_',num2str(acqP.nrep));
        
    case 'T1var'
        seqname=strcat('TSE_T1var_',T2prep,'_',num2str(round(1000*acqP.TE)),'_',num2str(floor(1000*acqP.TEeff)),'_',num2str(100*acqP.sliceGAP),'_',num2str(acqP.NSlices));
    case 'MTvar'
        seqname=strcat('TSE_MTvar_',T2prep,'_',num2str(acqP.MTfac),'_',num2str(round(1000*acqP.TE)),'_',num2str(floor(1000*acqP.TEeff)),'_',num2str(100*acqP.sliceGAP),'_',num2str(acqP.NSlices));
    case 'TEvar'
        TEvals=[1 1 2 3 4 7 10 13 16];
        if(kseq)>14, quitflag=1; return; end
        acqP.TEeff=TEvals(kseq)*acqP.TE;
        acqP.TEprep=acqP.TEeff;
        
        
        %         if(kseq<7), T2prep='on'; acqP.nTE=1;
        %         else T2prep='off'; acqP.nTE=round(acqP.TEeff/acqP.TE); end
        seqname=strcat('TSE_',T2prep,'_',num2str(round(1000*acqP.TE)),'_',num2str(floor(1000*acqP.TEeff)),'_',num2str(100*acqP.sliceGAP),'_',num2str(acqP.NSlices));
        allname=strcat('TSEvar_',num2str(round(1000*acqP.TE)),'_',num2str(100*acqP.sliceGAP),'_',num2str(acqP.NSlices));
    case 'TEall'
        TEvals=[10:10:220];
        if(kseq)>22, quitflag=1; return; end
        acqP.TEeff=TEvals(kseq)*10^-3;
        acqP.TEprep=acqP.TEeff;
        T2prep='on';
        T2prep='off'; acqP.nTE=round(acqP.TEeff/acqP.TE);
        seqname=strcat('TSE_',T2prep,'_',num2str(round(1000*acqP.TE)),'_',num2str(floor(1000*acqP.TEeff)),'_',num2str(100*acqP.sliceGAP),'_',num2str(acqP.NSlices));
        allname=strcat('TSEall_',num2str(round(1000*acqP.TEeff)),'_',num2str(100*acqP.sliceGAP),'_',num2str(acqP.NSlices));
        
    case 'accvar'
        if(kseq)>8, quitflag=1; return; end
        N_var=[0.7 1];
        acc_var=[1.2 1.4 2 4];
        if(2*floor(kseq/2)==kseq),
            spiral.N=N_var(2);
        else
            spiral.N=N_var(1);
        end
        seg.n_acc=acc_var(floor(kseq/2));
        seqname=strcat(seqtype,num2str(round(1000*acqP.TE)),'_',num2str(floor(1000*acqP.TEeff)),'_',num2str(acqP.flipref),'_',num2str(10*seg.n_acc),'_',num2str(10*spiral.N),'_',num2str(acqP.NSlices));
    case 'fsvar',
        fatsat='on'
        if(kseq)>7, quitflag=1; return; end
        fac=[1 1 1.1 1.2 1.3 1.4 1.5];
        acqP.sat_freq=fac(kseq)*acqP.sat_ppm*1e-6*B0*system.gamma;
        seqname=strcat(seqtype,num2str(round(1000*acqP.TE)),'_',num2str(floor(1000*acqP.TEeff)),'_',num2str(acqP.flipref),'_',num2str(floor(abs(acqP.sat_freq))),'_',num2str(acqP.NSlices));
        
end
%return



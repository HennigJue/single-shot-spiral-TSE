function [PEorder,PEindex,NS,IPS] = PE_order(PEtype,NY,Necho,k0,PE_steps,HFfac,Nref,GRfac)
if(isempty(PE_steps))
    if(strcmp(PEtype,'faise')),
        disp('PE-mode changed to centric encoding');
        PEtype='centric';
    end
    PE_first=floor(1+HFfac*NY/2);
    PE_ref=floor(NY/2-Nref/2);
    GRPElow=[PE_ref:-GRfac:PE_first];
    PE_steps=[GRPElow(end:-1:1) PE_ref+1:PE_ref+Nref-1 PE_ref+Nref:GRfac:NY];
    nPE=length(PE_steps);
    NS=floor(nPE/Necho);
    if(NS==0), NS=1; Necho=length(PE_steps);

    else
        PE_steps=PE_steps(1:NS*Necho);
    end
end
nPE=length(PE_steps);
NS=floor(nPE/Necho);
if(NS==0), NS=1; Necho=length(PE_steps);
else
    PE_steps=PE_steps(1:NS*Necho);
end
ind=find(PE_steps<NY/2);
nPE0=ind(end)+1;
PEindex=PE_steps;


switch PEtype
    case 'linear'
        %PE0=ind(end)+1;
        PEorder=reshape(PE_steps,[NS Necho])'-NY/2;

        [ke,ks]=find(PEorder==0);
        PEorder=circshift(PEorder,[k0-ke(1) 0]);


    case {'centric'}
        if(strcmp(PEtype,'centric'))
            if(NS/2-floor(NS/2)>0)
                disp('Warning! odd number of excitations may be not optimal for centric mode!')
            end
        end
        ind=0;
        kn=0*PE_steps;
        kn(1)=ind(1);
        for k=1:floor(length(PE_steps)/2)
            kn(2*(k-1)+2)=ind+k;
            kn(2*(k-1)+3)=ind-k;
        end
        NS=floor(length(PE_steps)/Necho);
        kn=kn(1:NS*Necho);
        kn=kn-min(kn)+1;
        knorder=reshape(kn,[NS Necho]);
        [ke,ks]=find(knorder'==nPE0);
        knorder=circshift(knorder',[-ke+k0 0]);
        PEorder=PE_steps(knorder)-NY/2;

        
    case 'faise'
        if(floor(NS/2)==NS/2)
            kn=0*PE_steps;
            ind=0;
            kn(1)=ind(1);
            for k=1:floor(length(PE_steps)/2)
                kn(2*(k-1)+2)=ind+k;
                kn(2*(k-1)+3)=ind-k;
            end
            NS=floor(length(PE_steps)/Necho);
            kn=kn(1:NS*Necho);
            kn=kn-min(kn)+1;
            knorder=reshape(kn,[NS Necho]);
            [ke,ks]=find(knorder'==nPE0);
            knorder=knorder+NS/2*(k0-ke);
            knorder=mod(knorder-1,length(PE_steps))+1;
            PEorder=PE_steps(knorder);
            PEorder=PEorder'-NY/2;
            [ke,ks]=find(knorder'==nPE0);

        else
            NYF=nPE;
            IPS=round(NS/4)+round((k0-1)*NS/2);
            for m=1:floor(NS/2)
                for n=1:(NYF/NS)
                    A(m,n)=mod(NYF/2-n*NS/2+m-1+IPS,NYF)+1;
                end
            end
            for m=floor(NS/2)+1:NS
                for n=1:(NYF/NS)
                    A(m,n)=mod(NYF/2+(n-2)*NS/2+m-1+IPS,NYF)+1;
                end
            end
            for m=1:NS
                for n=1:floor(NYF/NS)
                    k(m,n)=(A(m,n)-(NYF+1)/2);
                end
            end
            ks=k(:);
            [kss,indk]=sort(ks-min(ks)+1);
            PE_new(indk)=PE_steps;
            PE_new1=PE_new-PE_new(1)+ks(1);
            PEorder=reshape(PE_new1,[NS Necho]);
            PEorder=PEorder';
        end

    case 'paired'
        ind=0;
        kn=0*PE_steps;
        kn(1)=ind(1);
        for k=1:floor(length(PE_steps)/4)
            kn(4*(k-1)+2)=ind+2*(k-1)+1;
            kn(4*(k-1)+3)=ind+2*(k-1)+2;
            kn(4*(k-1)+4)=ind-2*(k-1)-1;
            kn(4*(k-1)+5)=ind-2*(k-1)-2;
        end
        kn=kn(1:length(PE_steps));
        [ke,ks]=find(kn==nPE0);
        knorder=circshift(kn,[-ke+k0 0])-min(kn)+1;
        PEorder=PE_steps(knorder)-NY/2;




end
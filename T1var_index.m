function [Sindex] = T1var_index(NS,TRfac),
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%Sindex=zeros([2*NS length(TRfac)]);
%Sindex=[1:NS];
Sindex=[];
for nTR=1:length(TRfac)
    %for nfac=1:TRfac(nTR)
        temp=[1:TRfac(nTR):NS];
        %tempr=temp;
        Sindex=[Sindex temp];
        for nrep=2:TRfac(nTR)
            tempr=temp+nrep-1;
            Sindex=[Sindex tempr];
        end
    %end
   %Sindex(1:length(temp),nTR)=temp;
end
%Sindex=[1:TRfac(nTR):NS Sindex];
ind=find(Sindex>NS);
Sindex(ind)=NS;
%Sindex=temp;

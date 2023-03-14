function [Sindex] = MTvar_index(sltemp,MTfac,Sind),

if(nargin<3),Sind=floor(length(sltemp)/2); end;
maxS=ceil(MTfac(end)/length(sltemp));
Sindex=repmat(sltemp,[1 maxS]);
ind=find(Sindex==Sind);
Sindex(ind)=length(sltemp)+1;
Sindex(MTfac)=Sind;



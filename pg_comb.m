% Y= pg_comb(mpg)
%   conversion between phasegraph-formats:
%   
%   typ 'm2p':
%	used in pg_cpmg to sort phasestates mpg (from pg_cpmg) into Y as a representation 
%   suitable for display as phasegraph:
%   Y(sm...1,:,1) will contain all dephasing transversal states mpg(2,:,1...sm),
%   Y(sm+1...2 sm,:,1) will contain all rephasing transversal states mpg(1,:,1...sm).
%   Y(:,:,2) same for longitudinal states.
%
%   typ 'p2m':
%   to convert phasegraphs generated by phasegraph.m into phasestates suitable for use in pg_cpmg
function Y= pg_comb(mpg,typ)
if (nargin<2),
    typ='m2p';
end;

sm=size(mpg);
switch typ,
case 'm2p',
    
Y=zeros([2*sm(2) sm(3) 2]);
Y((sm(2)+1):2*sm(2),:,1)=squeeze(mpg(1,:,:));
Y(sm(2):-1:1,:,1)=mpg(1,:,:);
Y((sm(2)+1):2*sm(2),:,2)=squeeze(-mpg(3,:,:));
Y(sm(2):-1:1,:,2)=mpg(3,:,:);


case 'p2m',
smm=sm(1)/2;    
Y=zeros([4 smm sm(3)]);
size(Y);
Y(1,1:smm-1,:)=squeeze(mpg(smm-1:-1:1,1,:));
Y(2,1:smm+1,:)=squeeze(mpg(smm:2*smm,1,:));
Y(3,1:smm,:)=squeeze(mpg(smm:-1:1,2,:));
Y(4,1:smm,:)=squeeze(mpg(smm+1:2*smm,2,:));

end;
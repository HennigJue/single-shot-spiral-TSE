% function rf = fliptraps(rflip,nflip,ca);
% calculates flip-angles for smooth echo amplitudes
%
%	requires:
%	rflip:      target flip angle of refocusing pulses at the beginning of the echo train 
%   nETL:       echo train length
%   lZero:      number of echo for hyperecho (only for hyp=1)
%   ftyp:       type for flip angle generation (90+a/2,opt, lin, sin)
%   hyp:        flag for hyperecho
%   traps:      flag for TRAPS(1 = linear flip angles, 2 = smooth flip angles)
%   fb:         no
%   tr_high:    flip angle at plateau
%   tr_end:     flip angle at end
%   tr_var:     echo number for [start to plateau; begin plateau; end plateau; reach tr_end]
%   example:
%   [rf,pow] = fliptraps_new(60,20,10,'opt',0,2,0,90,40,[7 10 15 20]);
function [rf,pow] = fliptraps(rflip,nETL,lZero,ftyp,hyp,traps,fb,tr_high,tr_end,tr_var);

rf=ones([1 nETL])*rflip;
switch ftyp,
    
case '90+a/2',
rf(1)=90+rf(1)/2;

case 'opt',
ex1=2; ca=0.4; cb=2; cc=2;
rf(1)=90+rflip/2+ca*((cc-1)/cc)^ex1*(90-rflip/2);
dr=rf(1)-rflip;
    for k=2:nETL,
  % rf(1,k)=rflip+(cb^(nflip+1-k)-1)*(rf(1,1)-rflip)/cb^nflip;
        rf(k)=rflip+dr/cb^(k-0.5);
    end;

case 'lin',
   dinc=(180-rflip)/(tr_var(1)+1)
    for k=1:tr_var(1),
        rf(k)=180-k*dinc;
    end;

case 'sin',
  for k=1:tr_var(1),        
        rf(k)=180-(180-rflip)*sin(pi*k/(2*(tr_var(1)-1)));
  end; 

end;



if (traps==1),
    dinc=(tr_high-rf(tr_var(1)))/(tr_var(2)-tr_var(1));
    for k=tr_var(1)+1:tr_var(2),
        if (k>1),        rf(k)=rf(k-1)+dinc;    end;
    end;
    
    for k=tr_var(2)+1:tr_var(3),
        rf(k)=rf(k-1);
    end;
    dinc=(tr_end-rf(tr_var(3)))/(tr_var(4)-tr_var(3));
    for k=tr_var(3)+1:tr_var(4),
        rf(k)=rf(k-1)+dinc;
    end;
    for k=tr_var(4)+1:nETL,
        rf(k)=rf(k-1);
    end;
end;

if (traps==2),
    dinc=(tr_high-rf(tr_var(1)))/(tr_var(2)-tr_var(1));
    for k=tr_var(1)+1:tr_var(2),
        if (k>1),        
        rf(k)=tr_high-(tr_high-rf(tr_var(1)))*cos(pi*(k-tr_var(1))/(2*(tr_var(2)-tr_var(1))));
        end;
    end;
    rf(tr_var(2))=(tr_high+rf(tr_var(2)-1))/2;
    for k=tr_var(2)+1:tr_var(3),
        rf(k)=rf(k-1);
    end;
    dinc=(tr_end-rf(tr_var(3)))/(tr_var(4)-tr_var(3));
    for k=tr_var(3)+1:tr_var(4),
        rf(k)=tr_high-(tr_high-tr_end)*cos(pi*(k-tr_var(4))/(2*(tr_var(4)-tr_var(3))));
    end;
    for k=tr_var(4)+1:nETL,
        rf(k)=rf(k-1);
    end;
end;


if (hyp==1) &(traps==0), 
k1=floor(lZero/2);
if (k1==lZero/2), k1=k1-1; end;
for k=nETL:-1:(2*k1+2),
            rf(k)=-rf(k-2*k1-1);
end;
rf(k1+1)=180;
for k=(k1+2):(2*k1+1),
    rf(k)=-rf(2*(k1+1)-k);
end;
if (k1+1==lZero/2),
    for k=nETL:-1:2,
        rf(k)=rf(k-1);
    end
    rf(1)=180;
end;
end;

if (hyp==1) &(traps==1), 

k1=lZero;

for k=nETL:-1:(2*k1+1),
            rf(k)=rf(k-2*k1-1);
end;
rf(k1)=180;
for k=(k1+2):(2*k1+1),
    rf(k)=rf(2*(k1+1)-k);
end;

end;
pow=sum((rf./180).^2)/nETL;





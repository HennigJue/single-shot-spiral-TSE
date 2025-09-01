% function rf = fliptraps(rflip,nflip,ca);
% calculates flip-angles for smooth echo amplitudes
%
%	requires:
%	rflip:      target flip angle of refocusing pulses at the beginning of the echo train 
%   nETL:       echo train length
%   ftyp:       type for flip angle generation (90+a/2,opt, lin, sin)
%   traps:      flag for TRAPS(1 = linear flip angles, 2 = smooth flip angles)
%   rf_high:    flip angle at plateau
%   rf_end:     flip angle at end
%   tr_var:     echo number for [start to plateau; begin plateau; end plateau; reach rf_end]
%   example:
%   [rf,pow] = fliptraps_new(60,20,10,'opt',0,2,0,90,40,[7 10 15 20]);
function [rf,pow] = TRAPS_flip(rflip,nETL,k0,ftyp,traps,rf_high,rf_end,tr_var);

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
    end

case 'lin',
   dinc=(180-rflip)/(tr_var(1)+1)
    for k=1:tr_var(1),
        rf(k)=180-k*dinc;
    end

case 'sin'
  for k=1:tr_var(1),        
        rf(k)=180-(180-rflip)*sin(pi*k/(2*(tr_var(1)-1)));
  end; 

end;



if (traps==1),
    dinc=(rf_high-rf(tr_var(1)))/(tr_var(2)-tr_var(1));
    for k=tr_var(1)+1:tr_var(2),
        if (k>1),        rf(k)=rf(k-1)+dinc;    end;
    end;
    
    for k=tr_var(2)+1:tr_var(3),
        rf(k)=rf(k-1);
    end;
    dinc=(rf_end-rf(tr_var(3)))/(tr_var(4)-tr_var(3));
    for k=tr_var(3)+1:tr_var(4),
        rf(k)=rf(k-1)+dinc;
    end;
    for k=tr_var(4)+1:nETL,
        rf(k)=rf(k-1);
    end;
end;

if (traps==2),
    dinc=(rf_high-rf(tr_var(1)))/(tr_var(2)-tr_var(1));
    for k=tr_var(1)+1:tr_var(2),
        if (k>1),        
        rf(k)=rf_high-(rf_high-rf(tr_var(1)))*cos(pi*(k-tr_var(1))/(2*(tr_var(2)-tr_var(1))));
        end;
    end;
    rf(tr_var(2))=(rf_high+rf(tr_var(2)-1))/2;
    for k=tr_var(2)+1:tr_var(3),
        rf(k)=rf(k-1);
    end
    dinc=(rf_end-rf(tr_var(3)))/(tr_var(4)-tr_var(3));
    for k=tr_var(3)+1:tr_var(4),
        rf(k)=rf_high-(rf_high-rf_end)*cos(pi*(k-tr_var(4))/(2*(tr_var(4)-tr_var(3))));
        %rf(k)=(rf(tr_var(3))+rf_end)/2-(rf(tr_var(3))-rf_end)*0.5*cos(pi*(k-tr_var(4))/((tr_var(4)-tr_var(3)))); 
    end;
    
    for k=tr_var(4)+1:nETL,
        rf(k)=rf(k-1);
    end
end;

%plot(rf)
pow=sum((rf./180).^2)/nETL;





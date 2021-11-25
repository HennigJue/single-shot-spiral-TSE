% bloch_evo_r.m: calculates time-evolution of spins
%	usage: [w_array,t_array]= bloch_evo_r(mag,tn,winc,te,r1,r2)
%	input: 	
%	mag:		magnetization vector after excitation
%	tn:		number of time points;
%  winc:    max offresonance frequency (in Hz) 
%	te:		Duration of time-evolution interval (in ms);
%	r1, r2:	relaxation rates (in 1/s)
% 	written by J.Hennig;
%	tested on 02, 17.02.99;	
function m_out= pg_evo_cpmg(mag,te,r1,r2)

if(nargin<4),
   r1=0; r2=0;
end;
if(nargin<2),
   te=1;
end;
sm=size(mag);
sk=sm(2);
if (mag(1,sm(2))~=0),sk=sk+1; end;
m_out=zeros([sm(1) sk]);
m_out(1,2:sk)=mag(1,1:(sk-1));
m_out(2,1:(sm(2)-1))=mag(2,2:sm(2));
m_out(1,1)=mag(2,1);
m_out(3,1:sm(2))=mag(3,:);
m_out(4,1:sm(2))=mag(4,:);
m_out(1:2,:)=m_out(1:2,:)*exp(-te*r2/1000);
m_out(3:4,:)=m_out(3:4,:)*exp(-te*r1/1000);

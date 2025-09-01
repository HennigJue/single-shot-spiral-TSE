% calculate RF-power in a TSE-sequence during sequence construction
function [RFP,RFPtr,RFPrel] = calcRFP(rfex,rfref,acqP,ampref);
if (nargin<4),
    ampref=1;
end;
nEcho=length(acqP.flip);
dtEx=rfex.t(2)-rfex.t(1);
RFPex=sum(abs(rfex.signal).^2)*dtEx;
dtRef=rfref.t(2)-rfref.t(1);
RFPref0=sum(abs(rfref.signal).^2)*dtRef;
RFP=RFPex;
RFP0=RFP;
for k=1:nEcho
    RFP0=RFP0+RFPref0;
    RFP=RFP+RFPref0*acqP.flip(k).^2/180.^2;
end

% sprintf('%4.2f',RFP)
% sprintf('%1.4f',RFP/RFP0)
RFPtr=acqP.NSlices*RFP/acqP.TR;
RFPrel=RFP/RFP0;
end
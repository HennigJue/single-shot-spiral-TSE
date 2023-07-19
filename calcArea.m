function [area,tint,Gint] = calcArea(G)
system.gradRasterTime=1e-5;
tint=G.tt(1)+system.gradRasterTime/2:system.gradRasterTime:G.shape_dur-system.gradRasterTime/2;
Gint=interp1(G.tt,G.waveform,tint);
area=sum(Gint)*system.gradRasterTime;

np=length(G.tt);
area=0;
for k=2:np,
    area=area+(G.waveform(k)+G.waveform(k-1))/2*(G.tt(k)-G.tt(k-1));
end
end


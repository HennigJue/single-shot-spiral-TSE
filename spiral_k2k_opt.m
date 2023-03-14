function [Graster, traster, Gtran, tG, tdur] = spiral_k2k_opt(kBegin, kEnd, GBegin, GEnd, maxGrad, maxSlew, tsp, dW,pflag);
%function [Graster, traster, tdur] = spiral_k2k_opt(kBegin, kEnd, GBegin, GEnd, maxGrad, maxSlew, tsp, dW,pflag);

%%
if(nargin<9), pflag=0; end
k=1;
[GRasterx, tRasterx, GPoix, tPoix, toutx] = k2k_grads(kBegin(k),kEnd(k),GBegin(k),GEnd(k),maxGrad,maxSlew,tsp,dW, pflag);
k=2;
[GRastery, tRastery, GPoiy, tPoiy, touty] = k2k_grads(kBegin(k),kEnd(k),GBegin(k),GEnd(k),maxGrad,maxSlew,tsp,dW, pflag);
if((toutx||touty)>tsp),
    if(toutx>touty),
        k=2; [GRastery, tRastery, GPoiy, tPoiy, touty] = k2k_grads(kBegin(k),kEnd(k),GBegin(k),GEnd(k),maxGrad,maxSlew,toutx,dW, pflag);
    else
        k=1; [GRasterx, tRasterx, GPoix, tPoix, toutx] = k2k_grads(kBegin(k),kEnd(k),GBegin(k),GEnd(k),maxGrad,maxSlew,touty,dW, pflag);
    end
end
Graster=[GRasterx; GRastery]';
traster=tRasterx;

if (length(GPoix)>length(GPoiy)), GPoiy=[GPoiy GPoiy(end)];tPoiy=[tPoiy tPoiy(end)];  end
if (length(GPoiy)>length(GPoix)), GPoix=[GPoix GPoix(end)];tPoix=[tPoix tPoix(end)];  end
Gtran=[GPoix; GPoiy]';
tG=[tPoix; tPoiy]';
tdur=toutx;



        
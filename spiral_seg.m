% function [kseg, kseg2, nkseg, nkseg2, nkt] = spiral_traj(np,da,ninc,kt,mode)
% creates segmented spiral trajectory.
% parameters:
%
% Input:
%   kt:             k-space trajectory
%   ninc:           min. number of points in each segment
%   sp:             x-coordinate of starting point (y=0)
%   ntot:           total number of points per segment (must be > ninc)
%                   (if not available) trajectory will be created using spiral_calc
%   segmode            cont - continuous tight
%                   alt - alternating
% Output:
%   kseg:           segmented trajectories.
%   nkseq:          number of points in each segment
%   ikseg:          first and last index of full trajectory for each segment
%                   ikseg(k,2)-ikseg(k,1)=nkseg(k)+1; (except for k=1)
%   nkt:            number of segments.


function [kseg, nkseg, ikseg, nkt, i1] = spiral_seg(kt,ninc,sp,ntot,segmode)
if(nargin<5),
    segmode='cont'
end
%%
%close all
tinc=1e-5;  % precision of equispacing of datapoints
nc=100;      % number of points on tangent.
nsp=size(kt);
sp1=[-sp 0];
[tanxy1,itanxy10]=spiral_find_tangents(kt,sp1,nc,0);
ntan1=floor(length(itanxy10)/2);

sp2=[sp 0];
[tanxy2,itanxy20]=spiral_find_tangents(kt,sp2,nc,0);
ntan2=floor(length(itanxy20)/2);

itanxy1=itanxy10(ntan1+1:end);
itanxy2=itanxy20(ntan2+1:end);



k=1;
nkt=1;
kseg=zeros([ntot ceil(nkt)]);
ksegn=0*kseg;
nkta=1;
col=['b-';'r-'];
m=1;

ind=find(itanxy2(2:end)<(k-1)*ninc+ninc/2)+1;   % first segment half length
nkte=ind(end);
nktlast=itanxy2(nkte);
kt1=[kt(nkta:nktlast,1)';kt(nkta:nktlast,2)']';
nkseg(k)=length(kt1);
ikseg(k,1)=1;ikseg(k,2)=nktlast;
i1(k)=floor(ntot/2);
kseg(i1+1:i1+nkseg(k),k)=kt1(:,1)+i*kt1(:,2);

figure
pta(m)=itanxy1(nkta);ipta(m)=nkta;pte(m)=itanxy2(nkte);ipte(m)=nkte;lkt(m)=pte(m)-pta(m)+1;m=m+1;
sp1(1);
plot([sp1(1) ;kt1(:,1); sp2(1)],[0 ;kt1(:,2); 0],col(1,:))
col=flipud(col);
%col
hold on
plot([kt1(1,1) kt1(end,1)],[kt1(1,2) kt1(end,2)],'ro');
axis([-sp sp -sp sp])
axis square

%%

while((nsp(1)-nktlast)>ninc),
%for mtest=1:5,
    try,
        k=k+1;
        nkt=nkt+1;
        clear kt1 kt2
        
        if(strcmp(segmode,'cont')|(k~=2*floor(k/2))),
            
            itanxy1=itanxy10(ntan1+1:end);
            itanxy2=itanxy20(ntan2+1:end);
            
            ind1=find(itanxy1<nktlast);
            if(ind1(end)==nkta), disp('readout time too short'), nkt=nkt-2; return; end
            nkta=ind1(end);
            
            ind=find(itanxy2(2:end)<(itanxy1(nkta)+ninc));
            nkte=ind(end);
            ikseg(k,1)=itanxy1(nkta);ikseg(k,2)=itanxy2(nkte);
            kt1=[kt(itanxy1(nkta):itanxy2(nkte),1)';kt(itanxy1(nkta):itanxy2(nkte),2)']';
            nkseg(k)=length(kt1);
            i1(k)=floor((ntot-nkseg(k))/2);
            kseg(i1(k)+1:i1(k)+nkseg(k),k)=kt1(:,1)+i*kt1(:,2);
            nktlast=itanxy2(nkte);
            
        end
        
        if(strcmp(segmode,'alt')&&(k==2*floor(k/2))),
            
                        itanxy1=itanxy10(1:ntan1);
                        itanxy2=itanxy20(1:ntan2);
            %             ind1=find(itanxy1<nktlast);
            %             if(ind1(end)==nkta), disp('readout time too short'), nkt=nkt-2; return; end
            %             nkta=ind1(end);
                        ind1=find(itanxy2<(nktlast));
                        nkte=ind1(end);
            
                        ind=find(itanxy1<(itanxy2(nkte)+ninc));
                        nkta=ind(end);
            
                        ikseg(k,1)=itanxy1(nkta);ikseg(k,2)=itanxy2(nkte);
                        kt1=[kt(itanxy1(nkta):-1:itanxy2(nkte),1)';kt(itanxy1(nkta):-1:itanxy2(nkte),2)']';
                        nkseg(k)=length(kt1);
                        i1(k)=floor((ntot-nkseg(k))/2);
                        kseg(i1(k)+1:i1(k)+nkseg(k),k)=kt1(:,1)+i*kt1(:,2);
                        nktlast=itanxy1(nkta);
            
        end
        
        %figure
        pta(m)=itanxy1(nkta);ipta(m)=nkta;pte(m)=itanxy2(nkte);ipte(m)=nkte;lkt(m)=pte(m)-pta(m)+1;m=m+1;
        col(1,:);
        plot([sp1(1) ;kt1(:,1); sp2(1)],[0 ;kt1(:,2); 0],col(1,:))
        col=flipud(col);
        hold on
        plot([kt1(1,1) kt1(end,1)],[kt1(1,2) kt1(end,2)],'ro');
        axis([-sp sp -sp sp])
        axis square
    catch
        ninc=2*ninc; nkt=nkt-1; k=k-1; kseg=(kseg(1:ntot,:));
    end
    
end


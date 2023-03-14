% function bigx=csi2jpg(x,ny,nx,zoom,ninc);
% HAS BEEN MODIFIED FOR B/W-IMAGES, SCALING FACTOR SHOULD BE SET TO 0 FOR
% COLOUR IMAGES
% produces big image from multiple small images for easy transfer and save
% x:    original data matrix (4D with nx,ny,ns,nf)
% ny:   number of images in a column
% nx:   number of images in a row
% zoom: relative size of image to matrix
% ninc: image increment (default=1)

function bigx=im_mosaic(x,ny,nx,zoom,ifac,rsize,scale);

six=size(x);
if(nargin<7),
    scale=1;
end
if(nargin<6),
    rsize=[six(1) six(2)];
end


if(nargin<5),
    ifac=ones([six(3) 1]);
end
        
if(nargin<4),
    zoom=100;
end

if (ndims(x)==3),
    [sx,sy,nf]=size(x);
    x=reshape(x,[sx sy 1 nf]);
end

%mov = avifile(strcat(fn,'.avi'),'fps',4,'quality',100,'compression','Cinepak');
fig=figure;
ax1 = axes('position',[0,0,1,1]);
ax3 = axes('position',[0,0,1,1]);

colormap(gray)

% x=imread(files(1).name,'jpg');
[sx,sy,sn,nf]=size(x);
asp=sx/sy;
set(gcf,'Position',[100 100 round(nx*zoom) (asp*ny*zoom)]);
%bigx=uint8(zeros([ny*sx nx*sy sn]));
bigx=double(zeros([ny*sx nx*sy sn]));
size(bigx);
globsc=max(abs(x(:)));
%globsc=1;
kf=1;


    

for ky=1:ny,
        indy=(ky-1)*sx;
%         for kx=1:nx,
%     indx=(kx-1)*sy;
for kx=1:nx,
    indx=(kx-1)*sy;        
        %locsc=1;
      try, 
          xtemp=squeeze(x(:,:,:,kf));
          locsc=max(abs(xtemp(:)));
          scfac=1;
          %scfac=globsc/locsc;
          if(scfac>500), scfac=1; end;
     %     scfac=1;
     %     scfac=ifac(kf)
         % xtemp=xtemp./max(xtemp(:));
          
%           if(scale==1),
%           imax=max(xtemp(:))
%           imin=min(xtemp(:))
%           xtemp=(xtemp-imin)./(imax-imin);
%           xtemp=255*xtemp;
%           max(xtemp(:))
%           end
            %bigx(indy+1:indy+sx,indx+1:indx+sy,:)=globsc*squeeze(x(:,:,:,kf)); 
        bigx(indy+1:indy+sx,indx+1:indx+sy,:)=xtemp*scfac; 
        catch disp(strcat('image #',num2str(kf),'  is no valid image'));
            
        kf
        end
        kf=kf+1;
    end
end
try,
    %imagesc(real(bigx),[0 ifac*globsc]);
    imagesc(real(bigx));
catch
    %imagesc(abs(uint8(bigx)));
end
axis equal



% function Y= pg_pulse(mp,al,phi,profile)
%   calculates effect of a refocusing pulse via the transition matrix
%   formalism according to EPG in Notation compliant with Weiwei
%   requires:
%   mp      matrix of magnetization states [F1 F-1 Z1 Z-1,....]
%   al      flip angle
%   phi     phase of the pulse
%
%
function Y= pg_pulse_w(mp,al,phi,profile)
sm=size(mp); 
np=sm(2);
if (nargin<4),
   profile=1;
end;
if (nargin<3),
   phi=0;
end;

rmp=real(mp);
imp=imag(mp);

ca=cos(pi*phi/180);
sa=sin(pi*phi/180);
% 
% ca1=cos(pi*(phi-90)/180);
% sa1=sin(pi*(phi-90)/180);
% 
% c2a=cos(2*pi*phi/180);
% s2a=sin(2*pi*phi/180);
% 
% c2a1=cos(2*pi*(phi-90)/180);
% s2a1=sin(2*pi*(phi-90)/180);


% tc=[(cos(al*pi/360)).^2 -(sin(al*pi/360)).^2 sin(al*pi/180) 0;
% -(sin(al*pi/360)).^2 (cos(al*pi/360)).^2 0 -sin(al*pi/180);
% -0.5*sin(al*pi/180) -0.5*sin(al*pi/180) cos(al*pi/180) 0;
% 0.5*sin(al*pi/180) 0.5*sin(al*pi/180) 0 cos(al*pi/180)];

tc=[(cos(al*pi/360)).^2 -(sin(al*pi/360)).^2 sin(al*pi/180) 0;
-(sin(al*pi/360)).^2 (cos(al*pi/360)).^2 0 -sin(al*pi/180);
-0.5*sin(al*pi/180) -0.5*sin(al*pi/180) cos(al*pi/180) 0;
0.5*sin(al*pi/180) 0.5*sin(al*pi/180) 0 cos(al*pi/180)];

ts=[(cos(al*pi/360)).^2 (sin(al*pi/360)).^2 sin(al*pi/180) 0;
(sin(al*pi/360)).^2 (cos(al*pi/360)).^2 0 sin(al*pi/180);
-0.5*sin(al*pi/180) 0.5*sin(al*pi/180) cos(al*pi/180) 0;
0.5*sin(al*pi/180) -0.5*sin(al*pi/180) 0 cos(al*pi/180)];


% tc=[(cos(al*pi/360)).^2 -(sin(al*pi/360)).^2 sin(al*pi/180)/sqrt(2) 0;
% -(sin(al*pi/360)).^2 (cos(al*pi/360)).^2 0 -sin(al*pi/180)/sqrt(2);
% -sin(al*pi/180)/sqrt(2) -sin(al*pi/180)/sqrt(2) cos(al*pi/180) 0;
% sin(al*pi/180)/sqrt(2) sin(al*pi/180)/sqrt(2) 0 cos(al*pi/180)];
% 
% ts=[(cos(al*pi/360)).^2 (sin(al*pi/360)).^2 sin(al*pi/180/sqrt(2)) 0;
% (sin(al*pi/360)).^2 (cos(al*pi/360)).^2 0 sin(al*pi/180/sqrt(2));
% -sin(al*pi/180)/sqrt(2) sin(al*pi/180)/sqrt(2) cos(al*pi/180) 0;
% sin(al*pi/180)/sqrt(2) -sin(al*pi/180)/sqrt(2) 0 cos(al*pi/180)];

Y=0*mp;

np;
   for k=1:np,
        ry(:,k)=ca*rmp(:,k)-sa*imp(:,k);
        iy(:,k)=sa*rmp(:,k)+ca*imp(:,k);
        r1=ts*ry(:,k);
        i1=tc*iy(:,k);
  
   ry(:,k)=ca*r1+sa*i1;
   iy(:,k)=-sa*r1+ca*i1;
   end;
   %Y=iy+i*ry;
   Y=ry+i*iy;


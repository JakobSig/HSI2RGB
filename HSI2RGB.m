function [RGB,XYZ]=HSI2RGB(wY,Y,ydim,xdim,d,thresholdRGB)
% wY: wavelengths in nm
% Y : HSI as a (#pixels x #bands) matrix,
% dims: x & y dimension of image
% d: 50, 55, 65, 75, determines the illuminant used, if in doubt use d65
% thresholdRGB : set to small value greater than 0 if thesholding should be done to increase contrast 
%
%
% If you use this method, please cite the following paper:
%  M. Magnusson, J. Sigurdsson, S. E. Armansson, M. O. Ulfarsson, 
%  H. Deborah and J. R. Sveinsson, 
%  "Creating RGB Images from Hyperspectral Images Using a Color Matching Function", 
%  IEEE International Geoscience and Remote Sensing Symposium, Virtual Symposium, 2020
%
%  @INPROCEEDINGS{hsi2rgb, 
%  author={M. {Magnusson} and J. {Sigurdsson} and S. E. {Armansson} 
%  and M. O. {Ulfarsson} and H. {Deborah} and J. R. {Sveinsson}}, 
%  booktitle={IEEE International Geoscience and Remote Sensing Symposium}, 
%  title={Creating {RGB} Images from Hyperspectral Images using a Color Matching Function}, 
%  year={2020}, volume={}, number={}, pages={}}
%
% Paper is available at
% https://www.researchgate.net/profile/Jakob_Sigurdsson
%
%



load D_illuminants;
w=wxyz(:,1);
x=wxyz(:,2);
y=wxyz(:,3);
z=wxyz(:,4);

switch d    
    case 65
        i=1;
    case 50
        i=2;
    case 55 
        i=3;
    case 75
        i=4 ;   
end
I=D(:,i+1);
wI=D(:,1);




%interpolate to image wavelengths
I=interp1(wI,I,wY,'pchip','extrap')';
x=interp1(w,x,wY,'pchip','extrap')';
y=interp1(w,y,wY,'pchip','extrap')';
z=interp1(w,z,wY,'pchip','extrap')';

% Truncate at 780nm
i=find(wY>780);i=i(1);
Y=Y(:,1:i);
wY=wY(1:i);
I=I(1:i);
x=x(1:i);
y=y(1:i);
z=z(1:i);

% compute k
k=1/(trapz(wY,y.*I));
% Compute X,Y & Z for image
X=k*trapz(wY,Y*diag(I.*x),2);
Z=k*trapz(wY,Y*diag(I.*z),2);
Y=k*trapz(wY,Y*diag(I.*y),2);

XYZ=[X Y Z]';



%if useM
%    sRGB=xyz2rgb(XYZ','ColorSpace','sRGB','WhitePoint','d65')';
%else

% Convert to RGB
M=[ 3.2404542 -1.5371385 -0.4985314;
   -0.9692660  1.8760108  0.0415560;
    0.0556434 -0.2040259  1.0572252];
sRGB=M*XYZ;




% Correct gamma 
gamma_map = (sRGB >  0.0031308);
sRGB(gamma_map) = 1.055 * power(sRGB(gamma_map), (1. / 2.4)) - 0.055;
sRGB(~gamma_map) = 12.92*sRGB(~gamma_map);
sRGB(sRGB>1)=1;
sRGB(sRGB<0)=0;

if (thresholdRGB>0)
thres=thresholdRGB;
for idx=1:3
    y=sRGB(idx,:);
    [a,b]=hist(y(y>0),100);
    a=cumsum(a)/sum(a);
    th=b(1);
    i=find(a<thres);
    if ~isempty(i)
    th=b(i(end));
    end
    y=max(0,y-th);

    
    [a,b]=hist(y,100);
    a=cumsum(a)/sum(a);    
    i=find(a>1-thres);
    th=b(i(1));
    y(y>th)=th;
    y=y./th;
    sRGB(idx,:)=y;
end
end

RGB(:,:,1)=reshape(sRGB(1,:),ydim,xdim);
RGB(:,:,2)=reshape(sRGB(2,:),ydim,xdim);
RGB(:,:,3)=reshape(sRGB(3,:),ydim,xdim);


end

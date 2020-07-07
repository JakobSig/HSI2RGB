load washington_hsi.mat
% Y holds the HSI (uint8 to save space
% wl holds the wavelengths of the spectral bands
% The image is 307x307 pixels
% The illuminant used is D65
% 


% load data and convert to double
Y=double(Y)/255;
[ydim,xdim,zdim]=size(Y);

% reorder data so that each column holds the spectra of of one pixel
for i=1:zdim
    z=Y(:,:,i);z=z(:);
    Z(:,i)=z;
end

% use the D65 illuminant
illuminant=65;

% do minor thresholding
dothreshold=1;

%Create the RBG image, 
[RGB,XYZ]=HSI2RGB(wl,Z,ydim,xdim,illuminant,dothreshold);
imagesc(RGB)

function [ imagerot ] = imagerot(image,theta)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
imagepad = image;
[nrows ncols] = size(imagepad);
midx=ceil((ncols+1)/2);
midy=ceil((nrows+1)/2);

Mr = [cos(theta) sin(theta); -sin(theta) cos(theta)]; % rotation matrix

% rotate about center
[X Y] = meshgrid(1:ncols,1:nrows);
XYt = [X(:)-midx Y(:)-midy]*Mr;
XYt = bsxfun(@plus,XYt,[midx midy]);

xout = round(XYt(:,1)); yout = round(XYt(:,2)); % nearest neighbor!
outbound = yout<1 | yout>nrows | xout<1 | xout>ncols;
%zout=repmat(cat(3,1,2,3),nrows,ncols,1); zout=zout(:);
xout(xout<1) = 1; xout(xout>ncols) = ncols;
yout(yout<1) = 1; yout(yout>nrows) = nrows;
xout_2 = repmat(xout,[2 1]); yout_2 = repmat(yout,[2 1]);
imagerot = imagepad(sub2ind(size(imagepad),yout,xout));%zout(:))); % lookup
imagerot = reshape(imagerot,size(imagepad));
imagerot(repmat(outbound,[1 1])) = 0; % set background value to [0 0 0] (black)
end


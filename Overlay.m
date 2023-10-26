function Overlay(Slice1, Slice2, Alpha, Limits, op1)

% figure
baseImage=Slice1;
parameterROIImage=Slice2;
baseImage = baseImage/(max(baseImage(:))); % normalize base (anatomical) image
rgbSlice  = baseImage(:,:,[1 1 1]);        % converting to RGB (ignore colormaps)
if ~strcmp(version('-release'),'2018b')
    imagesc(parameterROIImage)
else
    imshow(parameterROIImage, []);             % show parametric image
end
colormap('jet');                           % apply colormap
hold on;
h = imshow(rgbSlice);                      % superimpose anatomical image
set(h, 'AlphaData', Alpha);      % make pixels in the ROI transparent
if nargin==5
    if op1==1
        colorbar;
    end
else
    colorbar;  
end
% caxis([0 30])
caxis(Limits)
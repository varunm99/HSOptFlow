function [Ix, Iy] = ImageGradient(image)
    % Boundary condition: Extend image by 1 pixel in each direction
    imageExt = zeros(size(image,1)+2, size(image, 2) + 2);
    imageExt(2:end-1,2:end-1) = image; 
    imageExt(2:end-1,1) = image(:,1);
    imageExt(2:end-1,end) = image(:,end);
    imageExt(1,2:end-1) = image(1,:);
    imageExt(end,2:end-1) = image(end,:);
    
    
    % Forward Difference Gradient
%     Ix = conv2(imageExt, [1 -1], 'same');
%     Iy = conv2(imageExt, [1;-1], 'same');
    
    % Central Difference Gradient
    Ix = conv2(imageExt, 0.5*[1 0 -1], 'same');
    Iy = conv2(imageExt, 0.5*[1;0;-1], 'same');

    % Remove excess lines from boundary convolution result
    Ix = Ix(2:end-1, 2:end-1);
    Iy = Iy(2:end-1, 2:end-1);
end
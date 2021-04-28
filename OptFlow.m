function [u, v] = OptFlow(currFrame, prevFrame, alpha)
u = zeros(size(currFrame));
v = zeros(size(currFrame));

% Gaussian filter preprocessing to smooth image
currFrame = imgaussfilt(currFrame, 1);
prevFrame = imgaussfilt(prevFrame, 1);

% Calculate image derivatives
[Ix, Iy] = ImageGradient(prevFrame);
It = prevFrame - currFrame;

avg_kernel =0.25*[0 1 0; 1 0 1; 0 1 0]; % Averaging kernel for Laplacian

iterations = 50;
for i=1:iterations
    % Compute average vector components for Laplacian
    u0=conv2(u,avg_kernel,'same');
    v0=conv2(v,avg_kernel,'same');
    
    % Update linear system solution across entire image
    u = u0 - (Ix.*(Ix.*u0 + Iy.*v0 + It))./(alpha + Ix.^2 + Iy.^2); 
    v = v0 - (Iy.*(Ix.*u0 + Iy.*v0 + It))./(alpha + Ix.^2 + Iy.^2);
end
u(isnan(u))=0;
v(isnan(v))=0;
end
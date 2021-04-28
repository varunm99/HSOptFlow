close all
%% Read in Frames

prevFrame = ones(480, 640);
l = 30;
w = 30;
x0 = 480/2;
y0 = 640/2;
posX = x0;
posY = y0;
prevFrame(x0:x0+l, y0:y0+w) = 0;


T = 40;
xtrue = zeros(4,T);
x = zeros(4,T);
xhat = zeros(4,T+1);
xhat(:,1) = [x0; y0; 0; 0];
Q = diag([1, 1, 1000, 1000]);
R = diag([1, 1, 1000,1000]);

P = Q;




%% Opt Flow Algorithm

for i = 1:T
    currFrame = ones(480, 640);
    
    velX = round(5*sin(pi/8*i));
    velY = round(5*cos(pi/8*i));
%     velX = 4;
%     velY = 0;
    posX = velX + posX;
    posY = velY + posY;
    
    xtrue(:,i) = [posX; posY; velX; velY];
    currFrame(posX:posX+l, posY:posY+w) = 0;
    
    
    alpha = 1;
    [u,v]= OptFlow(prevFrame, currFrame, alpha);
   
    %% Block Pos Estimation
    matrix = (1-currFrame);
    matrix=matrix/sum(matrix(:));
    [m,n]=size(matrix);
    [I,J]=ndgrid(1:m,1:n);
    c=[dot(I(:),matrix(:)),  dot(J(:),matrix(:))] - [l/2, w/2]
    x(1:2,i)= c';
    velxy = [sum(v(:).*(1-prevFrame(:)))/(w*l); sum(u(:).*(1-prevFrame(:)))/(w*l)]*10
    x(3:4, i) = velxy;
    
    %% Kalman Filtering
    % state vector = [x, y, vx, vy]^T
    
    % Measurement is [x, y, vx, vy]
    measure = @(x) x; % identity
    H = eye(4);
    
    % Predict
    F = [1 0 1 0; 0 1 0 1; 0 0 1 0; 0 0 0 1];
    xhat(:,i+1) = F*xhat(:,i);
    P = F*P*F' + Q;
    error = measure(x(:,i)) - H*xhat(:,i+1);
    S = H*P*H' + R;
    K = P*H'*inv(S);
    xhat(:,i+1) = xhat(:,i+1) + K*error;
    P = (eye(4) - K*H)*P;
    
    
    

    
    %% Plotting
    imshow(currFrame);
    hold on;
    step = 5;
    A = 1:step:size(prevFrame,1);
    B = 1:step:size(prevFrame,2);
    [X, Y] = meshgrid(B,A);
    h = quiver(X, Y, u(A,B), v(A,B), 'color', 'b');
    visscale = 2;
    b = quiver(xhat(2,i) + l/2,  xhat(1,i) + w/2, visscale*xhat(4,i), visscale*xhat(3,i), 'color', 'r');
    drawnow
    prevFrame = currFrame;
end

figure
subplot(2,2,1)
plot(x(1,:));
hold on
plot(xhat(1,2:end));
plot(xtrue(1,:), '--');
legend("X", "Xhat", 'True');
title("X Position");
subplot(2,2,2)
plot([0 diff(x(1,:))]);
hold on
plot(x(3,:));
plot(xhat(3,2:end));
plot(xtrue(3,:), '--');
legend("Pos Diff", "X", "Xhat", 'True');
title("X Velocity");
ylim([-6, 6]);


subplot(2,2,3)
plot(x(2,:));
hold on
plot(xhat(2,2:end));
plot(xtrue(2,:), '--');
legend("X", "Xhat", "True");
title("Y Position");
ylim([200, 400]);
subplot(2,2,4)
plot([0 diff(x(2,:))]);
hold on
plot(x(4,:));
plot(xhat(4,2:end));
plot(xtrue(4,:), '--');
legend("Pos Diff", "X", "Xhat", "True");
title("Y Velocity");
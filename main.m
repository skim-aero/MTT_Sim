%% Multi target tracking with GNN and JPDA
% Kernel GNN and kernel JPDA are included (DA in kernel space)
% Sukkeun Kim (Sukkeun.Kim@cranfield.ac.uk)

clc; clear; close all;
addpath("functions")

rng(42);

iter = 5;   % iteration number

%% Define a Scenario
% Time
T = 1.0;    % sampling time
t = 0:T:100;    % time

% Parameters
n = 6; % no of state

xr = 5000; yr = 3000; % radar position

d2r = pi/180; r2d = 1/d2r; 
wclttuer = 1;   % 1 if with clutter
no = 10;    % number of measurement + number of clutter

% Target 1
sig_r1 = 10;    % standard deviation of range measurement noise 10
sig_az1 = 3*d2r;    % standard deviation of azimuth measurement noise 3.5

% Target 2
sig_r2 = 10;    % standard deviation of range measurement noise
sig_az2 = 3*d2r;    % standard deviation of azimuth measurement noise

% Clutter
sig_rc = 50;
sig_azc = 5*d2r;

% Process noise covariance matrix
q = 0.1; Q = q^2;

% Measurement noise covariance matrix
R1 = [sig_r1^2,0; 0,sig_az1^2]; % for target 1
R2 = [sig_r2^2,0; 0,sig_az2^2]; % for target 2

% System/measurement equations
F1 = [1,T,T^2/2; 0,1,T; 0,0,1];
F = [F1, zeros(size(F1)); zeros(size(F1)),F1];

G1 = [T^3/6; T^2/2; T];
G = [G1; G1];

%% For evaluation
error = zeros(4,1);

RMSE = zeros(4,iter+1);
SRMSE = zeros(4,1);
ARMSE = zeros(4,1);

%% True target trajectory generation
x1 = zeros(1,length(t));
y1 = zeros(1,length(t));
x2 = zeros(1,length(t));
y2 = zeros(1,length(t));

x1(1) = 0; y1(1) = 0; x2(1) = 0; y2(1) = 3000; % departure position
vx1 = 2; ax1 = 0; vy1 = 0.2; ay1 = 0.01; % x&y velocity and acceleration
vx2 = 2; ax2 = 0; vy2 = -2; ay2 = 0.05; % x&y velocity and acceleration

% Trajectory generation by a constant acc motion
for i = 1:length(t)-1
    x1(i+1) = x1(i)+vx1*t(i)+0.5*ax1*t(i)^2;
    y1(i+1) = y1(i)+vy1*t(i)+0.5*ay1*t(i)^2;

    x2(i+1) = x2(i)+vx2*t(i)+0.5*ax2*t(i)^2;
    y2(i+1) = y2(i)+vy2*t(i)+0.5*ay2*t(i)^2;
end

for it = 1:iter
    disptemp = ['Iteration: ', num2str(it)];
    disp(disptemp)

    % Noise-corrupted radar measurements in spherical coordinate
    r1 = sqrt((x1-xr).^2+(y1-yr).^2);   % range target 1
    az1 = atan2(y1-yr,x1-xr);   % azimuth target 1

    r2 = sqrt((x2-xr).^2+(y2-yr).^2);   % range target 2
    az2 = atan2(y2-yr,x2-xr);   % azimuth target 2
    
    temp = zeros(no,length(t));
    
    % Noise
    for i = 1:no
        temp(i,:) = randn(length(t),1);
    end
    
    % Target 1
    rm1 = r1+sig_r1*temp(1,:);    % noise-corrupted range measurement
    azm1 = az1+sig_az1*temp(2,:); % noise-corrupted azimuth measurement
    
    % Target 2
    rm2 = r2+sig_r2*temp(3,:);  % noise-corrupted range measurement
    azm2 = az2+sig_az2*temp(4,:);   % noise-corrupted azimuth measurement
    
    % Clutter    
    rmc1 = [r1+sig_rc*temp(5,:); r1+sig_rc*temp(6,:); r1+sig_rc*temp(7,:)];
    azmc1 = [az1+sig_azc*temp(5,:); az1+sig_azc*temp(6,:); az1+sig_azc*temp(7,:)];
    
    rmc2 = [r2+sig_rc*temp(8,:); r2+sig_rc*temp(9,:); r2+sig_rc*temp(10,:)];
    azmc2 = [az2+sig_azc*temp(8,:); az2+sig_azc*temp(9,:); az2+sig_azc*temp(10,:)];
    
    rmc = [rm1; rm2; rmc1; rmc2];
    azmc = [azm1; azm2; azmc1; azmc2];
    
    [xtemp,ytemp] = pol2cart(azmc,rmc);
    
    xnew = xtemp+xr;
    ynew = ytemp+yr;
    
    %% Initialisation of saving variable
    % State and error covariance for target 1
    X1 = zeros(n,length(t),2);
    P_f1 = zeros(n,n,length(t),2);
    
    % State and error covariance for target 2
    X2 = zeros(n,length(t),2);
    P_f2 = zeros(n,n,length(t),2);
    
    %% Initial condition
    for j = 1:4
        X1(:,1,j) = [x1(1)-50,0,0,y1(1)-50,0,0]';
        P_f1(:,:,1,j) = 1e3*eye(n,n);

        X2(:,1,j) = [x2(1)-50,0,0,y2(1)+50,0,0]';
        P_f2(:,:,1,j) = 1e3*eye(n,n);
    end
    
    %% Main loop
    for i = 1:length(t)-1
        % Make the absolute value of angles under 180deg
        if abs(azmc(i))>pi
            azmc(i) = azmc(i)-2*pi*sign(azmc(i));
        end
    
        %% Prediction step
        for j = 1:4
            X1(:,i+1,j) = F*X1(:,i,j);
            P_f1(:,:,i+1,j) = F*P_f1(:,:,i,j)*F'+G*Q*G';

            X2(:,i+1,j) = F*X2(:,i,j);
            P_f2(:,:,i+1,j) = F*P_f2(:,:,i,j)*F'+G*Q*G';
        end

        %% Data Association and update step
        % 1) GNN
        [rm1,azm1,rm2,azm2] = GNN(xr,yr,rmc(:,i),azmc(:,i),X1(:,i+1,1),X2(:,i+1,1),no,wclttuer);
        [X1(:,i+1,1),P_f1(:,:,i+1,1)] = EKF(xr,yr,rm1,azm1,X1(:,i+1,1),P_f1(:,:,i+1,1),R1);
        [X2(:,i+1,1),P_f2(:,:,i+1,1)] = EKF(xr,yr,rm2,azm2,X2(:,i+1,1),P_f2(:,:,i+1,1),R2);

        % 2) KGNN
        [rm11,azm11,rm21,azm21] = KGNN(xr,yr,rmc(:,i),azmc(:,i),X1(:,i+1,2),X2(:,i+1,2),no,wclttuer);
        [X1(:,i+1,2),P_f1(:,:,i+1,2)] = EKF(xr,yr,rm11,azm11,X1(:,i+1,2),P_f1(:,:,i+1,2),R1);
        [X2(:,i+1,2),P_f2(:,:,i+1,2)] = EKF(xr,yr,rm21,azm21,X2(:,i+1,2),P_f2(:,:,i+1,2),R2);
    
        % 3) JPDA with assumption that at least one measurement associated
        [X1(:,i+1,3),P_f1(:,:,i+1,3),X2(:,i+1,3),P_f2(:,:,i+1,3)] = ...
            JPDA(xr,yr,rmc(:,i),azmc(:,i),X1(:,i+1,3),P_f1(:,:,i+1,3),...
            X2(:,i+1,3),P_f2(:,:,i+1,3),R1,R2);

        % 4) KJPDA with assumption that at least one measurement associated
        [X1(:,i+1,4),P_f1(:,:,i+1,4),X2(:,i+1,4),P_f2(:,:,i+1,4)] = ...
            KJPDA(xr,yr,rmc(:,i),azmc(:,i),X1(:,i+1,4),P_f1(:,:,i+1,4),...
            X2(:,i+1,4),P_f2(:,:,i+1,4),R1,R2,no);
    end

    % For evaluation 
    for i = 1:length(t)
        for j = 1:4
            error(j) = error(j)+(X1(1,i,j)-x1(i))^2+(X1(4,i,j)-y1(i))^2;
        end
    end
    
    for j = 1:4
        error(j) = sqrt(error(j)/(length(t)-1));
        RMSE(j,it) = error(j);
        SRMSE(j) = SRMSE(j)+error(j);
    end
end

%% Evaluation
for j = 1:4
    ARMSE(j) = SRMSE(j)/iter;
    RMSE(j,it+1) = ARMSE(j);
end

txt_1 = ['GNN: ', num2str(ARMSE(1))];
txt_2 = ['Kernel GNN: ', num2str(ARMSE(2))];
txt_3 = ['JPDA: ', num2str(ARMSE(3))];
txt_4 = ['Kernel JPDA: ', num2str(ARMSE(4))];

%% Plot: trajectoreis
figure;

% True trajectories
plot(x1(1:end-1),y1(1:end-1),':','LineWidth',4); 
xlabel('East, m'), ylabel('North, m'); 
hold on;
plot(x2(1:end-1),y2(1:end-1),':','LineWidth',4);

% Start and end points
plot(x1(1),y1(1),'o'); text(x1(1),y1(1)+300,'Start target 1','fontsize',12)
plot(x1(end-1),y1(end-1),'x'); text(x1(end-1),y1(end-1)+300,'End target 1','fontsize',12)
plot(x2(1),y2(1),'o'); text(x2(1),y2(1)+300,'Start target 2','fontsize',12)
plot(x2(end-1),y2(end-1),'x'); text(x2(end-1),y2(end-1)+300,'End target 2','fontsize',12)

% Estiamted tracks GNN
plot(X1(1,:,1),X1(4,:,1),'b','LineWidth',2); 
plot(X2(1,:,1),X2(4,:,1),'r','LineWidth',2); 

% Estiamted tracks KGNN
plot(X1(1,:,2),X1(4,:,2),'b-o','LineWidth',2); 
plot(X2(1,:,2),X2(4,:,2),'r-o','LineWidth',2); 

% Estiamted tracks JPDA
plot(X1(1,:,3),X1(4,:,3),'b-*','LineWidth',2); 
plot(X2(1,:,3),X2(4,:,3),'r-*','LineWidth',2); 

% Estiamted tracks KJPDA
plot(X1(1,:,4),X1(4,:,4),'b-+','LineWidth',2); 
plot(X2(1,:,4),X2(4,:,4),'r-+','LineWidth',2); 

% Measurements including clutters
scatter(xnew(1,1:end-1),ynew(1,1:end-1),'square','filled',color='#4DBEEE');
scatter(xnew(2,1:end-1),ynew(2,1:end-1),'m','square','filled');

if wclttuer == 1
    scatter(xnew(3:8,1:end-1),ynew(3:8,1:end-1),10,'k',"^",'filled');
end

% Radar site
plot(xr,yr,'x');
text(xr+40,yr,'Radar 1','fontsize',12);

title('Result');
legend('','','','','','',...
       'T1-GNN','T2-GNN',...
       'T1-KGNN','T2-KGNN',...
       'T1-JPDA','T2-JPDA',...
       'T1-KJPDA','T2-KJPDA');
axis([-1000 11000 -1000 4000]);

%% Plot: RMSE
figure;
plot(RMSE(1,1:iter),'b','LineWidth',2); 
xlabel('Iteration'), ylabel('RMSE, m'); 
hold on;

plot(RMSE(2,1:iter),'g-o','LineWidth',2); 
plot(RMSE(3,1:iter),'r-*','LineWidth',2); 
plot(RMSE(4,1:iter),'k-+','LineWidth',2); 

txt = {'Avg RMSE',txt_1,txt_2,txt_3,txt_4};

title('RMSE');
legend('GNN','Kernel GNN','JPDA','Kernel JPDA');
text(2.5,250,txt,'FontSize',16)

function [X, P_f] = EKF(xr,yr,xm,ym,X_s,P_s,R)
% Simple EKF

% Nonlinear measurement eq
h = hk(xr,yr,X_s);

% Jacobian of nonlinear measurement eq
H = JH(xr,yr,X_s);
v = ([xm; ym]-h);

% Make the absolute value of angles under 180deg
if abs(v(2))>pi
    v(2) = v(2)-2*pi*sign(v(2));
end

% Update step
S = H*P_s*H'+R;
K = P_s*H'/S;   % Kalman gain
X = X_s+K*(v);
P_f = (eye(size(K*H))-K*H)*P_s;

% h(X): Nonlinear measurement eq
function h = hk(xr,yr,X_s)
x = X_s(1);
y = X_s(4);

h = [sqrt((x-xr)^2+(y-yr)^2); atan2(y-yr,x-xr)];

% Partial_h/partial_X: Jacobian of measurement eq
function H = JH(xr,yr,X_s)
x = X_s(1);
y = X_s(4);

temp1 = sqrt((x-xr)^2+(y-yr)^2);
temp2 = cos(atan2(y-yr,x-xr))^2;

H = [(x-xr)/temp1,0,0,(y-yr)/temp1,0,0;
     temp2*(-(y-yr)/(x-xr)^2),0,0,temp2/(x-xr),0,0];

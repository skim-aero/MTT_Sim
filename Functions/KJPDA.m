function [X1, P_f1, X2, P_f2] = KJPDA(xr,yr,xm,ym,...
                                           X_s1,P_s1,X_s2,P_s2,R1,R2,no)
% Kernel JPDA with EKF (currently only for two targets)
% The data association is done in kernel space

% Parameter setting
Pfa = 0.02; % probability of the false alarm
Pd = 0.98;  % probability of the detection
gate = 0.8; % gating 5.9915;

% Nonlinear measurement eq
h1 = hk(xr,yr,X_s1);
h2 = hk(xr,yr,X_s2);

% Jacobian of nonlinear measurement eq.
H1 = JH(xr,yr,X_s1);
H2 = JH(xr,yr,X_s2);

% Update step
S1 = H1*P_s1*H1'+R1;
S2 = H2*P_s2*H2'+R2;
K1 = P_s1*H1'/S1;   % Kalman gain
K2 = P_s2*H2'/S2;

% Innovation
vj1 = zeros(length(h1),length(xm(:)));
vj2 = zeros(length(h2),length(xm(:)));

for j = 1:length(xm(:))
    vj1(:,j) = ([xm(j); ym(j)]-h1);
    vj2(:,j) = ([xm(j); ym(j)]-h2);
    
    % Make the absolute value of angles under 180deg
    if abs(vj1(2,j))>pi
        vj1(2,j) = vj1(2,j)-2*pi*sign(vj1(2,j));
    end
    if abs(vj2(2,j))>pi
        vj2(2,j) = vj2(2,j)-2*pi*sign(vj2(2,j));
    end
end

% For kernel JPDA
temp1 = zeros(3,2);
temp2 = zeros(no-2,2);
temp3 = zeros(3,no-2);

temp1(:,1) = funker(X_s1(1),X_s1(4));
temp1(:,2) = funker(X_s2(1),X_s2(4));

temp2(:,1) = (xm.*cos(ym))+xr;
temp2(:,2) = (xm.*sin(ym))+yr;

for i = 1:no-2
    temp3(:,i) = funker(temp2(i,1),temp2(i,2));
end

% Innovation in kernel space
vjk1 = zeros(length(temp1),length(xm(:)));
vjk2 = zeros(length(temp1),length(xm(:)));

for j = 1:length(xm(:))
    vjk1(:,j) = (temp3(:,j)-temp1(:,1));
    vjk2(:,j) = (temp3(:,j)-temp1(:,2));
end

% Measurement covariance in kernel space
% Jacobian of nonlinear transfrom eq
T1 = JT(xr,yr,h1);
T2 = JT(xr,yr,h2);

K = 0.1;
Sk1 = (T1*S1*T1'+K*eye(3))*1E14;
Sk2 = (T2*S2*T2'+K*eye(3))*1E14;

% Probability calculation
D1 = zeros(1,length(xm(:)));
D2 = zeros(1,length(xm(:)));
g1 = zeros(1,length(xm(:)));
g2 = zeros(1,length(xm(:)));

for j = 1:length(xm(:))
    D1(j) = vjk1(:,j)'/Sk1*vjk1(:,j);
    D2(j) = vjk2(:,j)'/Sk2*vjk2(:,j);

    g1(j) = exp(-D1(j)/2)/((2*pi)^(length(temp1)/2)*sqrt(norm(Sk1)));
    g2(j) = exp(-D2(j)/2)/((2*pi)^(length(temp1)/2)*sqrt(norm(Sk2)));

    if D1(j) > gate
        D1(j) = inf;
    end

    if D2(j) > gate
        D2(j) = inf;
    end
end

[D1sort,idx1] = sort(D1);
[D2sort,idx2] = sort(D2);

g1 = g1(idx1);
g2 = g2(idx2);

vj1(:,:) = vj1(:,idx1);
vj2(:,:) = vj2(:,idx2);

pj1 = zeros(1,length(xm(:)));
pj2 = zeros(1,length(xm(:)));

% For the time being no consideration of no assignment.
for j = 1:length(xm(:))
    if D1sort(j) == inf
        break
    end
    
    for K = 1:length(xm(:))
        if D2sort(K) == inf
        end

        pj1(j) = pj1(j)+g1(j)*Pd*g2(K)*Pd*Pfa^6;
    end
end

for j = 1:length(xm(:))
    if D2sort(j) == inf
        break
    end
    
    for K = 1:length(xm(:))
        if D1sort(K) == inf
        end

        pj2(j) = pj2(j)+g2(j)*Pd*g1(K)*Pd*Pfa^6;
    end
end

pj1 = normalize(pj1,"norm",1);
pj2 = normalize(pj2,"norm",1);

v1 = zeros(length(h1),1);
v2 = zeros(length(h1),1);

for j = 1:length(xm(:))
    v1 = v1+pj1(j)*vj1(:,j);
    v2 = v2+pj2(j)*vj2(:,j);
end

p_01 = 0.0;
p_02 = 0.0;

% State update
X1 = X_s1+K1*(v1);
X2 = X_s2+K2*(v2);

% Error covariance correction term
temp1 = 0;
temp2 = 0;

for j = 1:length(xm(:))
    temp1 = temp1+pj1(j)*vj1(:,j)*vj1(:,j)';
    temp2 = temp2+pj2(j)*vj2(:,j)*vj2(:,j)';
end

P_k1 = K1*(temp1-v1*v1')*K1';
P_k2 = K2*(temp2-v2*v2')*K2';

% Error covariance update
P_f1 = P_s1-(1-p_01)*K1*H1*P_s1+P_k1;
P_f2 = P_s2-(1-p_02)*K2*H2*P_s2+P_k2;

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

% partial_h/partial_X: Jacobian of transform matrix (to kernel space)
function T = JT(xr,yr,h)
r=h(1);
t=h(2);

T=[2*r*(cos(t)^2),-r^2*sin(2*t);
   2*r*(sin(t)^2),r^2*sin(2*t);
   sqrt(2)*2*r*cos(t)*sin(t),sqrt(2)*r^2*cos(2*t)];

% Kerenl funciton
function kerval = funker(x,y)
kerval = [x^2;y^2;sqrt(2)*x*y];


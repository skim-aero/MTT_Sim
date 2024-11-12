function [X1, P_f1, X2, P_f2] = JPDA(xr,yr,xm,ym,...
                                           X_s1,P_s1,X_s2,P_s2,R1,R2)
    % JPDA with EKF (currently only for two targets)
    
    % Parameter setting
    Pfa = 0.02; % probability of the false alarm
    Pd = 0.98;  % probability of the detection
    gate = 5.9915;  % gating
    
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
    
    % Probability calculation
    D1 = zeros(1,length(xm(:)));
    D2 = zeros(1,length(xm(:)));
    g1 = zeros(1,length(xm(:)));
    g2 = zeros(1,length(xm(:)));
    
    for j = 1:length(xm(:))
        D1(j) = vj1(:,j)'/S1*vj1(:,j);
        D2(j) = vj2(:,j)'/S2*vj2(:,j);
    
        g1(j) = exp(-D1(j)/2)/((2*pi)^(length(h1)/2)*sqrt(norm(S1)));
        g2(j) = exp(-D2(j)/2)/((2*pi)^(length(h2)/2)*sqrt(norm(S2)));
    
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
        end
        
        for k = 1:length(xm(:))
            if D2sort(k) == inf
            end
    
            pj1(j) = pj1(j)+g1(j)*Pd*g2(k)*Pd*Pfa^6;
        end
    end
    
    for j = 1:length(xm(:))
        if D2sort(j) == inf
        end
        
        for k = 1:length(xm(:))
            if D1sort(k) == inf
            end
    
            pj2(j) = pj2(j)+g2(j)*Pd*g1(k)*Pd*Pfa^6;
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
end

% h(X): Nonlinear measurement eq
function h = hk(xr,yr,X_s)
    x = X_s(1);
    y = X_s(4);
    
    h = [sqrt((x-xr)^2+(y-yr)^2); atan2(y-yr,x-xr)];
end

% Partial_h/partial_X: Jacobian of measurement eq
function H = JH(xr,yr,X_s)
    x = X_s(1);
    y = X_s(4);
    
    temp1 = sqrt((x-xr)^2+(y-yr)^2);
    temp2 = cos(atan2(y-yr,x-xr))^2;
    
    H = [(x-xr)/temp1,0,0,(y-yr)/temp1,0,0;
         temp2*(-(y-yr)/(x-xr)^2),0,0,temp2/(x-xr),0,0];
end
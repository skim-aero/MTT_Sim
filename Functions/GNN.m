function [rm1, azm1, rm2, azm2] = GNN(xr,yr,rmc,azmc,X1,X2,no,wclttuer)
% GNN data association (currently only for two targets)

temp1 = zeros(2,2);
d = zeros(2,no-2);

temp1(:,1) = funmeas(xr,yr,X1);
temp1(:,2) = funmeas(xr,yr,X2);

for j = 1:2
    for k = 1:no-2
        d(j,k) = fundist(temp1(1,j),temp1(2,j),rmc(k),azmc(k));
    end
end

if wclttuer == 1
    % W/ clutter
    dtemp1 = d(1,:);
    dtemp2 = d(2,:);
else       
    % W/O clutter
    dtemp1 = d(1,1:2);
    dtemp2 = d(2,1:2);
end

[dsort1,idx1] = sort(dtemp1);
[dsort2,idx2] = sort(dtemp2);

% GNN (only for two targets)
if idx1(1) ~= idx2(1)
    rm1 = rmc(idx1(1));
    azm1 = azmc(idx1(1));

    rm2 = rmc(idx2(1));
    azm2 = azmc(idx2(1));
else
    dsum1 = dsort1(1) + dsort2(2);
    dsum2 = dsort1(2) + dsort2(1);

    if dsum1 < dsum2
        rm1 = rmc(idx1(1));
        azm1 = azmc(idx1(1));
    
        rm2 = rmc(idx2(2));
        azm2 = azmc(idx2(2));
    else
        rm1 = rmc(idx1(2));
        azm1 = azmc(idx1(2));
    
        rm2 = rmc(idx2(1));
        azm2 = azmc(idx2(1));
    end
end

% Measurement function
function meas = funmeas(xr,yr,X1)
x = X1(1);
y = X1(4);
meas = [sqrt((x-xr)^2+(y-yr)^2); atan2(y-yr,x-xr)];

% Distance function
function dist = fundist(a,b,c,d)
dist = sqrt(a^2+c^2-2*a*c*cos(d-b));



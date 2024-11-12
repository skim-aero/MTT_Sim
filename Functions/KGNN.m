function [rm1, azm1, rm2, azm2] = KGNN(xr,yr,rmc,azmc,X1,X2,no,wclttuer)
    % Kernel GNN data association (currently only for two targets)
    % The data association is done in kernel space
    
    temp1 = zeros(3,2);
    temp2 = zeros(no-2,2);
    temp3 = zeros(3,no-2);
    d = zeros(2,no-2);
    
    temp1(:,1) = funker(X1(1),X1(4));
    temp1(:,2) = funker(X2(1),X2(4));
    
    temp2(:,1) = (rmc.*cos(azmc))+xr;
    temp2(:,2) = (rmc.*sin(azmc))+yr;
    
    for i = 1:no-2
        temp3(:,i) = funker(temp2(i,1),temp2(i,2));
    end
    
    for j = 1:2
        for k = 1:no-2
            d(j,k) = fundist2(temp1(1,j),temp1(2,j),temp1(3,j),...
                              temp3(1,k),temp3(2,k),temp3(3,k));
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
end

% Distance function
function dist2 = fundist2(a,b,c,d,e,f)
    dist2 = sqrt((d-a)^2+(e-b)^2+(f-c)^2);
end

function kerval = funker(x,y)
    kerval = [x^2; y^2; sqrt(2)*x*y];
end
function [disp,vel,accl] = parabolicEqn( n, M,invM,C,K,F,t,x0,v0)

% 0 - the constant-average accelaration method (stable)
% 1 - the linear accelaration method (conditionally stable)
% 2 - the central difference method (conditionally stable)
% 3 - the Galerkin method (stable)
% 4 - the backward difference method (stable)
    
    switch n
        case 0
            a = 0.5;
            b = 0.5;
        case 1
            a = 0.5;
            b = 1/3;
        case 2
            a = 0.5;
            b = 0;
        case 3
            a = 3/2;
            b = 8/5;
        case 4
            a = 3/2;
            b = 2;
    end
    
    nD = length(M);%  Dimension of problem
    disp = zeros(nD,length(t));
    vel = zeros(nD,length(t));
    accl = zeros(nD,length(t));
    disp(:,1) = x0;
    vel(:,1) = v0;
    
    
    accl(:,1) = invM * ( F(:,1) - K * x0 - C * v0);
    
    for i = 2:length(t)
        if i == size(F,2) + 1
            dF = -F(:,size(F,2));
        elseif i > size(F,2) + 1
            dF = zeros(nD,1);
        else
            dF = F(:,i) - F(:,i-1);
        end
        dT = t(i) - t(i-1);
        KK = (2/(b * dT^2)) * M + (2 * a/(b * dT)) * C + K;
        FF = dF + ( (2/(b * dT)) * M + (2*a/b) * C) * vel(:,i-1) + ( (1/b) * M  + dT * ( 1 - a/b) *  C ) * accl(:,i-1);
        
        dU = KK\FF;
        disp(:,i) = disp(:,i-1) + dU;
        vel(:,i) = vel(:,i-1) + dT * (1 - a/b) * accl(:,i-1) + (2*a/(b*dT)) * dU  - (2*a/b) * vel(:,i-1);
        accl(:,i) = accl(:,i-1) + (2/(b*dT^2))*dU - (2/(b*dT))*vel(:,i-1) - (1/b) * accl(:,i-1);
    end
end


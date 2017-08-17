function C = fluid_radiating_matrix( localization_matrix,node1,node2,node3,node4, side_number, c, thickness )
    gauss_X = [-0.7745966692	0	0.7745966692];
    gauss_W = [0.5555555556	0.8888888889	0.5555555556];
    switch side_number
        case 1
            ee = gauss_X;
            nn = [-1,-1,-1];
            den = [1;0];
        case 2
            ee = [1,1,1];
            nn = gauss_X;
            den = [0;1];
        case 3
            ee = gauss_X;
            nn = [1,1,1];
            den = [1;0];
        case 4
            ee = [-1,-1,-1];
            nn = gauss_X;
            den = [0;1];
    end
    C = zeros(4,4);
    for i=1:3
        [J,N,~] = jacobian_matrix(node1, node2, node3, node4,ee(i),nn(i));
        dL = sqrt((J*den)'*(J*den));
        C = C + thickness*(1/c)*gauss_W(i)*(N'*N)*dL;
    end
    C = localization_matrix'*C*localization_matrix;
    return
end
function [J,N,dN] = jacobian_matrix(node1, node2, node3, node4,e,n)
    X = [node1; node2; node3; node4];
    P = 0.25*[-(1-n),(1-n),(1+n),-(1+n);
                -(1-e),-(1+e),(1+e),(1-e)];
    J = P*X;
    dN = inv(J)*P;
    N = [(1-e)*(1-n)/4,(1+e)*(1-n)/4,(1+e)*(1+n)/4,(1-e)*(1+n)/4];
end

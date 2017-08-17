function [ S,H ] = fluid_quad_element( node1, node2, node3, node4, c, thickness )
    gauss_X = [-0.7745966692	0	0.7745966692];
    gauss_W = [0.5555555556	0.8888888889	0.5555555556];
    S = zeros(4,4);
    H = zeros(4,4);
    for i=1:3
        for j=1:3
            [J,N,dN] = jacobian_matrix(node1, node2, node3, node4,gauss_X(i),gauss_X(j));
            H = H + thickness*gauss_W(i)*gauss_W(j)*det(J)*(dN'*dN);
            S = S + thickness*gauss_W(i)*gauss_W(j)*det(J)*(N'*N)/c^2;
        end
    end
end
function [J,N,dN] = jacobian_matrix(node1, node2, node3, node4,e,n)
    X = [node1; node2; node3; node4];
    P = 0.25*[-(1-n),(1-n),(1+n),-(1+n);
                -(1-e),-(1+e),(1+e),(1-e)];
    J = P*X;
    dN = inv(J)*P;
    N = [(1-e)*(1-n)/4,(1+e)*(1-n)/4,(1+e)*(1+n)/4,(1-e)*(1+n)/4];
end
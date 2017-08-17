function [ M, K ] = solid_quadrilateral_element( node1, node2, node3, node4, density, modulus_elasticity, poisson_ratio, thickness, TYPE )
    gauss_X = [-0.7745966692	0	0.7745966692];
    gauss_W = [0.5555555556	0.8888888889	0.5555555556];
    if strcmp(TYPE,'strain')
        D = [1-poisson_ratio,poisson_ratio,0;
            poisson_ratio,1-poisson_ratio,0;
            0,0,0.5*(1-2*poisson_ratio)];
        D = (modulus_elasticity/((1+poisson_ratio)*(1-2*poisson_ratio))) * D;
    elseif strcmp(TYPE,'stress')
        D = [1.0,poisson_ratio,0.0;
            poisson_ratio,1.0,0.0;
            0.0,0.0,0.5*(1-poisson_ratio)];
        D = (modulus_elasticity/(1-poisson_ratio^2))*D;
    end
    K = zeros(8,8);
    M = zeros(8,8);
    weight = 0;
    for i=1:3
        for j=1:3
            [Jdet,B] = jacobian_matrix(gauss_X(i),gauss_X(j),node1, node2, node3, node4);
            K = K + thickness*gauss_W(i)*gauss_W(j)*Jdet*B'*D*B;
            M = M + thickness*density*gauss_W(i)*gauss_W(j)*mass_interpolation(gauss_X(i),gauss_X(j))*Jdet;
            weight = weight + thickness*density*gauss_W(i)*gauss_W(j)*Jdet;
        end
    end
end
function M = mass_interpolation(s,t)
    N = [(1-s)*(1-t)/4,0.0,(1+s)*(1-t)/4,0.0,(1+s)*(1+t)/4,0.0,(1-s)*(1+t)/4,0.0;
        0.0,(1-s)*(1-t)/4,0.0,(1+s)*(1-t)/4,0.0,(1+s)*(1+t)/4,0.0,(1-s)*(1+t)/4];
    M = N'*N;
end
function [Jdet,B] = jacobian_matrix(e,n,node1, node2, node3, node4)
    X = [node1(1),node2(1),node3(1),node4(1)];
    Y = [node1(2),node2(2),node3(2),node4(2)];
    dE = 0.25*[-(1-n),(1-n),(1+n),-(1+n)];
    dN = 0.25*[-(1-e),-(1+e),(1+e),(1-e)];
    J11 = sum(dE.*X);
    J12 = sum(dE.*Y);
    J21 = sum(dN.*X);
    J22 = sum(dN.*Y);
    Jdet = J11*J22-J12*J21;
    dNx = (J22*dE-J12*dN)/Jdet;
    dNy = (-J21*dE+J11*dN)/Jdet;
    B = [dNx(1),0,dNx(2),0,dNx(3),0,dNx(4),0;0,dNy(1),0,dNy(2),0,dNy(3),0,dNy(4);dNy(1),dNx(1),dNy(2),dNx(2),dNy(3),dNx(3),dNy(4),dNx(4)];
    return
end
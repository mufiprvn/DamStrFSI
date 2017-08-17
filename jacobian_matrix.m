function [J,Jdet,B] = jacobian_matrix(e,n,node1, node2, node3, node4)
    X = [node1(1),node2(1),node3(1),node4(1)];
    Y = [node1(2),node2(2),node3(2),node4(2)];
    dE = 0.25*[-(1-n),(1-n),(1+n),-(1+n)];
    dN = 0.25*[-(1-e),-(1+e),(1+e),(1-e)];
    J = [dE;dN]*[X',Y'];
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


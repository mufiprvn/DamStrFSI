function Q = Qmatrix( Fnode1, Fnode2, Fnode3, Fnode4, Snode1, Snode2, Snode3, Snode4, Bnode1, Bnode2,normal)
    if size(normal,1) == 1
        normal = normal';
    end
    BC_normal = (Bnode2-Bnode1)./sqrt(sum((Bnode2-Bnode1).^2));
    gauss_X = [-0.7745966692	0	0.7745966692];
    gauss_W = [0.5555555556	0.8888888889	0.5555555556];
    if BC_normal(1) ~= 0
        m = BC_normal(2)/BC_normal(1);
        c = Bnode1(2)-m*Bnode1(1);
        X = Bnode1(1)*(1-gauss_X)+Bnode2(1)*gauss_X;
        Y = m*X+c;
        jacobian = sqrt(1+m^2)*(Bnode2(1)-Bnode1(1));
    else
        X = [Bnode1(1),Bnode1(1),Bnode1(1)];
        Y = Bnode1(2)*(1-gauss_X)+Bnode2(2)*gauss_X;
        jacobian = Bnode2(2)-Bnode1(2);
    end
    Q = zeros(8,4);
    for i=1:3
        N = cartesian_interpolation( X(i), Y(i), Snode1, Snode2, Snode3, Snode4 );
        Ns = [N(1),0.0,N(2),0.0,N(3),0.0,N(4),0.0;
            0.0,N(1),0.0,N(2),0.0,N(3),0.0,N(4)];
        Nf = cartesian_interpolation( X(i), Y(i), Fnode1, Fnode2, Fnode3, Fnode4 );
        Q = Q + gauss_W(i)*jacobian*( (Ns'*normal)*Nf );
    end
%     if  isinf(m)
%         ss
%     else
%         for i=1:3
%             N = cartesian_interpolation( X(i), Y(i), Snode1, Snode2, Snode3, Snode4 );
%             Ns = [N(1),0.0,N(2),0.0,N(3),0.0,N(4),0.0;
%                 0.0,N(1),0.0,N(2),0.0,N(3),0.0,N(4)];
%             Nf = cartesian_interpolation( X(i), Y(i), Fnode1, Fnode2, Fnode3, Fnode4 );
%         end
%     end
    
end
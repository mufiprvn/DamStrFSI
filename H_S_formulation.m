function [ H, S ] = H_S_formulation( nodes, elements, boundary_nodes, free_nodes, c, thickness )
    boundary_dof = boundary_nodes;
    free_dof = free_nodes ;
    Nf = length(free_dof);
    H = zeros(Nf,Nf);
    S = zeros(Nf,Nf);
    for i=1:size(elements)
        n1 = elements(i,1);
        n2 = elements(i,2);
        n3 = elements(i,3);
        n4 = elements(i,4);
        global_all_dof = [n1,n2, n3,n4];
        global_free_dof = index_of_vector_in_vector(free_dof,setdiff(global_all_dof,boundary_dof));
        local_free_dof = index_of_vector_in_vector(global_all_dof,setdiff(global_all_dof,boundary_dof));
        [ Se,He ] = fluid_quad_element( nodes(n1,:), nodes(n2,:), nodes(n3,:), nodes(n4,:), c, thickness );
        H(global_free_dof,global_free_dof) = H(global_free_dof,global_free_dof) + He(local_free_dof,local_free_dof);
        S(global_free_dof,global_free_dof) = S(global_free_dof,global_free_dof) + Se(local_free_dof,local_free_dof);
    end
end
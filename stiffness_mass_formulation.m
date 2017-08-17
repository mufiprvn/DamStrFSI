function [ Kff,Kbf,Mff,Mbf ] = stiffness_mass_formulation( nodes, elements, boundary_nodes, free_nodes, density, thickness, modulus_elasticity, poisson_ratio, plane_stress_strain )
    boundary_dof = convert_nodes_dof( boundary_nodes );
    free_dof = convert_nodes_dof( free_nodes );
    Nb = length(boundary_dof);
    Nf = length(free_dof);
    Kff = zeros(Nf,Nf);
    Kbf = zeros(Nb,Nf);
    Mff = zeros(Nf,Nf);
    Mbf = zeros(Nb,Nf);
    for i=1:size(elements)
        n1 = elements(i,1);
        n2 = elements(i,2);
        n3 = elements(i,3);
        n4 = elements(i,4);
        global_all_dof = [2*n1-1,2*n1,2*n2-1,2*n2,2*n3-1,2*n3,2*n4-1,2*n4];
        global_free_dof = index_of_vector_in_vector(free_dof,setdiff(global_all_dof,boundary_dof));
        local_free_dof = index_of_vector_in_vector(global_all_dof,setdiff(global_all_dof,boundary_dof));
        global_boundary_dof = index_of_vector_in_vector(boundary_dof, setdiff(global_all_dof,free_dof));
        local_boundary_dof = setdiff(1:8,local_free_dof);
        [ Me, Ke ] = solid_quadrilateral_element( nodes(n1,:), nodes(n2,:), nodes(n3,:), nodes(n4,:), density, modulus_elasticity, poisson_ratio, thickness, plane_stress_strain );
%        [ Me , Ke, status ] = element_matrix( nodes(n1,:), nodes(n2,:),
%        nodes(n3,:), density, thickness, modulus_elasticity, poisson_ratio, plane_stress_strain);
        Kff(global_free_dof,global_free_dof) = Kff(global_free_dof,global_free_dof) + Ke(local_free_dof,local_free_dof);
        Kbf(global_boundary_dof,global_free_dof) = Kbf(global_boundary_dof,global_free_dof) + Ke(local_boundary_dof,local_free_dof);
        Mff(global_free_dof,global_free_dof) = Mff(global_free_dof,global_free_dof) + Me(local_free_dof,local_free_dof);
        Mbf(global_boundary_dof,global_free_dof) = Mbf(global_boundary_dof,global_free_dof) + Me(local_boundary_dof,local_free_dof);
    end
end


function C = fluid_radiating_matrix_formulation( nodes, elements, free_dof, boundary_dof, radiating_elements, c, thickness )
    figure
    C = zeros(length(free_dof),length(free_dof));
    for i=1:length(radiating_elements)
        F = elements(radiating_elements(i),:);
        X = [nodes(F(1),1), nodes(F(2),1), nodes(F(3),1), nodes(F(4),1), nodes(F(1),1)];
        Y = [nodes(F(1),2), nodes(F(2),2), nodes(F(3),2), nodes(F(4),2), nodes(F(1),2)];
        for p=1:4
            normal = normal_to_line( [X(p),Y(p)], [X(p+1),Y(p+1)],'R');
            if (normal(2) == 0) && (normal(1) < 0)
                side = p;
                plot([X(p),X(p+1)],[Y(p),Y(p+1)],'g')
                break
            end
        end
        global_all_dof = [F(1), F(2), F(3), F(4)];
        global_free_dof = index_of_vector_in_vector(free_dof,setdiff(global_all_dof,boundary_dof));
        local_free_dof = index_of_vector_in_vector(global_all_dof,setdiff(global_all_dof,boundary_dof));
        Ce = fluid_radiating_matrix( eye(4),nodes(F(1),:),nodes(F(2),:),nodes(F(3),:),nodes(F(4),:), side, c, thickness );
        C(global_free_dof,global_free_dof) = C(global_free_dof,global_free_dof) + Ce(local_free_dof,local_free_dof);
    end
end


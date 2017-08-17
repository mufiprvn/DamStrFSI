function Q = Qmatrix_formulation( nodes, fluid_nodes, fluid_free_nodes, fluid_boundary_nodes, dam_nodes, dam_free_nodes, dam_boundary_nodes, interface, dam_side )
    dam_free_dof = convert_nodes_dof( dam_free_nodes );
    dam_boundary_dof = convert_nodes_dof( dam_boundary_nodes );
    fluid_free_dof = fluid_free_nodes;
    fluid_boundary_dof = fluid_boundary_nodes;
    Q = zeros(length(dam_free_dof), length(fluid_free_dof));
    figure;
    for i=1:length(interface)
        Fs = interface{i}.dam_nodes;
        dam_global_dof_all = convert_nodes_dof( Fs );
        dam_global_dof_free = index_of_vector_in_vector(dam_free_dof, setdiff(dam_global_dof_all,dam_boundary_dof));
        dam_local_dof_free = index_of_vector_in_vector(dam_global_dof_all, setdiff(dam_global_dof_all,dam_boundary_dof));
        Ff = interface{i}.fluid_nodes;
        fluid_global_dof_all = Ff;
        fluid_global_dof_free = index_of_vector_in_vector(fluid_free_dof, setdiff(fluid_global_dof_all,fluid_boundary_dof));
        fluid_local_dof_free = index_of_vector_in_vector(fluid_global_dof_all, setdiff(fluid_global_dof_all,fluid_boundary_dof));
        Bnode1 = nodes(interface{i}.interface_line(1),:);
        Bnode2 = nodes(interface{i}.interface_line(2),:);
        normal =  normal_to_line( Bnode1, Bnode2, dam_side );
        Fnode1 = fluid_nodes(Ff(1),:); 
        Fnode2 = fluid_nodes(Ff(2),:); 
        Fnode3 = fluid_nodes(Ff(3),:);
        Fnode4 = fluid_nodes(Ff(4),:);
        Snode1 = dam_nodes(Fs(1),:);
        Snode2 = dam_nodes(Fs(2),:); 
        Snode3 = dam_nodes(Fs(3),:);
        Snode4 = dam_nodes(Fs(4),:);
        Qc = Qmatrix( Fnode1, Fnode2, Fnode3, Fnode4, Snode1, Snode2, Snode3, Snode4, Bnode1, Bnode2,normal);
        Q(dam_global_dof_free,fluid_global_dof_free) = Qc(dam_local_dof_free,fluid_local_dof_free);
        temp = cross([normal,0],[0,1,0]);
        if temp(3) < 0
            temp_node = Bnode1;
            Bnode1 = Bnode2;
            Bnode2 = temp_node;
            normal = -normal;
        end
        plot([Bnode1(1),Bnode2(1)],[Bnode1(2),Bnode2(2)],'r')
        CC = Bnode2 + 1*normal;
        hold all
        plot([Bnode2(1),CC(1)],[Bnode2(2),CC(2)],'g')
    end
    axis equal
end


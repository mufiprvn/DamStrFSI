%%
clear all; close all; clc;
%% Material propeties
modulus_elasticity = 5000*sqrt(30)*10^8; %Pa
poisson_ratio = 0.2;
density = 2500; % kg/m^3
c = 1481; %m/s
density_water = 1000; %kg/m^3
%% reading mesh file
load accl;
upstream_boundary = -200000;
fname = 'Mesh1.msh';
[ nodes, elements, nodes_cat, elements_cat ] = gmsh_processor( fname);
nodes = nodes/1000;
%% Processing data
% Categorize fluid and solids nodes
if (sum(nodes(nodes_cat{3}{1})) < 0)
    fluid_nodes = nodes_cat{3}{1};
    solid_nodes = nodes_cat{3}{2};
else
    solid_nodes = nodes_cat{3}{1};
    fluid_nodes = nodes_cat{3}{2};
end
% Categorize fluid and solids elements
if sum(reshape(elements{3}.nodes(elements_cat{3}{1},:),length(elements_cat{3}{1})*4,1) == 9)
    fluid_elements = elements_cat{3}{1};
    solid_elements = elements_cat{3}{2};
else
    fluid_elements = elements_cat{3}{2};
    solid_elements = elements_cat{3}{1};
end
F_N = fluid_nodes;
S_N = solid_nodes;
F_E = fluid_elements;
S_E = solid_elements;
%% Solids nodes and elements
% finding boundary of dams, fixity nodes
for i=1:length(nodes_cat{1})
    S1 = sum(nodes_cat{1}{i} == 1);
    S2 = sum(nodes_cat{1}{i} == 2);
    if ((S1+S2) == 2)
        index = i;
        break
    end
end
solid_boundary_nodes = index_of_vector_in_vector(solid_nodes,nodes_cat{1}{index} );
solid_free_nodes = setdiff(1:length(solid_nodes),solid_boundary_nodes);
% solid_nodes
solid_elements2 = zeros(length(solid_elements),4);
for i=1:length(solid_elements)
    F = elements{3}.nodes(solid_elements(i),:);
    solid_elements2(i,:) = index_of_vector_in_vector(solid_nodes, F );
end
solid_elements = solid_elements2;
solid_nodes = nodes(solid_nodes,:);

figure
axis equal
hold all
for i=1:size(solid_elements,1)
    F = solid_elements(i,:);
    XX = [solid_nodes(F(1),1),solid_nodes(F(2),1),solid_nodes(F(3),1),solid_nodes(F(4),1),solid_nodes(F(1),1)];
    YY = [solid_nodes(F(1),2),solid_nodes(F(2),2),solid_nodes(F(3),2),solid_nodes(F(4),2),solid_nodes(F(1),2)];
    plot(XX,YY,'c')
end
scatter(solid_nodes(:,1),solid_nodes(:,2))
scatter(solid_nodes(solid_boundary_nodes,1),solid_nodes(solid_boundary_nodes,2), 'r*')
[ Kff,Kbf,Mff,Mbf ] = stiffness_mass_formulation( solid_nodes, solid_elements, solid_boundary_nodes, solid_free_nodes, density, 1.0, modulus_elasticity, poisson_ratio, 'strain' );
Cs = zeros(size(Kff));
% [e_vect,e_val] = eig(Kff,Mff);
% [a,b] = sort(diag(e_val));
% for i=1:6
%     figure
%     modes = postprocessing( 10, solid_nodes, solid_elements, solid_free_nodes, solid_boundary_nodes, e_vect(:,b(i)) );
% end
%% Fluid nodes and elements
for i=1:length(nodes_cat{1})
    S1 = sum( nodes_cat{1}{i} == 6 );
    S2 = sum( nodes_cat{1}{i} == 8 );
    if  (S1+S2) == 2
        fluid_surf_dof = nodes_cat{1}{i};
    end
end
fluid_surf_nodes = index_of_vector_in_vector(fluid_nodes, fluid_surf_dof );
fluid_free_nodes = setdiff(1:length(fluid_nodes),fluid_surf_nodes);

% Finding elements in radiating elements
for i=1:length(elements_cat{1})
    S1 = sum( nodes_cat{1}{i} == 9 );
    S2 = sum( nodes_cat{1}{i} == 8 );
    if  (S1+S2) == 2
        fluid_radiating_nodes = nodes_cat{1}{i};
        fluid_radiating_line_elements = elements_cat{1}{i};
        break
    end
end

fluid_elements2 = zeros(length(fluid_elements),4);
fluid_upstream_elements = zeros(1,length(fluid_radiating_line_elements));
count = 1;
for i=1:length(fluid_elements)
    F = elements{3}.nodes(fluid_elements(i),:);
    sum1 = 0;
    for j=1:4
        sum1 = sum1 + sum(fluid_radiating_nodes==F(j));
    end
    fluid_elements2(i,:) = index_of_vector_in_vector(fluid_nodes, F );
    if sum1 == 2
        fluid_upstream_elements(count) = i;
        count = count + 1;
    end
end
fluid_elements = fluid_elements2;
fluid_nodes = nodes(fluid_nodes,:);
[ H, S ] = H_S_formulation( fluid_nodes, fluid_elements, fluid_surf_nodes, fluid_free_nodes, c, 1.0 );
Cf = fluid_radiating_matrix_formulation( fluid_nodes, fluid_elements, fluid_free_nodes, fluid_surf_nodes, fluid_upstream_elements, c, 1.0 );

figure
axis equal
hold all
for i=1:size(fluid_elements,1)
    F = fluid_elements(i,:);
    XX = [fluid_nodes(F(1),1),fluid_nodes(F(2),1),fluid_nodes(F(3),1),fluid_nodes(F(4),1),fluid_nodes(F(1),1)];
    YY = [fluid_nodes(F(1),2),fluid_nodes(F(2),2),fluid_nodes(F(3),2),fluid_nodes(F(4),2),fluid_nodes(F(1),2)];
    plot(XX,YY,'c')
end
for i=1:length(fluid_upstream_elements)
    F = fluid_elements(fluid_upstream_elements(i),:);
    XX = [fluid_nodes(F(1),1),fluid_nodes(F(2),1),fluid_nodes(F(3),1),fluid_nodes(F(4),1),fluid_nodes(F(1),1)];
    YY = [fluid_nodes(F(1),2),fluid_nodes(F(2),2),fluid_nodes(F(3),2),fluid_nodes(F(4),2),fluid_nodes(F(1),2)];
    plot(XX,YY,'r')
end
scatter(fluid_nodes(:,1),fluid_nodes(:,2))
scatter(fluid_nodes(fluid_surf_nodes,1),fluid_nodes(fluid_surf_nodes,2), 'r*')

%% Coupling matrix
for i=1:length(nodes_cat{1})
    S1 = sum( nodes_cat{1}{i} == 1 );
    S2 = sum( nodes_cat{1}{i} == 7 );
    S3 = sum( nodes_cat{1}{i} == 6 );
    if  (S1+S2) == 2
        index1 = i;
    end
    if  (S3+S2) == 2
        index2 = i;
    end
end
interface = cell(1,length(elements_cat{1}{index1})+length(elements_cat{1}{index2}));
for i=1:length(elements_cat{1}{index1})
    interface{i}.interface_line = [nodes_cat{1}{index1}(i),nodes_cat{1}{index1}(i+1)];
    for j=1:length(F_E)
        F = elements{3}.nodes(F_E(j),:);
        S1 = sum(F == elements{1}.nodes(elements_cat{1}{index1}(i),1));
        S2 = sum(F == elements{1}.nodes(elements_cat{1}{index1}(i),2));
        if (S1+S2) == 2
            interface{i}.fluid_nodes = fluid_elements(j,:);
            break;
        end
    end
    for j=1:length(S_E)
        F = elements{3}.nodes(S_E(j),:);
        S1 = sum(F == elements{1}.nodes(elements_cat{1}{index1}(i),1));
        S2 = sum(F == elements{1}.nodes(elements_cat{1}{index1}(i),2));
        if (S1+S2) == 2
            interface{i}.dam_nodes = solid_elements(j,:);
            break;
        end
    end
end

for i=1:length(elements_cat{1}{index2})
    interface{length(elements_cat{1}{index1})+i}.interface_line = [nodes_cat{1}{index2}(i),nodes_cat{1}{index2}(i+1)];
    for j=1:length(F_E)
        F = elements{3}.nodes(F_E(j),:);
        S1 = sum(F == elements{1}.nodes(elements_cat{1}{index2}(i),1));
        S2 = sum(F == elements{1}.nodes(elements_cat{1}{index2}(i),2));
        if (S1+S2) == 2
            interface{length(elements_cat{1}{index1})+i}.fluid_nodes = fluid_elements(j,:);
            break;
        end
    end
    for j=1:length(S_E)
        F = elements{3}.nodes(S_E(j),:);
        S1 = sum(F == elements{1}.nodes(elements_cat{1}{index2}(i),1));
        S2 = sum(F == elements{1}.nodes(elements_cat{1}{index2}(i),2));
        if (S1+S2) == 2
            interface{length(elements_cat{1}{index1})+i}.dam_nodes = solid_elements(j,:);
            break;
        end
    end
end

Q = Qmatrix_formulation( nodes, fluid_nodes, fluid_free_nodes, fluid_surf_nodes, solid_nodes, solid_free_nodes, solid_boundary_nodes, interface, 'R' );

figure
axis equal
hold all
for i=1:length(interface)
    F = interface{i}.fluid_nodes;
    XX = [fluid_nodes(F(1),1),fluid_nodes(F(2),1),fluid_nodes(F(3),1),fluid_nodes(F(4),1),fluid_nodes(F(1),1)];
    YY = [fluid_nodes(F(1),2),fluid_nodes(F(2),2),fluid_nodes(F(3),2),fluid_nodes(F(4),2),fluid_nodes(F(1),2)];
    plot(XX,YY,'c')
    F = interface{i}.dam_nodes;
    XX = [solid_nodes(F(1),1),solid_nodes(F(2),1),solid_nodes(F(3),1),solid_nodes(F(4),1),solid_nodes(F(1),1)];
    YY = [solid_nodes(F(1),2),solid_nodes(F(2),2),solid_nodes(F(3),2),solid_nodes(F(4),2),solid_nodes(F(1),2)];
    plot(XX,YY,'r')
end
%% Final Matrix
MM = [Mff,zeros(2*length(solid_free_nodes), length(fluid_free_nodes));density_water*Q',S];
CC = [Cs,zeros(2*length(solid_free_nodes), length(fluid_free_nodes));zeros(length(fluid_free_nodes),2*length(solid_free_nodes)), Cf];
KK = [Kff,-Q;zeros(length(fluid_free_nodes),2*length(solid_free_nodes)), H];
%% Force Transfer function
Is = zeros(2*length(solid_free_nodes),2);
Is(1:2:end,1) = 1;
IsUg = Is*[1;0];
F = [-Mff*IsUg;density*Q'*IsUg];
x0 = zeros(size(F));
v0 = zeros(size(F));
FF = 9.81*F*accl(:,2)';
time = accl(:,1);
[disp,vel,accl] = parabolicEqn( 0, MM,inv(MM),CC,KK,FF,accl(:,1),x0,v0);
matfile = strcat(fname(1:end-4));
save(matfile)


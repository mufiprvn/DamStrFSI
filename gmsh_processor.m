function [ nodes, elements, nodes_cat, elements_cat ] = gmsh_processor( fname)
    %% gmsh node number
    gmsh_nodes = [2	3	4	4	8	6	5	3	6	9	10	27	18	14	1	8	20	15	13	9	10	12	15	15	21	4	5	6	20	35	56];
    %%
    fid = fopen(fname,'r');
    nodes_start = 0;
    elements_start = 0;
    elements = cell(1,31);
    for i=1:31
        elements{i}.num_nodes = gmsh_nodes(i);
        elements{i}.nodes = zeros(1,gmsh_nodes(i));
        elements{i}.element_num = [];
        elements{i}.num_tags = [];
        elements{i}.tags = zeros(1,5);
    end
    while ~feof(fid)
        %% Reading nodes
        str = fgets(fid);
        if (strcontains( '$Nodes' , str ) > 0)
            nodes_start = 1;
            str = fgets(fid);
            num_nodes = str2double(str);
            nodes = zeros(num_nodes,2);
            count = 1;
            continue
        end
        if (strcontains( '$EndNodes' , str ) > 0)
            nodes_start = 0;
            continue
        end
        if (nodes_start > 0)
            A = sscanf(str,'%f %f %f');
            nodes(count,1) = A(2);
            nodes(count,2) = A(3);
            count = count +1;
        end
        %% Reding elements
        if (strcontains( '$Elements' , str ) > 0)
            elements_start = 1;
            str = fgets(fid);
            display(str);
%            num_elements = str2double(str);
            count = 1;
            continue
        end
        if (strcontains( '$EndElements' , str ) > 0)
            elements_start = 0;
            continue
        end
        if (elements_start > 0)
            A = sscanf(str,'%f');
            elements{A(2)}.element_num(end+1) = A(1);
            elements{A(2)}.num_tags(end+1) = A(3);
            elements{A(2)}.tags(end+1,1:A(3)) = A(4:3+A(3));
            elements{A(2)}.nodes(end+1,:) = A(4+A(3):end);
        end
    end
    fclose(fid);
    for i=1:31
        elements{i}.nodes = elements{i}.nodes(2:end,:);
        elements{i}.tags = elements{i}.tags(2:end,:);
    end
    
    %% Nodes scatter
%    scatter(nodes(:,1),nodes(:,2));
    figure;
    hold all;
    axis equal;
    %% Plot elements
    elem_color = ['r','g'];
    for i=[3,1]
        elements{i}.categories = categorize_vector( elements{i}.tags(:,2) );
    end
    line_color = rand(length(elements{1}.categories),3);
    for i=1:size(elements{3}.nodes,1)
        F = elements{3}.nodes(i,:);
        XX = [nodes(F(1),1),nodes(F(2),1),nodes(F(3),1),nodes(F(4),1),nodes(F(1),1)];
        YY = [nodes(F(1),2),nodes(F(2),2),nodes(F(3),2),nodes(F(4),2),nodes(F(1),2)];
        [~,b] = max(elements{3}.categories == elements{3}.tags(i,2));
        plot(XX,YY,elem_color(b))
    end
    for i=1:size(elements{1}.nodes,1)
        F = elements{1}.nodes(i,:);
        XX = [nodes(F(1),1),nodes(F(2),1)];
        YY = [nodes(F(1),2),nodes(F(2),2)];
        [~,b] = max(elements{1}.categories == elements{1}.tags(i,2));
        plot(XX,YY,'color',line_color(b,:), 'LineWidth', 3);
    end
    %% Processing the data
    % Categorize nodes and elements
    nodes_cat = cell(1,3);
    elements_cat = cell(1,3);
    for i=[1,3]
        nodes_cat{i} = cell(1,length(elements{i}.categories));
        elements_cat{i} = cell(1,length(elements{i}.categories));
        for j=1:length(elements{i}.categories)
            nodes_cat{i}{j} = [];
            elements_cat{i}{j} = [];
        end
    end
    PPP = [2,0,4];
    for p=[1,3]
        for i=1:size(elements{p}.nodes,1);
            [~,b] = max(elements{p}.categories == elements{p}.tags(i,2));
            elements_cat{p}{b}(end+1) = i;
            for k=1:PPP(p)
                if ~(sum(nodes_cat{p}{b} == elements{p}.nodes(i,k)))
                    nodes_cat{p}{b}(end+1) = elements{p}.nodes(i,k);
                end
            end
        end
    end
    cc = {'c','y*'};
    for i=1:length(elements{3}.categories)
        scatter(nodes(nodes_cat{3}{i},1),nodes(nodes_cat{3}{i},2), cc{i})
    end
end

% 
% 
% 
% 1
% 
%     2-node line. 
% 2
% 
%     3-node triangle. 
% 3
% 
%     4-node quadrangle. 
% 4
% 
%     4-node tetrahedron. 
% 5
% 
%     8-node hexahedron. 
% 6
% 
%     6-node prism. 
% 7
% 
%     5-node pyramid. 
% 8
% 
%     3-node second order line (2 nodes associated with the vertices and 1 with the edge). 
% 9
% 
%     6-node second order triangle (3 nodes associated with the vertices and 3 with the edges). 
% 10
% 
%     9-node second order quadrangle (4 nodes associated with the vertices, 4 with the edges and 1 with the face). 
% 11
% 
%     10-node second order tetrahedron (4 nodes associated with the vertices and 6 with the edges). 
% 12
% 
%     27-node second order hexahedron (8 nodes associated with the vertices, 12 with the edges, 6 with the faces and 1 with the volume). 
% 13
% 
%     18-node second order prism (6 nodes associated with the vertices, 9 with the edges and 3 with the quadrangular faces). 
% 14
% 
%     14-node second order pyramid (5 nodes associated with the vertices, 8 with the edges and 1 with the quadrangular face). 
% 15
% 
%     1-node point. 
% 16
% 
%     8-node second order quadrangle (4 nodes associated with the vertices and 4 with the edges). 
% 17
% 
%     20-node second order hexahedron (8 nodes associated with the vertices and 12 with the edges). 
% 18
% 
%     15-node second order prism (6 nodes associated with the vertices and 9 with the edges). 
% 19
% 
%     13-node second order pyramid (5 nodes associated with the vertices and 8 with the edges). 
% 20
% 
%     9-node third order incomplete triangle (3 nodes associated with the vertices, 6 with the edges) 
% 21
% 
%     10-node third order triangle (3 nodes associated with the vertices, 6 with the edges, 1 with the face) 
% 22
% 
%     12-node fourth order incomplete triangle (3 nodes associated with the vertices, 9 with the edges) 
% 23
% 
%     15-node fourth order triangle (3 nodes associated with the vertices, 9 with the edges, 3 with the face) 
% 24
% 
%     15-node fifth order incomplete triangle (3 nodes associated with the vertices, 12 with the edges) 
% 25
% 
%     21-node fifth order complete triangle (3 nodes associated with the vertices, 12 with the edges, 6 with the face) 
% 26
% 
%     4-node third order edge (2 nodes associated with the vertices, 2 internal to the edge) 
% 27
% 
%     5-node fourth order edge (2 nodes associated with the vertices, 3 internal to the edge) 
% 28
% 
%     6-node fifth order edge (2 nodes associated with the vertices, 4 internal to the edge) 
% 29
% 
%     20-node third order tetrahedron (4 nodes associated with the vertices, 12 with the edges, 4 with the faces) 
% 30
% 
%     35-node fourth order tetrahedron (4 nodes associated with the vertices, 18 with the edges, 12 with the faces, 1 in the volume) 
% 31
% 
%     56-node fifth order tetrahedron (4 nodes associated with the vertices, 24 with the edges, 24 with the faces, 4 in the volume) 

clear all; close all; clc;
%% 
water_level = 80000;
C = 1624.584;
D = 39014.400;
M = D/C;
load NODE;
load ELEMENT;
scatter(NODE(:,2),NODE(:,3))
hold all
NODE1 = [0.000	0.000];
NODE2 = [70189.344	0.000];
NODE3 = [22006.560	66446.400];
NODE4 = [16407.384	103022.400];
NODE5 = [1624.584	103022.400];
NODE6 = [1624.584	39014.400];
NODES = [NODE1;NODE2;NODE3;NODE4;NODE5;NODE6;NODE1];
%% Fluid nodes and elements
plot(NODES(:,1),NODES(:,2),'g')
FLUID_NODES = [];
DAM_NODES = [];
for i=1:size(NODE,1)
    if NODE(i,2) <= C
        if (NODE(i,3) <= water_level) && (NODE(i,2)<=(NODE(i,3)/M+1))
            scatter(NODE(i,2),NODE(i,3),'r*')
            FLUID_NODES = [FLUID_NODES;NODE(i,:)];
        end
    end
end
FLUID_ELEMENT = [];
for i=1:size(ELEMENT,1)
    F = ELEMENT(i,:);
    G = zeros(1,4);
    DD = 0;
    for j=1:length(F)
        [a,b] = max(FLUID_NODES(:,1) == F(j));
        DD = DD + a;
        G(j) = b;
    end
    if DD == 4
        FLUID_ELEMENT = [FLUID_ELEMENT;G];
%        X = [NODE(F(1),2),NODE(F(2),2),NODE(F(3),2),NODE(F(4),2),NODE(F(1),2)];
%        Y = [NODE(ELEMENT(i,1),3),NODE(ELEMENT(i,2),3),NODE(ELEMENT(i,3),3),NODE(ELEMENT(i,4),3),NODE(ELEMENT(i,1),3)];
%        plot(X,Y,'c');
    end
end
for i=1:size(FLUID_ELEMENT,1)
    F = FLUID_ELEMENT(i,:);
    X = [FLUID_NODES(F(1),2),FLUID_NODES(F(2),2),FLUID_NODES(F(3),2),FLUID_NODES(F(4),2),FLUID_NODES(F(1),2)];
    Y = [FLUID_NODES(F(1),3),FLUID_NODES(F(2),3),FLUID_NODES(F(3),3),FLUID_NODES(F(4),3),FLUID_NODES(F(1),3)];
    text(sum(X)/5,sum(Y)/5,num2str(i));
    plot(X,Y,'c');
end
%% Boundary
boundary_nodes = [];
for i=1:size(NODE,1)
    if (NODE(i,2) >= 0) && (NODE(i,2) <= C)
        if NODE(i,2) == C
            if (NODE(i,3) <= water_level) && (NODE(i,3) >= D)
                boundary_nodes = [boundary_nodes,i];
                scatter(NODE(i,2),NODE(i,3),'c*')
            end
        else
            if (NODE(i,2)>=(NODE(i,3)/M-1)) && (NODE(i,2)<=(NODE(i,3)/M+1))
                boundary_nodes = [boundary_nodes,i];
                scatter(NODE(i,2),NODE(i,3),'g','fill')
            end
        end
    end
end
%% DAM NODES and ELEMENTS
DAM_NODES = [];
for i=1:size(NODE,1)
    F1 = sum(FLUID_NODES(:,1) == NODE(i,1));
    F2 = sum(boundary_nodes == NODE(i,1));
    if (F1 < 1)
        scatter(NODE(i,2),NODE(i,3),'r','fill')
        DAM_NODES = [DAM_NODES;NODE(i,:)];
    end
    if F2 > 0
        scatter(NODE(i,2),NODE(i,3),'r','fill')
        DAM_NODES = [DAM_NODES;NODE(i,:)];
    end
end
DAM_ELEMENT = [];
for i=1:size(ELEMENT,1)
    F = ELEMENT(i,:);
    G = zeros(1,4);
    DD = 0;
    for j=1:length(F)
        [a,b] = max(DAM_NODES(:,1) == F(j));
        DD = DD + a;
        G(j) = b;
    end
    if DD == 4
        DAM_ELEMENT = [DAM_ELEMENT;G];
%        X = [NODE(F(1),2),NODE(F(2),2),NODE(F(3),2),NODE(F(4),2),NODE(F(1),2)];
%        Y = [NODE(ELEMENT(i,1),3),NODE(ELEMENT(i,2),3),NODE(ELEMENT(i,3),3),NODE(ELEMENT(i,4),3),NODE(ELEMENT(i,1),3)];
%        plot(X,Y,'c');
    end
end
for i=1:size(DAM_ELEMENT,1)
    F = DAM_ELEMENT(i,:);
    X = [DAM_NODES(F(1),2),DAM_NODES(F(2),2),DAM_NODES(F(3),2),DAM_NODES(F(4),2),DAM_NODES(F(1),2)];
    Y = [DAM_NODES(F(1),3),DAM_NODES(F(2),3),DAM_NODES(F(3),3),DAM_NODES(F(4),3),DAM_NODES(F(1),3)];
    plot(X,Y,'y');
end
%% Fluid Boundary Elements
FLUID_BOUNDARY_ELEMENT = [];
FLUID_UPSTREAM_ELEMENT = [];
for i=1:size(FLUID_ELEMENT,1)
    F = FLUID_ELEMENT(i,:);
    X = [FLUID_NODES(F(1),2),FLUID_NODES(F(2),2),FLUID_NODES(F(3),2),FLUID_NODES(F(4),2),FLUID_NODES(F(1),2)];
    Y = [FLUID_NODES(F(1),3),FLUID_NODES(F(2),3),FLUID_NODES(F(3),3),FLUID_NODES(F(4),3),FLUID_NODES(F(1),3)];
    sum = 0;
    sum2 = 0;
    for j=1:length(F)
        YY = FLUID_NODES(F(j),1);
        [a,b] = max(boundary_nodes == YY);
        sum = sum + a;
        if sum == 2
            FLUID_BOUNDARY_ELEMENT = [FLUID_BOUNDARY_ELEMENT,i];
            fill(X,Y,'r')
            break
        end
    end
    for j=1:length(F)
        YY = FLUID_NODES(F(j),2);
        c = (YY<=-200000);
        sum2 = sum2 + c;
        if sum2 == 2
            FLUID_UPSTREAM_ELEMENT = [FLUID_UPSTREAM_ELEMENT,i];
            fill(X,Y,'g')
            break
        end
    end
end
%% Dam Boundary Element
DAM_BOUNDARY_ELEMENT = [];
for i=1:size(DAM_ELEMENT,1)
    F = DAM_ELEMENT(i,:);
    X = [DAM_NODES(F(1),2),DAM_NODES(F(2),2),DAM_NODES(F(3),2),DAM_NODES(F(4),2),DAM_NODES(F(1),2)];
    Y = [DAM_NODES(F(1),3),DAM_NODES(F(2),3),DAM_NODES(F(3),3),DAM_NODES(F(4),3),DAM_NODES(F(1),3)];
    sum = 0;
    for j=1:length(F)
        YY = DAM_NODES(F(j),1);
        [a,b] = max(boundary_nodes == YY);
        sum = sum + a;
        if sum == 2
            DAM_BOUNDARY_ELEMENT = [DAM_BOUNDARY_ELEMENT,i];
            fill(X,Y,'y')
            break
        end
    end
end
DAM_FIX_NODES = [];
for i=1:size(DAM_NODES,1)
    if DAM_NODES(i,3) == 0
        DAM_FIX_NODES = [DAM_FIX_NODES,i];
        scatter(DAM_NODES(i,2),DAM_NODES(i,3),'g','fill')
    end
end
FLUID_SURF_NODES = [];
for i=1:size(FLUID_NODES,1)
    if FLUID_NODES(i,3) == 80000
        FLUID_SURF_NODES = [FLUID_SURF_NODES,i];
        scatter(FLUID_NODES(i,2),FLUID_NODES(i,3),'b','fill')
    end
end
%% Upstream Boundarty
save('models_processed')
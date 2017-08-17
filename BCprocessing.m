clear all; close all; clc;
load accl;
%time = accl(:,1);
time = linspace(0,0.1,21);
frequency = 10;
sin_wave = sin(2*pi*frequency*time);
accl = [time',sin_wave'];
%% BC Mesh 1
load Mesh_BC5.msh
key = 'Mesh_BC5';
fluid_pressure = [zeros(length(fluid_surf_nodes),length(time));disp(2*length(solid_free_nodes)+1:end,:)];
N = length(time);
F_NODE = [fluid_nodes(fluid_surf_nodes,:);fluid_nodes(fluid_free_nodes,:)];
tri = delaunay(F_NODE(:,1),F_NODE(:,2));

h = figure;
trisurf(tri, F_NODE(:,1),F_NODE(:,2) ,max(abs(fluid_pressure')),'EdgeColor','none','FaceColor','interp')
colorbar()
axis tight
axis equal
az = 0;
el = 90;
view(az, el);
figname = sprintf('%s_max.png',key);
saveas(h,figname);
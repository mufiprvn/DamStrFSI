clear all; close all; clc;
load Mesh1

[eig_vect, eig_val] = eig(KK,MM);
[a,b] = sort(diag(eig_val));
for i=1:4
    figure
    modes = postprocessing( 10^8, solid_nodes, solid_elements, solid_free_nodes, solid_boundary_nodes, eig_vect(1:2*length(solid_free_nodes),b(i)) );
    title(sprintf('%f',a(i)));
    axis equal
end
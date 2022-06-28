function [shape_array_3d] = create_3d_structure(array_2d,axial_height,axial_step)

if nargin~=3
    disp('Check input: Requires 2d array, axial height and step size');
    return
end

if axial_height == 0
    disp('No height selected, returning 2D shape');
    shape_array_3d = array_2d;
end

N=ceil(axial_height/axial_step);

size_2d=size(array_2d);
shape_array_3d = zeros([size_2d(1),size_2d(2),N]);

for ii=1:N
    shape_array_3d(:,:,ii) = array_2d;
end


end
function [weight_matrix] = get_normalised_weights(xpoints,ypoints)

    assert(length(xpoints)==length(ypoints), ... 
            ' X and Y must be the same length'); 
        
    %Load model
    addpath('R:\OMG\declan\PrintingSummaryProgram\Runfiles');
    addpath('R:\declan\PrintingSummaryProgram\Runfiles');
    load('slmLUT_191121.mat');
    
    xx=10e-6;
    yy=xx;
    dim=100;
%     [Xq,Yq] = meshgrid(linspace(coords(1),xlength+coords(1),dim)....
%             ,linspace(coords(2),ylength+coords(2),dim));
    [Xq,Yq] = meshgrid(linspace(-xx,xx,dim)...
            ,linspace(-yy,yy,dim));

    powerSurface =fittedmodel(Xq,Yq);
    powerSurface = smoothdata(powerSurface);
    P0Norm=max(max(powerSurface));
    powerSurface = powerSurface./P0Norm;
    
        
    %1D power scaling - assume x/y device symmetry
    rd = sqrt(xpoints.^2 + ypoints.^2);
    scale_values = fittedmodel(rd,zeros(size(rd)))./P0Norm;
    checkNAN = length(find(isnan(scale_values))>0);
    if sum(checkNAN) > 0
        disp('NAN Values found!')
        return
    end
    
    weight_matrix = scale_values;
    
    
%2D surface function, though WGS algorithm takes in 1D input 
%So assume equal x-y variation now
%     weight_matrix = ones([q,length(ypoints)]);
%     
%     for i=1:length(xpoints)
%         for j=1:length(ypoints)
%             
%             locX = xpoints(i);
%             locY = ypoints(j);
%             value = fittedmodel(locX,locY)./P0Norm;
%             weight_matrix(i,j) = value;
%             
%         end
%     end
%     
end
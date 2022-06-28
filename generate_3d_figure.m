function [output_success] = generate_3d_figure(bitmap)


if ndims(bitmap)==3
    %Number of elements
    xn=numel(bitmap(:,1,1));
    yn=numel(bitmap(1,:,1));
    zn=numel(bitmap(1,1,:));

    % %generate grid
    x=linspace(0,1,xn);
    y=linspace(0,1,yn);
    z=linspace(0,1,zn);

    for ii=1:zn
        slice_2d = bitmap(:,:,ii);
        [zrow,zcol]=find(slice_2d==1); %Take elevated points

        figure(12345),
        plot3(x(zrow),y(zcol),ones(size(zcol)).*z(ii),'+','color','k','linewidth',1.7); hold on;
        xlabel('X'); ylabel('Y');
        zlabel('Object Height'); 
        title('3D Bitmap Object');
        grid on;

    end
    
elseif ndims(bitmap)==2
    %Number of elements
    xn=numel(bitmap(:,1,1));
    yn=numel(bitmap(1,:,1));
        
    x=linspace(0,1,xn);
    y=linspace(0,1,yn);
    
    [zrow,zcol]=find(bitmap==1); %Take elevated points
    
    figure(12345),
    plot(x(zrow),y(zcol),'+','color','k','linewidth',1.7); hold on;
    xlabel('X'); ylabel('Y');
    title('2D Bitmap Object');
    grid on;
else
    disp(['ERROR: Input must be 2D or 3D. ']);
    output_success = false;
end
    
output_success = true;

end
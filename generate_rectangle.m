function [rect_bit,rect_polygon] = generate_rectangle(coords,width,height,dim)

%Input coordinates for bottom left corner
%shape is generated anti-clockwise

if numel(coords)~=2
    disp('First variable must contain x0,y0 input.');
    return
end

x0=coords(1); y0=coords(2);

xside = linspace(x0,x0+width,dim);
yside = linspace(y0,y0+height,dim);


rect_bit = ones([dim,dim]);

x_polygon = [xside(1),xside(1),xside(end),xside(end)];
y_polygon = [yside(end),yside(1),yside(1),yside(end)];

rect_polygon = polyshape(x_polygon,y_polygon);

%Plotting of grid
[X,Y] = meshgrid(xside,yside);
figure(),
plot(X,Y,'color','k'); hold on;
plot(rect_polygon,'facecolor','b'); grid on;
axis tight;
%plot(Y,X,'color','k'); hold on;
% xlim([(xside(1)-0.1*xside(1)), xside(end)+0.1*xside(end)]);
% xlim([(yside(1)-0.1*yside(1)), yside(end)+0.1*yside(end)]);


end
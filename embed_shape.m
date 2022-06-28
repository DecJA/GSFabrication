function [shapeScaled] =  embed_shape(polygon,xPoints,yPoints)

    %Input: polygon shape, xgrid(allRows), ygrid(allCols)
    
    shp = polygon;
    [xc,yc] = shp.centroid;
    xShift = xPoints - xc;
    yShift = yPoints - yc;

%     figure(),
%     plot(xgrid(allRows),ygrid(allCols),'+','color','r');
%     hold on;
%     plot(xShift,yShift,'+','color','b');
%     grid on;

    allDist = sqrt(xShift.^2 + yShift.^2);
    dx = max(abs(xShift)) - min(allDist);
    dy = max(abs(yShift)) - min(allDist);
    
    scaleRatio = (max(xShift-dx)/(max(xShift)));
    poly2 = scale(shp,scaleRatio,[xc,yc]);
    
    figure(7651),
    plot(shp); hold on;
    plot(poly2); grid on;
    
    shapeScaled = poly2;

    
end

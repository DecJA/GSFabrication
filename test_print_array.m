%load 2D pattern
dim =30;
coords=[0,0];
xlength=5e-6;
ylength=5e-6;
xl=linspace(0,xlength,dim);
yl=xl;

%shape = generate_triangle([-10e-6,-10e-6],xlength,ylength,dim);
[shape,~] = generate_rectangle([-10e-6,-10e-6],xlength,ylength,dim);
shape_bitmap_3d = create_3d_structure(shape,1e-6,200e-9);
shape2D=squeeze(shape_bitmap_3d(:,:,1)).';

[zrow,zcol]=find(shape2D==1); %Take elevated points
figure(12345),
plot(xl(zrow),yl(zcol),'+','color','k','linewidth',1.7); hold on;
xlabel('X'); ylabel('Y');
title('2D Bitmap Object');
grid on

load('slmLUT_191121.mat')
[Xq,Yq] = meshgrid(linspace(coords(1),xlength+coords(1),dim)....
        ,linspace(coords(2),ylength+coords(2),dim));
powerSurface =fittedmodel(Xq,Yq);%-Yq?
powerSurface = powerSurface./(max(max(powerSurface)));
powerSurface = smoothdata(powerSurface);
P0Max = min(min(powerSurface));
reducePowerBy=0;

lin_idx = find(shape2D>0);
cc = 1;
patternSet = {};
patternMulti = {};
spotLoc = {};
noSpots =5;

startPoints = zeros([1 noSpots]);
for xx=1:noSpots
  if xx==1
    startPoints(xx) = 1;
  else
    startPoints(xx) = ceil(length(lin_idx)*((xx-1)/noSpots)) +1;
  end
end

%%
counter=1;
xx=zeros([1 noSpots]);
yy=xx;
while cc<=length(lin_idx)
%   startPoints+counter
  for vv=1:noSpots
%     [row,col]=ind2sub(size(shape2D),lin_idx(cc+(vv-1))); 
%      [row,col]=ind2sub(size(shape2D),lin_idx(test_points(vv)));
    if startPoints(vv) +  counter>length(lin_idx)
      continue
    else
     [row,col]=ind2sub(size(shape2D),lin_idx(startPoints(vv)+counter));
     scalex = return_scale(xl(row));
     scaley = return_scale(yl(col));
     
      P0temp = powerSurface(row,col);
      powerFactor = 0;%reducePowerBy + ((1-P0Max) - (1-P0temp));
      [~,ditheredPattern] = modulateAmplitude(powerFactor,256,sz);
      
      weight(vv) = powerFactor;
      spotLoc{vv} = [scalex,-scaley];
      
      xx(vv) = scalex;
      yy(vv) = -scaley;
      counter=counter+1;
    end
  end
  cc=cc+noSpots;
  phase =  gs_spot_array(xx,yy,abs(amp)); %N spot pattern
  
%   ftImage = fftshift(fft2(fftshift(amp.*exp(1i.*phase))));
%   figure(1111),
%   subplot(1,2,1)
%   imagesc(phase);
%   subplot(1,2,2)
%   imagesc(abs(ftImage));
%   pause(0.5);
  slm.show(phase);
end

%%
for ii=1:dim
  for jj=1:dim
    
     loc = shape2D(ii,jj);
    if loc>0
      
      scalex = return_scale(xl(ii));
      scaley = return_scale(yl(jj));
      
      phase =  gs_spot_array(scalex,-scaley,abs(beam));
      
      P0temp = powerSurface(ii,jj);
      powerFactor = reducePowerBy + ((1-P0Max) - (1-P0temp));
      [~,ditheredPattern] = modulateAmplitude(powerFactor,128,sz);
      
      allPatterns = phase + ditheredPattern + zeroCorrection;
      
      slm.show(allPatterns);

    end
    
  end
end

%%
noSpots = 5;
startSpot = 1;
for tt=1:noSpots
  
      scalex = return_scale(xl(tt));
      scaley = return_scale(yl(tt));
      
      phase =  gs_spot_array(scalex,-scaley,abs(beam));
      
      P0temp = powerSurface(ii,jj);
      powerFactor = reducePowerBy + ((1-P0Max) - (1-P0temp));
      [~,ditheredPattern] = modulateAmplitude(powerFactor,128,sz);
      
      weight(tt) = powerFactor;
      
      patternSet{tt} = phase;
end
%%
lin_idx = find(shape2D>0);
cc = 1;
patternSet = {};
patternMulti = {};
spotLoc = {};
noSpots =10;


startPoints = zeros([1 ceil(length(lin_idx)/noSpots)]);
for xx=1:ceil(length(lin_idx)/noSpots)
  if xx==1
    startPoints(xx) = 1;
  else
    startPoints(xx) = startPoints(xx-1) + noSpots;
  end
end


savePatterns = zeros([length(startPoints),sz(1),sz(2)]);
xx=zeros([1 noSpots]);
yy=xx;
ccc=1;
cc=1;

for i=1:length(startPoints)
  counter=0;
  xx=zeros([1 noSpots]);
  yy=xx;
  
  if i==length(startPoints)
      noSpots = max(lin_idx) - max(startPoints) + 1;
  end
    
  for vv=1:noSpots
     [row,col]=ind2sub(size(shape2D),lin_idx(startPoints(i)+counter));
     scalex = return_scale(xl(row));
     scaley = return_scale(yl(col));
      
      weight(vv) = powerFactor;
      
      xx(vv) = scalex;
      yy(vv) = -scaley;
      counter=counter+1;
  end
  phase =  gs_spot_array(xx,yy,abs(amp)); %N spot pattern
  savePatterns(i,:,:) = phase;
  figure(1234),
  imagesc(phase);
  
  disp(['Found spot: ',num2str(i),' of ', num2str(length(startPoints))]);

end

%%

for j=1:size(savePatterns,1)

	phaseTemp = squeeze(savePatterns(j,:,:));
  ftImage = fftshift(fft2(fftshift(amp.*exp(1i.*(phaseTemp+grate)))));
  figure(1111),
  subplot(1,2,1)
  imagesc(phase);
  subplot(1,2,2)
  imagesc(abs(ftImage));
  pause(1);
  
	slm.show(phaseTemp);

end
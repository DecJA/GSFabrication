kyr =(0:10:100);% linspace(0,100,30);
kxr = ones(size(kyr)).*100;%kxr; %set constant y-offset
CCam =  1.424e-7;

%setup camera
noPics = 1;
picSize = [0 0 1280 1024];
vid = videoinput('ni', 1, 'img0');
vid.ROIPosition = picSize;

px_sizeY = picSize(4);
px_sizeX = picSize(3);
x1=ceil(px_sizeX/2)-1; x2=ceil(px_sizeX/2)+1;
y1=ceil(px_sizeY/2)-1; y2=ceil(px_sizeY/2)+1;

triggerconfig(vid, 'manual');
vid.FramesPerTrigger = 1;

%%
preview(vid);
    
addpath('../')
% for ii=1:length(kxr)
   for jj=1:length(kyr) %change 2 kyr locations for grid
    
    phasePattern = gs_spot_array(kxr(1),kyr(jj),abs(beam));
    slm.show(phasePattern+zeroCorrection);
    pause(0.2);
    %Take image and save
    cd('GSimageCalibration');
    for kk=1:noPics
      image = getsnapshot(vid);
      pause(0.2);
      imName = strcat('X_',num2str(kxr(1)),'_Y_',num2str(kyr(jj)),'_Pic_',num2str(kk),'.bmp');
      imwrite(image,imName,'bmp');
      
      figure(48561),
      imagesc(image)
    end
    cd('../');
    
   end
% end
closepreview(vid)
%%
%Scan through images
allIms = dir('*.bmp');
allImsName = struct2cell(allIms);

%zero order location
x0 = 737;
y0 = 571;

xDistScale = zeros([1,length(kxr)]);
for nn=1:length(kxr)

  singleIm = allImsName(1,nn);
  Im = imread(char(singleIm));
  
  [mr,mc] = find(Im == max(max(Im)));
  if length(mr)>1
    mr = mr(floor(length(mr)/2));
    mc = mc(floor(length(mc)/2));
  end
  dist_to_zero = sqrt((mr-x0)^2 + (mc-y0)^2);
  dx(nn)=(mr);
  dy(nn)=mc;
  
  imSub = Im(mr-20:mr+20,mc-20:mc+20);
  figure(813),
  imagesc(Im); hold on;
  plot(mc,mr,'x','color','r','markersize',14);
  pause(0.1);
%   dx_to_zero = sqrt(mr-x0)^2;
%   
%   xDistScale(nn) = dx_to_zero;
%   xDistM = CCam*xDistScale(nn);
%   
end
figure(),
plot(dx);
% figure(),
% plot(kxr,xDistScale);
%   
% figure(),
% plot(kxr,xDistM,'-+','linewidth',2.5,'markersize',13); grid on; grid minor;
% xlabel('Lateral Scale Value'); ylabel('Displacement from 0th Order');
%Locate zero order
%Locate diffracted spot
%Calculate x/y distance to spot from zero order with CCAM
%Find no pixels between zero and 1st order

%Get converstion
%

%%
for nn=1:length(kyr)

  singleIm = allImsName(1,nn);
  Im = imread(char(singleIm));
  
  figure(813),
  imagesc(Im); hold on;
  
  max_intensity=ginput;
  row(nn)=max_intensity(1);
  col(nn)=max_intensity(2);
  
end
col=sort(col);
pxDistCol = (col-min(col)).*CCam;

row = sort(row);
pyDistRow = (row-min(row)).*CCam;


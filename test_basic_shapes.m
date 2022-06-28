%% GS - Many spots (see other: random holograms, runscript (1 spot, stage + slm) 
%Print using many individual voxels - use minimum distnaces found from
%uniformity programs to reject pattern arrays where spots are too close to
%gether 
rng('default');
slm_check = exist('sz','var');
if slm_check==0  
    addpath('R:\declan\TableD');
    addpath('R:\declan\PrintingSummaryProgram');
    run('setup_slm.m');
    slm_check = 1;
    offset = -(3/2)*pi;
    load('ampSmoothAveraged_objectivePlane_15Dim.mat');
    load('zeroSmooth_ObjectivePlane.mat');
    beam=amp;
    zeroCorrection=zeroSmooth;
end

%% visualise pattern to compare

cutoff = 1.0;
x=linspace(-1,1,sz(1));
[X,Y] = meshgrid(x,x);
R = sqrt(X.^2 + Y.^2);
disk = R<cutoff;
beamCut = disk.*amp;
o = 150;              % Region of interest size in output
padding = 500;        % Padding for FFT
zoom = @(im) im(round(size(im, 1)/2)+(-o:o), round(size(im, 2)/2)+(-o:o));
visualize = @(pattern,amp) zoom(abs(otslm.tools.visualise(pattern, ...
'method', 'fft', 'padding', padding, 'incident', beamCut)).^2);

%% Setup slm and load corrections/amplitudes

%load 2D pattern
densityScale = 2;
coords=[0,0];
xlength=5e-6;
ylength=5e-6;
voxSpacing = 200e-9;
dim  = densityScale*ceil(xlength/voxSpacing);% - mod(dim,10);

xgrid=linspace(0,xlength,dim);
ygrid=xgrid;
[X,Y]=meshgrid(xgrid,ygrid);

disp(['Voxel Spacing is: ',num2str(xlength/dim)]);
shape = generate_triangle([-10e-6,-10e-6],xlength,ylength,dim);
% [shape,~] = generate_rectangle([0,0],xlength,ylength,dim);
shape_bitmap_3d = create_3d_structure(shape,1e-6,200e-9);
shape_bitmap=squeeze(shape_bitmap_3d(:,:,1)).';

% structure = 'image';%Botton left grid cornder locdddaation
% im_name = 'squareTrap2';
% [shape_bitmap,shape_polygon] = shape_map_2d_parse(dim,'coords',coords,...
%     'xdim',xlength,'ydim',ylength,'shape',structure,'image',im_name);
% height_struct = 5e-6;
% res_struct = 500e-9;
% [shape_bitmap_3d] = create_3d_structure(shape_bitmap,height_struct,res_struct);
% [~] = generate_3d_figure(shape_bitmap_3d);


%% 
noTraps = 30;
disp(['Using: ',num2str(noTraps),' traps.']);
exposedPoints = length(find(squeeze(shape_bitmap(:,:))==1));
randIdx = randperm(exposedPoints); 
bitmap_plane = shape_bitmap(:,:);
lin_idx = find(bitmap_plane>0);
[allRows,allCols]=ind2sub([dim,dim],lin_idx);


counter=1;
totalLight = zeros(201);
totalField=zeros(size(visualize(zeros(sz))));
noPatterns = ceil(exposedPoints/noTraps);
ditherPattern=zeros(sz);
savePatterns=zeros([sz(1),sz(2),ceil(dim/densityScale)]);
%% Random point method
for j=1:(ceil(dim/densityScale))
  
    if ((counter+noTraps)>exposedPoints)
      nPrev = noTraps;
      noTraps = noTraps  - ((counter+noTraps)-exposedPoints);
      disp(['No Traps Changed: ',num2str(noTraps)]);
      %Keep power in each trap constant
      powerReduction =  noTraps/nPrev;
      ditherPattern=modulateAmplitude(powerReduction,256,sz);
    end
    
    randLocs = lin_idx(randIdx(counter:counter+noTraps));
    [ir,ic] = ind2sub([dim,dim],randLocs);

    
    for i=1:noTraps
      
      xScale(i) = return_scale_wgs(xgrid(ir(i)));
      yScale(i) = return_scale_wgs(ygrid(ic(i)));
      
    end
    
    %Feed coordinates into GS program - visualise
    phase =  gs_spot_array(xScale,yScale,beamCut,sz);
    totalField = totalField+visualize(phase,beamCut);
    if j==noPatterns
      continue
    else
      figure(4556),
      imagesc(totalField./max(max(totalField))); colormap jet;
      title(['Pattern ',num2str(j),' of ',num2str(ceil(xlength/voxSpacing))])
      hold on;
    end

    savePatterns(:,:,j) = phase+zeroSmooth+ditherPattern;
    
    counter=counter+noTraps;
end


%% Implment Halton code sequence
rng('default')
noTraps=30;
counter=1;
totalLight = zeros(201);
totalField=zeros(size(visualize(zeros(sz))));
noPatterns = ceil(exposedPoints/noTraps);
ditherPattern=zeros(sz);
savePatterns=zeros([sz(1),sz(2),noPatterns]);

disp(['Using: ',num2str(noTraps),' traps.']);
storedValues = [];
avaliablePoints = lin_idx;
[allRows,allCols] = ind2sub([dim,dim],lin_idx);
P=[xgrid(allRows);ygrid(allCols)].';

% k=convhull(xgrid(allRows).',ygrid(allCols).');
k=boundary(xgrid(allRows).',ygrid(allCols).');
polygon=polyshape(P(k,1),P(k,2));

p=haltonset(2);
p=scramble(p,'RR2');
haltonPoints = net(p,1e4).*xlength;
haltonCount = 2;
trapLocs = zeros([1 noTraps]);
trapCount=1; 

savePatternsHalton=zeros([sz(1),sz(2),noPatterns]);
totalIntensity=zeros(size(visualize(zeros(sz))));
figure(123),
plot(polygon); hold on; grid on; grid minor; axis tight;
plot(xgrid(allRows),ygrid(allCols),'.'); hold on;

for n=1:(ceil(dim/densityScale))
    
    if length(allRows)>=1 && length(allRows)<noTraps
        nPrev=noTraps;
        noTraps = noTraps - length(allRows);
        
        powerReduction =  noTraps/nPrev;
        ditherPattern=modulateAmplitude(powerReduction,256,sz);
    end
    
    trapLocs = zeros([1 noTraps]);
    trapCount=1; 
    
    clear P; %clear xSet; clear ySet; 
    while trapCount<(noTraps+1)
    
        %Test halton point interior, test dist from other points
        pTest=[haltonPoints(haltonCount,1),haltonPoints(haltonCount,2)];
        haltonCount=haltonCount+1;
        
        if isinterior(polygon,pTest)==1
            
            %todo: test distance from other spots

            %Find nearest discrete point to go to
            P=[xgrid(allRows);ygrid(allCols)].';
            idxClose = dsearchn(P,[pTest(1); pTest(2)].');
            closePoint = [P(idxClose,1); P(idxClose,2)];
            %Remvoe selected point from avaliable set
            [rowFind,colFind] = find(xgrid(allRows)== closePoint(1) & ygrid(allCols)==closePoint(2));
            
            if length(rowFind)>1
                rowFind(2) = [];
                disp('Picked one of two rows.');
            elseif length(colFind)>1
                colFind(2) = [];
                 disp('Picked one of two columns.');
            end
            allRows(rowFind(1)) = [];
            allCols(colFind(1)) = [];            
            
            xSet(trapCount) = closePoint(1);
            ySet(trapCount) = closePoint(2);
            
            trapCount=trapCount+1;
            clear closePoint; clear idxClose; 
   
        end 
   
    end
     scalex = return_scale_wgs(xSet);
     scaley = return_scale_wgs(ySet);
    [phi_noSpots] = gs_spot_array(scalex,scaley,beamCut,sz); %check if need to add phase
    totalIntensity=totalIntensity + visualize(phi_noSpots,beamCut);

    savePatternsHalton(:,:,n) = phi_noSpots+zeroSmooth+ditherPattern;
%     figure(111),
%     subplot(1,2,1)
%     imagesc(phi_noSpots); colormap jet;
%     subplot(1,2,2)
%     imagesc(totalIntensity); colormap jet;
    if n==noPatterns
      continue
     else
       figure(117),
       imagesc(totalIntensity./max(max(totalIntensity))); colormap jet;
       title(['Pattern ',num2str(n),' of ',num2str((ceil(dim/densityScale)))]) 
    end

    clear scalex; clear scaley;

end

%% Take images and overlap to see Experimental intensity
%rand, halton, w and w/out correction
devices = imaqfind;
lumenera = devices(1);
cam = get(lumenera);

vid = videoinput('winvideo', 1,'RGB16_1280x1024');
set(vid,'FramesPerTrigger',1);
src = getselectedsource(vid);
src.Exposure = -11;

picSize = [833 408 303 282];
vid.ROIPosition = picSize;
runningTotal=zeros(picSize(3),picSize(4));
preview(vid)
%Setup camera

%% Load images onto slm - random distibution
cd('R:\declan\PrintingSummaryProgram\GSPrinting\basicShapeImages\overDensityRandom');
% lin_shift=otslm.simple.linear(sz,20);
% lin_shift=otslm.tools.finalize(lin_shift);
pause(2.0);
for m=1:(ceil(dim/densityScale))
    temp1 = savePatterns(:,:,m);
    slm.show(temp1+lin_shift)
    pause(0.5);
    
   
   pic = getsnapshot(vid);
   pause(0.5)
   name = strcat(num2str(m),'_rand_',num2str(noTraps),'.jpg');
   imwrite(pic,name);
end
    
%% Load images onto slm - halton distibution
cd('R:\declan\PrintingSummaryProgram\GSPrinting\basicShapeImages\overDensityHalton');
pause(2.0);
for m=2:(ceil(dim/densityScale))
    temp1 = savePatternsHalton(:,:,m);
    slm.show(temp1+lin_shift)
    pause(0.5);
   
   pic = getsnapshot(vid);
    pause(0.5);
   name = strcat(num2str(m-1),'_halton_',num2str(noTraps),'_120nm','.jpg');
   imwrite(pic,name);
   pause(1.0);
end

%% Load images and get total intensity reading
seq1 = 'R:\declan\PrintingSummaryProgram\GSPrinting\basicShapeImages\overDensityRandom';
seq2 ='R:\declan\PrintingSummaryProgram\GSPrinting\basicShapeImages\overDensityHalton';

runningTotal=zeros(picSize(3),picSize(4));
cd(seq1);
for k=1:(ceil(dim/densityScale))
  
  files = dir(fullfile(seq1, '*.jpg'));
  fnames = {files.name};
  fnames=natsort(fnames);

    imname=fnames{k};
    im=imread(imname);
    im=rgb2gray(im).';
    
    runningTotal=runningTotal+im2double(im);
    
    figure(3538),
    imagesc(runningTotal./max(max(runningTotal))); colormap jet;
    colorbar
    title('Sequence 1: Random Distribution');
    pause(0.1);
end
clear fnames

    
%%
xres = linspace(-1,1, sz(1));
freq=128;
xrescos=cos(xres.*(freq));
yrescos=cos((xres.').*freq);
xy=1*(pi/2).*sign(xrescos.*yrescos); %1 for maximum diffraction from 0 order
pattern_scatter=xy;

slm.show(pattern_scatter);

slm.image_handle.Parent.CLim=[-pi,pi];
Cmap=slm.lookupTable.value;
Cmap=double(Cmap)./255;
slm.figure_handle.Colormap = Cmap;
slm.image_handle.CDataMapping = 'scaled';
slm.image_handle.CData = pattern_scatter;

%% Load patterns rapidly onto slm
lin_shift=otslm.simple.linear(sz,15,'angle',offset+pi);
lin_shift=otslm.tools.finalize(lin_shift);

dispTime=0.05;
addpath('R:\declan\PrintingSummaryProgram\Runfiles');
currentVoltage=inputSingleScan(d);
sendVz=dist_to_volt(300e-9);
% addAnalogOutputChannel(d,'Dev2',3,'Voltage'); %shutter

vsteps=1;
dPz=(linspace(0,.2.^1,vsteps));
RPB=0;
%%
pause();
for i=1:vsteps
%   [~,ditheredPattern] = modulateAmplitude(RPB+dPz(i),freq,sz);
%   outputSingleScan(d,[i*sendVz,1]);
  for k=1:noPatterns-1
    
    tic
    tempPattern=savePatternsHalton(:,:,k);%+ditheredPattern;
    %tempPattern=savePatterns(:,:,k);
    slm.show(tempPattern + lin_shift);
    b=toc;
    pause(dispTime-b);
    slm.show(tempPattern + lin_shift);
    k
  end
  if vsteps>1
    outputSingleScan(d,[i*sendVz,0]);
    disp(['Moving to step: ',num2str(i),' of ', num2str(vsteps)]);
    pause(1.5);
  end
end

% outputSingleScan(d,[0,0]);

%%
pause();

for k=1:noPatterns-1
  tic

    tempPattern=savePatternsHalton(:,:,k);
%      tempPattern=savePatterns(:,:,k);
   slm.show(tempPattern + lin_shift);

  b=toc;
  pause(dispTime-b);
  k
end

outputSingleScan(d,[0,0]);

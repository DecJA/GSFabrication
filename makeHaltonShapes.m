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
o = 200;              % Region of interest size in output
padding = 500;        % Padding for FFT
zoom = @(im) im(round(size(im, 1)/2)+(-o:o), round(size(im, 2)/2)+(-o:o));
visualize = @(pattern,amp) zoom(abs(otslm.tools.visualise(pattern, ...
'method', 'fft', 'padding', padding, 'incident', beamCut)).^2);

%% Setup slm and load corrections/amplitudes

%load 2D pattern
densityScale = 2;
coords=[0,0];
xlength=10e-6;
ylength=10e-6;
voxSpacing = 400e-9;
dim  = densityScale*ceil(xlength/voxSpacing);% - mod(dim,10);

xgrid=linspace(0,xlength,dim);
ygrid=xgrid;
xc = (xlength/2) + coords(1);
yc = (ylength/2) + coords(2);
[X,Y]=meshgrid(xgrid,ygrid);

disp(['Voxel Spacing is: ',num2str(xlength/dim)]);
% shape = generate_triangle([-10e-6,-10e-6],xlength,ylength,dim);
% [shape,~] = generate_rectangle([0,0],xlength,ylength,dim);
% shape_bitmap_3d = create_3d_structure(shape,1e-6,200e-9);
% shape_bitmap=squeeze(shape_bitmap_3d(:,:,1)).';

addpath('GSIllustratorPatterns')
structure = 'image';%Botton left grid cornder locdddaation
im_name = 'circleTrap';
coords = [0,0];
[shape_bitmap,shape_polygon] = shape_map_2d_parse(dim,'coords',coords,...
    'xdim',xlength,'ydim',ylength,'shape',structure,'image_name',im_name);
height_struct = 5e-6;
res_struct = 500e-9;
[shape_bitmap_3d] = create_3d_structure(shape_bitmap,height_struct,res_struct);
[~] = generate_3d_figure(shape_bitmap_3d);


%% 
noTraps = 10;
disp(['Using: ',num2str(noTraps),' traps.']);
exposedPoints = length(find(squeeze(shape_bitmap(:,:))==1));
randIdx = randperm(exposedPoints); 
bitmap_plane = shape_bitmap(:,:);
lin_idx = find(bitmap_plane>0);

%% Implment Halton code sequence
rng('default')
noTraps=10;
counter=1;
patternCount = 0;
totalLight = zeros(201);
totalField=zeros(size(visualize(zeros(sz))));
noPatterns = ceil(exposedPoints/noTraps);
savePatterns=zeros([sz(1),sz(2),noPatterns]);

disp(['Using: ',num2str(noTraps),' traps.']);
storedValues = [];
avaliablePoints = lin_idx;
[allRows,allCols] = ind2sub([dim,dim],lin_idx);
P=[xgrid(allRows); ygrid(allCols)].';

% k=convhull(xgrid(allRows).',ygrid(allCols).');
k=boundary(xgrid(allRows).',ygrid(allCols).');
polygon=polyshape(P(k,1),P(k,2));

p=haltonset(2);
p=scramble(p,'RR2');
haltonPoints = net(p,1e4).*xlength;
haltonCount = 1;
trapLocs = zeros([1 noTraps]);
trapCount=1; 

savePatternsHalton=zeros([sz(1),sz(2),noPatterns]);
totalIntensity=zeros(size(visualize(zeros(sz))));

figure(123),
plot(polygon); hold on; grid on; grid minor; axis tight;
plot(xgrid(allRows),ygrid(allCols),'.');

ditherPattern = zeros(sz);

multiBoundary = 1;
if multiBoundary ==1
    [scalePoly] = embed_shape(polygon,xgrid(allRows),ygrid(allCols));
end


%% Implment Halton method
for n=1:floor(length(lin_idx)/noTraps)
    
    %Modulate power if trap number is reduced
    if length(allRows)>=1 && length(allRows)<noTraps
        nPrev=noTraps;
        noTraps = noTraps - length(allRows);
        
        powerReduction =  noTraps/nPrev;
        ditherPattern=modulateAmplitude(powerReduction,256,sz);
    end
    
    trapLocs = zeros([1 noTraps]);
    trapCount=1; 
    
    clear P; clear xSet; clear ySet; 
    while trapCount<(noTraps+1)
        
        %Test halton point interior, test dist from other points
        pTest=[haltonPoints(haltonCount,1),haltonPoints(haltonCount,2)];
        haltonCount=haltonCount+1;
        
        isInShape = isinterior(polygon,pTest);
        isInShape2 = 0;
        if  multiBoundary == 1
            isInShape2 = isinterior(scalePoly,pTest);
            %Reject point if not in band
        end
        
        if isInShape == 1 && isInShape2 == 0
            
            %Find nearest discrete point to go to
            P=[xgrid(allRows);ygrid(allCols)].';
            idxClose = dsearchn(P,[pTest(1); pTest(2)].');
            closePoint = [P(idxClose,1); P(idxClose,2)];
            %Remvoe selected point from avaliable set
            [rowFind,colFind] = find(xgrid(allRows)== closePoint(1) & ygrid(allCols)==closePoint(2));
            
            allRows(rowFind(1)) = [];
            allCols(colFind(1)) = [];            
            
            if length(allRows)<noTraps
              disp('No Points Left')
              return
            end
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

    if n==noPatterns
      continue
     else
       figure(117),
       imagesc(totalIntensity./max(max(totalIntensity))); colormap jet;
       title(['Pattern ',num2str(n),' of ',num2str(length(lin_idx)/noTraps)]) 
    end

    clear scalex; clear scaley;
    patternCount = patternCount+1;

end

%% Rapidly load onto slm
if exist('d')==0
  d = daq.createSession('ni');
  send=1;
  addAnalogInputChannel(d,'Dev2',3,'Voltage'); %shutter
  addAnalogOutputChannel(d,'Dev2',3,'Voltage'); %shutter

  outputSingleScan(d,[send]);
  pause(0.05)
  outputSingleScan(d,[0]);
end

%% 2D Pattern
lin_shift = otslm.simple.linear(sz,8,'angle',offset);
lin_shift = otslm.tools.finalize(lin_shift);

pause();
outputSingleScan(d,[send]);
dispTime = 0.05;

for k=1:patternCount-30
  tic
  tempPattern=savePatternsHalton(:,:,k);
  slm.show(tempPattern + lin_shift);
  b=toc;
  pause(dispTime-b);
  k
end

outputSingleScan(d,[0]);

%% Axial Shift
% clear d chi idi ch id
if exist('d') == 0
   d = daq.createSession('ni');
  [chi,idi] = addAnalogInputChannel(d,'Dev2',[2,3],'Voltage'); %2 - shutter, 3-axial
  [ch,id] = addAnalogOutputChannel(d,'Dev2',[2,3],'Voltage'); %Add only z-axis
end

axialHeight = 5e-6;
axialResolution = 500e-9;
axialVoltage = dist_to_volt(axialResolution);
axialPoints = floor(axialHeight/axialResolution);

lin_shift = otslm.simple.linear(sz,6,'angle',offset);
lin_shift = otslm.tools.finalize(lin_shift);

dispTime = 0.05;

%% - continued axial shift
pause();

for h=1:axialPoints
  %Open Shutter - move stage
  sendVoltage = dist_to_volt((h-1)*axialResolution);
  outputSingleScan(d,[sendVoltage,1]);
  for k=1:patternCount
    tic
    tempPattern=savePatternsHalton(:,:,k);
    slm.show(tempPattern + lin_shift);
    b=toc;
    pause(dispTime-b);
    k
  end
end
  
outputSingleScan(d,[0,0]);


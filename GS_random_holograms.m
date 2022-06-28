addpath('R:\declan\PrintingSummaryProgram\Runfiles');
% run('setup_slm.m'); %Use beamCorrection for 510
addpath('../')
load('ampSmoothAveraged_objectivePlane_15Dim.mat')
load('zeroSmooth_ObjectivePlane.mat')

offset = -(3/2)*pi;
zeroCorrection=zeroSmooth;
Ebeam = abs(amp);
weight = 50;

%%
board_check = exist('board_check','var');
startAxialVoltage = 0;
if board_check==0
  d = daq.createSession('ni');
  %Make sure no voltage is output
  [chi,idi] = addAnalogInputChannel(d,'Dev2',2,'Voltage'); %Only z-axis
  [ch,id] = addAnalogOutputChannel(d,'Dev2',2,'Voltage'); %Add only z-axis
  outputSingleScan(d,startAxialVoltage); %Start with large positive voltage for axial dim.
  initialVoltage = inputSingleScan(d);
  board_check = 1;  
end

%% Load image
sz=[510,510]; offset=-(3/2)*pi;
addpath('R:\declan\TableD\GStestPatterns');
addpath('R:\declan\PrintingSummaryProgram\GSPrinting\GSIllustratorPatterns');
addpath('R:\declan\TableD\2PP_Testing\IllustratorImages');
% Ebeam=ampSmooth;
imScale = zeros(sz);
im=imread('Q2.PNG');
im=im(:,:,1);
im=im2bw(im);
im=abs(im-1);
im=imresize(im,0.32); %Use for scaing image/pattern size

im2 = padarray(im,[floor((sz(1)-size(im(:,1),1))/2)+1,floor((sz(2)-size(im(1,:),2))/2)+1],1,'both');
im2=im2(1:sz(1),1:sz(2));
im2=abs(1-im2);
lens = otslm.simple.spherical(sz,510).*(2*pi)*200;
Edoe = (Ebeam).*exp(1i.*lens);
targetAmp=im2;

figure(),
imagesc(targetAmp);
disp(['Size: ',num2str(size(targetAmp))]);

%Smooth out image
grating = otslm.simple.linear(sz,15,'angle',offset+pi/4);
grating = otslm.tools.finalize(grating);
%% Create phase patterns'
weight=50;
cutoff = 1.0;
x=linspace(-1,1,sz(1));
[X,Y] = meshgrid(x,x);
R = sqrt(X.^2 + Y.^2);
disk = R<cutoff;

noPatterns=60;
Esave = zeros([noPatterns,sz(1),sz(2)]);

for k=1:noPatterns
    phiRand = randn([sz(1),sz(2)]);
    Edoe = Ebeam + weight*exp(1i.*phiRand);
    for j=1:50
        Efocus = fftshift(fft2(fftshift(Edoe)));
        phifocus = angle(Efocus);
        targetE = abs(targetAmp).*exp(-1i.*phifocus);
        Einv = fftshift(ifft2(fftshift(targetE)));
        phihot = angle(Einv);
%         Edoe = disk.*(abs(Ebeam).*exp(-1i.*phihot));
        Edoe = ((abs(Ebeam)).*exp(-1i.*phihot));
    end
    k
    Esave(k,:,:) = phihot;
end

disp('done')

%% Display on SLM
if exist('d')==0
  d = daq.createSession('ni');
  send=1;
  addAnalogInputChannel(d,'Dev2',3,'Voltage'); %shutter
  addAnalogOutputChannel(d,'Dev2',3,'Voltage'); %shutter
end
pause();
for m=1:noPatterns-20
    epat = squeeze(Esave(m,:,:)).';
  
    slm.show((epat+grating+zeroSmooth));
    outputSingleScan(d,[send]);
    pause(0.05);
    outputSingleScan(d,[0]);
end

disp('done')


%% %% Scan in z

dispTime=0.05;

addpath('R:\declan\PrintingSummaryProgram\Runfiles');
currentVoltage=inputSingleScan(d);
sendVz=dist_to_volt(300e-9);
vsteps=20;

pause();
for i=1:vsteps
  for k=1:noPatterns-1
    tic
    epat = squeeze(Esave(k,:,:)).';
    slm.show((epat+grating+zeroSmooth.*0))
    b=toc;
    pause(dispTime-b);
  end
  disp(['Moving to step: ',num2str(i),' of ', num2str(vsteps)]);
  if vsteps>1
    outputSingleScan(d,i*sendVz);
  end
end

% outputSingleScan(d,0);
traps=15;
rows=1;
squareDist=25e-6;
spotY = ones([rows, traps]);
spotX = ones([rows traps]);

spotX(1,:) = linspace(-squareDist/2,squareDist/2,traps);
spotX(2,:) = linspace(-40,40,traps);
spotY(1,:) = -3e-6;
spotY(2,:) = 10e-6;

kxr = return_scale_wgs(spotY(1:rows,:));
kyr = return_scale_wgs(spotX(1:rows,:));

testPat = gs_spot_array(kxr,kyr,abs(ampSmooth),sz);

slm.show(testPat+zeroSmooth);

%%
addpath('R:\declan\PrintingSummaryProgram\Runfiles');
pause();

dx=200e-9;
totalDist = 200e-9;
for ii=1:100
  
  totalDist = totalDist + dx;
  grateScale = dist_to_scale(totalDist);
  
  grating=otslm.simple.linear([sz(1),sz(2)],grateScale,'angle_deg',180);
  grating=otslm.tools.finalize(grating);

  
  slm.show(testPat + grating+zeroSmooth.*0)
  pause(0.1);
end
  
function [phi_jx2] = gs_spot_array(kxr,kyr,beam,ampWeights,sz)


%   figure(1123),
%   plot(kxr(:),kyr(:),'x','markersize',10,'linewidth',2,'color','k'); 
%   grid on;
  
  [X,Y]=meshgrid(pi*1*linspace(-1,1,sz(1)),pi*linspace(-1,1,sz(2)));
  szx=size(X);

  % turn matrix into column vectors
  X=X(:);
  Y=Y(:);
  amp=ampWeights;%ones(size(kxr));%.*cos(linspace(0,pi,length(kxr))); %amplitude of each spot
  kx1=kxr(:)./1;
  ky1=kyr(:)./1;
  amp=amp(:);

  %%
  kx=kx1(:);
  ky=ky1(:);

  N=length(kx);
  indx=find(sqrt(X.^2+Y.^2)<pi);
  FN=exp(1i*(kx(:,end).*X'+ky(:,end).*Y')); %%all relevant modes for "fourier series"

  % % initialise V
  phi_jx=(-1).^(round([1:size(X,1)]/1))*pi/2;
%   phi_jx=zeros([1 512]);

  % create the spot intensity array:
  Ax=abs(beam(:)).';
  Vx=((Ax.*exp(1i*phi_jx))*(FN)').'/length(X);

  [phi_jx,Sx]=WGSalgorithm(Vx,amp,FN,Ax);

  phi_jx2 = reshape(phi_jx,[sz(1),sz(2)]);


%   figure(114),imagesc(phi_jx2)

end
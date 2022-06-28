%% initialisation and parameters for loop. algorithm in: 
function [phi_j,SLMout]=WGSalgorithm(V,amp,FN,A0)


wm=ones(size(V));
wfun=sqrt(abs(amp));

imax=50;
jmax=70;
vt=ones(size(wm));

for ii=1:imax
    
    ifunt=min(1,min(ii,jmax+0)/jmax);
    ifun=ifunt;
    
    % compute weights to apply to modes
    wm=wm.*(min(abs(V))./abs(V)).^ifun.*wfun;
    wm=wm/vecnorm(wm);
    
    % generate the current SLM phase pattern to try:
    SLMout = (FN.'*(wm.*V./abs(V)));
    phi_j=(angle(SLMout)');
    
    % compute the amplitude at the SLM.
    I_j=abs((FN.'*(wm.*V./abs(V))))';
    I_j=I_j/max(I_j);
    
    % calculate the complex amplitude at the spots.
    V=((A0.*exp(1i*phi_j))*(FN')).'/length(A0);
    aV=abs(V).^2;
    vt=aV/max(aV);
    
end
vt(end); %print out relative spot brightness.

%%
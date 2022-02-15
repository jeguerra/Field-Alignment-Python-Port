function [qtx,qty] = FA2DC4V2_simple(Xb,Y,niter,wt,lscale,mode)
%% Description:
% Field Alignment 2D version, Code branch 3, version 3
% options.
% Input Xb1
% Input Y1
% H1 -- locations where observed
% sq = size of problem domain, square
% niter = total available iterations.
% graphics -- should we display something
% wt -- what's the constraint weight.
% mode -- power law mode 2 or 4.
% debug --
%% Copyright  © Sai Ravela, MIT, 2003-2012.
% This code is part of the Field Alignment System and Testbed,
% developed by Ravela et al. at the Massachusetts Institute of Technology.
% All rights reserved.

convTol = 1.0E-6;
szx = size(Xb); 
nx = szx(1);
ny = szx(2);
pblmr1 = 1:nx;
pblmr2 = 1:ny;

% Compute the nominal pixel grid
[aa,bb] = meshgrid(pblmr2,pblmr1);

% Normalize Xb and Y -- not always necessary.
mxb = min(Xb(:)); 
dxb = max(Xb(:)) - mxb;
Xb = (Xb - mxb)./dxb;
Y = (Y - mxb)./dxb;

%% We will solve spectrally.
lscale=round((lscale*nx+nx)/2)*2;
[m,n] = meshgrid([-lscale/2:-1 1e-8 1:lscale/2-1], [-lscale/2:-1 1e-8 1:lscale/2-1]);

% Initialize forcing
ux = zeros(nx,ny); 
uy = ux; 
qtx = ux; 
qty = ux; 

w1 = 1; w2 = w1/3;% A "Newtonian Fluid."
iter = 0; 
errv0 = inf; 
derrv = inf;

while((iter<niter) && derrv > convTol)
    
    % Index displacements to region index (FFT solver changes dimensions)
    ssx = ux(pblmr1,pblmr2); 
    ssy = uy(pblmr1,pblmr2);
    
    %% p - q, interpolate.
    qx = aa - ssx;        
    qy = bb - ssy;
    qtx = qtx + ssx; 
    qty = qty + ssy; %% You can also interpolate with total disp.
    
    % Apply displacement with interpolation (suspect...)
    [bdx,bdy] = meshgrid(1:size(Xb,2),1:size(Xb,1));
    XXb1 = interp2(bdx,bdy,Xb,qx,qy,'makima');
    %XXb1 = bicutest(bdx,bdy,XXb1,bd+qx,bd+qy);
    Xb = XXb1;        
    
    %% Calc forcing here
    dq = (Xb - Y);
    dq(1,1:end) = 0; 
    dq(end,1:end) = 0; 
    dq (1:end,end) = 0; 
    dq(1:end,1) = 0;
    
    % Initiate force vectors and compute gradients (upgrade!)
    t1 = zeros(size(Xb)); 
    t2 = t1;
    [dXbx,dXby] = gradient(Xb(:,:));
    t1 = t1 - dXbx.*dq;            
    t2 = t2 - dXby.*dq;
    
    %% Solve system at current iteration. Spectral solution.
    % Mode = even, 2 or 4 is all that is needed in most problems.
    % Yang and Ravela SCA method 2009 (Yang's Masters Thesis)
    
    % Moving over to Fourier domain
    fFx = (fftshift(fft2(t1,lscale,lscale)));    
    fFy = (fftshift(fft2(t2,lscale,lscale)));
    % Avoid singularity
    fFx(lscale/2+1,lscale/2+1) = 0;     
    fFy(lscale/2+1,lscale/2+1) = 0;
    
    %Apply deformation filters
    fac0 = w1.*(n.^mode + m.^mode);
    fac1 = w2.*n.^2 + fac0;
    fac2 = w2.*m.*n;
    fac4=  w2.*m.^2 + fac0;
    
    xpy = (fac1 - fac2); 
    ypx = (fac4 - fac2);  
    zz = (fac1.*fac4 - fac2.*fac2) * wt/lscale;
    
    fux = (-fFx.*xpy + fac2.*fFy)./zz; ux = real(ifft2(ifftshift(fux)));
    fuy = (fFx.*fac2 - ypx.*fFy)./zz;  uy = real(ifft2(ifftshift(fuy)));
   
    % Change in the max norm of displacements signals minimum
    errv1 = max(max(max(abs(ux))), max(max(abs(uy))));
    derrv = abs(errv1 - errv0);
    errv0 = errv1;
    
    iter = iter+1;
end

sderrv = sprintf('%0.5e',derrv);
siter = sprintf('%i',iter);
disp([sderrv, ' ', siter]);
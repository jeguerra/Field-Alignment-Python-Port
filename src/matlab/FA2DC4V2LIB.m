function [qtx,qty,qq,errl] = FA2DC4V2LIB(Xb1,Y1,H1,Xu1,Yu1,Hu1,msk,niter,wt,lscale,mode)
%% Description:
% Field Alignment 2D version, Code branch 3, version 3
% options.
% Input Xb1, Xb2, Xb3 = three fields  normalized.
% Input Y1, Y2, Y2 = three fields  normalized.
% H1 -- locations where observed
% sq = size of problem domain, square
% niter = total available iterations.
% graphics -- should we display something
% wt -- what's the constraint weight.
% mode -- power law mode 2 or 4.
% debug --
%% Copyright  � Sai Ravela, MIT, 2003-2012.
% This code is part of the Field Alignment System and Testbed,
% developed by Ravela et al. at the Massachusetts Institute of Technology.
% All rights reserved.

szx = size(Xb1); 
nx = max(szx); %%% We just assume it is square set to max dimension
pblmr = 1:nx;
Xb  = zeros(nx); 
H = zeros(szx);
Hu = zeros(nx);
mxb = zeros(size(Xb1,3));
dxb = mxb;
for i = 1:size(Xb1,3)
    Xb = Xb1(:,:,i);
    Y = Y1(:,:,i);
    H(:,:,i) = H1(:,:,i); 
    H((Xb)==inf) = 0;H((Y)==inf) = 0;
    H((Xb)==-999999) = 0;H((Y)==-999999) = 0;
    Xb((Xb)==inf) = 0;Y((Y)==inf) = 0;
    Xb((Xb)==-999999) = 0; Y((Y)==-999999) = 0;
    
    % Normalize Xb and Y -- not always necessary.
    mxb(i) = min(Xb(:)); 
    dxb(i) = max(Xb(:)) - mxb(i);
    Xb = (Xb - mxb(i))./dxb(i);
    Y = (Y - mxb(i))./dxb(i);
    Xb1(:,:,i) = Xb;
    Y1(:,:,i) = Y;
end

if (~isempty(Hu1))
    Hu = Hu1;
    for i = 1:size(Xu1,3)
        Xu = squeeze(Xu1(:,:,i));
        Yu = squeeze(Yu1(:,:,i));
        Hu((Xu)==inf) = 0;Hu((Yu)==inf) = 0;
        Hu((Xu)==-999999) = 0;Hu((Yu)==-999999) = 0;
        Xu((Xu)==inf) = 0;Yu((Yu)==inf) = 0;
        Xu((Xu)==-999999) = 0; Yu((Yu)==-999999) = 0;

        % Normalize Xb and Y -- not always necessary.

        Xu1(:,:,i) = Xu;
        Yu1(:,:,i) = Yu;
    end
    mg  = sqrt(Xu1(:,:,1).*Xu1(:,:,1)+Xu1(:,:,2).*Xu1(:,:,2));
    mxmg = max(mg(:));
    Xu1 = Xu1./mxmg; Yu1 = Yu1./mxmg;
end

%% We will solve spectrally.
lscale=round((lscale*nx+nx)/2)*2;
[m,n] = meshgrid([-lscale/2:-1 1e-8 1:lscale/2-1], [-lscale/2:-1 1e-8 1:lscale/2-1]);

ux = zeros(nx,nx); uy = ux; qtx = ux; qty = ux; 
qq = zeros(szx); %qu = zeros(szu);

for i = 1:size(Xb1,3)
        qq(:,:,i) =  Xb1(:,:,i).*dxb(i)+mxb(i); 
end
w1 = 1; w2 = w1/3;% A "Newtonian Fluid."
iter = 0; err =[];
errl = inf; errv0 = inf; derrv = inf;
while((iter<niter) && derrv > 1e-6)% && errv > 1e-6)% && (derrl)<1e+1)
    % Clamp the displacement so it does not go wild.
    ssx = min(ux,nx/2); ssy = min(uy,nx/2);
    
    % Yank problem region displacement, but make consistent with NCEP
    % implementation  on Halo placement...to do 5/15/12 Sai
    ssx1 = msk.*ssx(pblmr,pblmr); 
    ssy1 = msk.*ssy(pblmr,pblmr);
    [aa,bb] = meshgrid(1:length(ssx1),1:length(ssx1));
    
    %% p - q, interpolate.
    qx = aa - ssx1;        
    qy = bb - ssy1;
    qtx = qtx + ssx1; 
    qty = qty + ssy1; %% You can also interpolate with total disp.
    t1 = zeros(size(Xb)); t2 = t1;dqe=0;
    for i = 1:size(Xb1,3)
        XXb1 = Xb1(pblmr,pblmr,i); %% make sure interp does not tank at edge.
        bd = min(25,size(XXb1,1));
        XXb1 = Border(XXb1,bd);
        [bdx,bdy]=meshgrid(1:size(XXb1,2),1:size(XXb1,1));
        XXb1 = interp2(bdx,bdy,XXb1,bd+qx,bd+qy,'makima');
        %XXb1 = bicutest(bdx,bdy,XXb1,bd+qx,bd+qy);

        Xb1(pblmr,pblmr,i) = XXb1;        
        qq(:,:,i) =  Xb1(pblmr,pblmr,i).*dxb(i)+mxb(i); 
        %% Calc forcing here
        dq = H(:,:,i).*(Xb1(:,:,i) - Y1(:,:,i));
        dq(1,1:end) = 0; 
        dq(end,1:end) = 0; 
        dq (1:end,end) = 0; 
        dq(1:end,1) = 0;
        [dXbx,dXby] = gradient(Xb1(:,:,i));
        t1 = t1-dXbx.*dq;            
        t2 = t2-dXby.*dq;
        dqe = median([dqe;abs(dq(H(:,:,i)>0))]);
    end
    
    if (~isempty(Hu1))    
        UUb1 = Xu1(pblmr,pblmr,1);
        bd = min(10,size(UUb1,1));
        UUb1 = Border(UUb1,bd); 
        [bdx,bdy]=meshgrid(1:size(UUb1,2),1:size(UUb1,1));
        Ub1 = bicutest(bdx,bdy,UUb1,bd+qx,bd+qy);

        VVb1 = Xu1(pblmr,pblmr,2);
        bd = min(10,size(VVb1,1));
        VVb1 = Border(VVb1,bd); 
        [bdx,bdy]=meshgrid(1:size(VVb1,2),1:size(VVb1,1));
        Vb1 = bicutest(bdx,bdy,VVb1,bd+qx,bd+qy);

        [qxx,qxy] = gradient(qx); [qyx,qyy] = gradient(qy);

        Vbt = qxx.*Vb1 -qyx .* Ub1;
        Ubt = qyy.*Ub1 -qxy .* Vb1;

        Xu1(:,:,1) = Ubt; 
        Xu1(:,:,2) = Vbt;

        dqu = Hu.*(Xu1(:,:,1) - Yu1(:,:,1));    
        [dUbx,dUby] = gradient(Xu1(:,:,1));  
        tu1 = -1.*dUbx.*dqu;    
        tu2 = -1.*dUby.*dqu;
        dqv = Hu.*(Xu1(:,:,2) - Yu1(:,:,2));    
        [dVbx,dVby] = gradient(Xu1(:,:,2));  
        tv1 = -1.*dVbx.*dqv;    
        tv2 = -1.*dVby.*dqv;
        
        %Calculate total force
        t1 = t1+tu1+tv1; 
        t2 = t2+tu2+tv2;
    end

    %%String along the residual (forcing) here.
    err = [err dqe];
    
    %% Solve system at current iteration. Spectral solution.
    % Mode = even, 2 or 4 is all that is needed in most problems.
    %Note: Tikhonov is no longer necessary, superseeded by Yang and
    %Ravela 2009, SCA method.
    
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
   
    % We expect convergence on ux, instantaneous deformation -> to zero. This
    % of course is an assumption and only valid for appropriately formulated
    % local optimization. Nevertheless, we'll track it and use it as a
    % criterion.
    errv1 = max(max(max(abs(ux))), max(max(abs(uy))));
    derrv = abs(errv1 - errv0);
    errv0 = errv1;
    
    % Logging here.
    if (iter > 2)
        errl = err(length(err));
        derrl = errl - err(length(err) -1);
       
    end
    
    iter = iter+1;
end

sderrv = sprintf('%0.5e',derrv);
siter = sprintf('%i',iter);
disp([sderrv, ' ', siter]);
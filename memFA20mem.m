clc;
close all;
clear;
%pkg load netcdf
t1 = clock;
%Field alignment
%input:
%---two fields---
%A and B
%output:
%--C--> qx, qy : vector field
%--D align B to A
%input files------------------
%dir = ...    filename = ...
tdex = 12;
fname = 'kasatochi_era5';
fiinput = ['test_data/' fname '.nc'];
%----read data---

nmem = length(ncread(fiinput, 'ens'));
lats = ncread(fiinput, 'y');
nlat = length(lats);
lons = ncread(fiinput, 'x');
nlon = length(lons);
hmem0 = ncread(fiinput, '__xarray_dataarray_variable__');
% Arrange matrix in (lat, lon, ensid)
hmem0 = squeeze(hmem0(:,:,tdex,:));
hmemT = permute(hmem0,[2,1,3]);
%============================== 
fprintf('Start Calculation');
%imem=1;    %++++++++++++++++++++where you need to change+++++++
%fprintf('This is the member',imem)
%----parameters for field alignment
iter = 2^10;
lscale = 1;   %change from 1~10,   1 is the best 

%------------------------
%==========Main algorithm=========
hmem_FA = zeros(nlon,nlat,nmem,'double');
for imem = 1:1:nmem  % reference member
	fprintf('reference mem is:');
	disp(imem);
	qxyT_mem = zeros(nlat,nlon,2,nmem,'double');

	%------derive the displacement vectors-------
	%hmemT(:,:,1),hmemT(:,:,2) ->->-> move 1 to 2
	%hmemT(:,:,i),hmemT(:,:,:) ->->-> move i to 1~20 member
	for i = 1:1:nmem   % reference mem to different members
		fprintf('imem is:'); disp(i);
		[qxyT_mem(:,:,1,i),qxyT_mem(:,:,2,i)] = ...
			FA2DImNoHLIB(hmemT(:,:,imem),hmemT(:,:,i),128,.1,iter,2,lscale);
	end

	%------calculate the mean of all DVs----
	qxyT = zeros(nlat,nlon,2,'double');
	qxy = zeros(nlon,nlat,2,'double');
	
    for i = 1:1:nmem
		qxyT(:,:,1) = qxyT(:,:,1) + qxyT_mem(:,:,1,i);
		qxyT(:,:,2) = qxyT(:,:,2) + qxyT_mem(:,:,2,i);
    end
    
	qxyT(:,:,1) = qxyT(:,:,1) ./ nmem;
	qxyT(:,:,2) = qxyT(:,:,2) ./ nmem;
	qxy(:,:,1) = transpose(qxyT(:,:,1));    %output after transopose:
	qxy(:,:,2) = transpose(qxyT(:,:,2));    %output after transopose:
	fprintf('Derive the displacement vector successfully!')

	%-----align the member i----
	%hmemT_FA = zeros(nlat,nlon,'double');
	hmemT_FA = advect(hmemT(:,:,imem),-qxyT(:,:,1),-qxyT(:,:,2));
	hmem_FA(:,:,imem) = transpose(hmemT_FA);    %output after transopose:
	fprintf('align member imem successfully!')

	%-----output data to new NC mean files----
	fileout = ['/test_data/' fname '_FA.nc'];
end
t2 = clock;

%% Plot original member data and aligned member data
[X,Y] = meshgrid(lats, lons);

% Visual window
xwin = [-0.5, +0.5];
ywin = [-77.8, -77.5];

% Plot original data into left axes
ax(1) = subplot(221);
for p = 1:nmem
	cn = contour(ax(1), X, Y, hmem0(:,:,p), 11);
    hold on;
end
hold off;
alpha(0.1)
%xlim(xwin); ylim(ywin);
title('Initial Ensemble');
	
% Plot aligned data into right axes
ax(2) = subplot(222);
for p = 1:nmem
	cn = contour(ax(2), X, Y, hmem_FA(:,:,p), 11);
    hold on;
end
hold off;
alpha(0.1)
%xlim(xwin); ylim(ywin);
title('Aligned Ensemble');

% Plot means
ax(3) = subplot(223);
cn = contour(ax(3), X, Y, mean(hmem0, 3), 11);
%xlim(xwin); ylim(ywin);
colorbar;
title('Arithmetic Mean');
ax(4) = subplot(224);
cn = contour(ax(4), X, Y, mean(hmem_FA, 3), 11);
%xlim(xwin); ylim(ywin);
colorbar;
title('Field-aligned Mean')
	
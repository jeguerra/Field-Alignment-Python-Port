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
zdex = 10;
fname = 'VeniMass';
fiinput = ['test_data/' fname '.nc'];
%----read data---

nmem = 21; % up to 13
lats = (ncread(fiinput, 'latitude'));
nlat = size(lats); nlat = nlat(2);
lons = (ncread(fiinput, 'longitude'));
nlon = size(lons); nlon = nlon(1); 
hmemAll = ncread(fiinput, '__xarray_dataarray_variable__');
osize = size(hmemAll);

raw_mean = zeros(nlon,nlat);
alg_mean = zeros(nlon,nlat);
pdex = [2,1,3];
for tt=1:1
    % Arrange matrix in (lat, lon, ensid)
    hmem0 = squeeze(hmemAll(:,:,end,:));
    hmemT = permute(hmem0,pdex);

    lats = permute(lats,pdex);
    lons = permute(lons,pdex);

    % Apply windowing function to the data
    win1 = tukeywin(nlat,0.25);
    win2 = tukeywin(nlon,0.25);
    win = (win1 * win2');
    K = 0.5*ones(3);
    for imem = 1:1:nmem
        % Take the ln and window
        %hmemT(:,:,imem) = log(hmemT(:,:,imem)) .* win;
        % Apply windowing only
        hmemT_nominal = hmemT(:,:,imem);
        hmemT(:,:,imem) = hmemT_nominal .* win;
        %hmemT_smooth = conv2(hmemT(:,:,imem),K,'same');
        %hmemT(:,:,imem) = hmemT_smooth .* win;
    end

    %============================== 
    fprintf('Start Calculation');
    %imem=1;    %++++++++++++++++++++where you need to change+++++++
    %fprintf('This is the member',imem)
    %----parameters for field alignment
    iter = 2^13;
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
                FA2DImNoHLIB(hmemT(:,:,imem),hmemT(:,:,i),256,1.0E-2,iter,4,lscale);
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
        hmemT_FA = advect(hmemT(:,:,imem),qxyT(:,:,1),qxyT(:,:,2));
        hmem_FA(:,:,imem) = transpose(hmemT_FA);    %output after transopose:
        fprintf('align member imem successfully!')

        %-----output data to new NC mean files----
        fileout = ['/test_data/' fname '_FA.nc'];
    end
    t2 = clock;
    
    raw_mean = raw_mean + mean(hmem0, 3);
    alg_mean = alg_mean + mean(hmem_FA, 3);
end

%% Plot original member data and aligned member data
%hmem_FA = exp(hmem_FA);
%[plons,plats] = meshgrid(lats, lons);

plats = lats';
plons = lons';
% Visual window
%xwin = [min(min(lons)), max(max(lons))];
%ywin = [min(min(lats)), max(max(lats))];
xwin = [215.0, 235.0];
ywin = [40.0, 65.0];

% Set a contour level
cl = 1.0 * [1.0, 1.0];

% Plot original data into left axes
ax(1) = subplot(221);
for p = 1:2:nmem
	[Mc, cn] = contour(ax(1), plons, plats, hmem0(:,:,p), cl);
    cn.LineWidth = 2;
    cn.LineColor = 'k';
    grid on;
    hold on;
end
hold off;
alpha(0.1)
xlim(xwin); ylim(ywin);
ax(1).FontSize = 12;
%xlabel('Deg. West');
ylabel('Deg. North');
title('Initial Ensemble');
	
% Plot aligned data into right axes
ax(2) = subplot(222);
for p = 1:2:nmem
	[Mc, cn] = contour(ax(2), plons, plats, hmem_FA(:,:,p), cl);
    cn.LineWidth = 2;
    cn.LineColor = 'k';
    grid on;
    hold on;
end
hold off;
alpha(0.1)
xlim(xwin); ylim(ywin);
ax(2).FontSize = 12;
title('Aligned Ensemble');

% Plot means
lim1 = 0.0;
lim2 = 20.0;
colormap('jet')
ax(3) = subplot(223);
cn = contourf(ax(3), plons, plats, mean(hmem0, 3), 21);
xlim(xwin); ylim(ywin);
caxis([lim1 lim2])
colorbar;
ax(3).FontSize = 12;
title('Arithmetic Mean');
xlabel('Deg. West');
ylabel('Deg. North');
ax(4) = subplot(224);
cn = contourf(ax(4), plons, plats, mean(hmem_FA, 3), 21);
xlim(xwin); ylim(ywin);
caxis([lim1 lim2])
colorbar;
ax(4).FontSize = 12;
title('Field-aligned Mean')
xlabel('Deg. West');
	
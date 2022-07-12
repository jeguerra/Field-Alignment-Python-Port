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
tdex = 1;
zdex = 10;
fname = 'reventador1';
fiinput = ['/media/jeguerra/FastDATA/Field-Alignment-Python-Port/test_data/' fname '.nc'];
%----read data---

meml = 1:10;
nmem = length(meml);
lats = (ncread(fiinput, 'latitude'));
nlat = size(lats); nlat = nlat(1);
lons = (ncread(fiinput, 'longitude'));
nlon = size(lons); nlon = nlon(1); 
hmemAll = ncread(fiinput, 'Flight_levelA');
osize = size(hmemAll);

raw_mean = zeros(nlon,nlat);
alg_mean = zeros(nlon,nlat);
pdex = [2,1,3];

%----parameters for field alignment
iter = 5;
lscale = 1;   %change from 1~10,   1 is the best
vmode = 4;  %hyperviscosity model

for tt=1:1
    % Arrange matrix in (lat, lon, ensid)
    hmem0 = squeeze(hmemAll(:,:,:));
    hmemT = permute(hmem0,pdex);

    lats = permute(lats,pdex);
    lons = permute(lons,pdex);

    % Apply windowing function to the data
    win1 = tukeywin(nlat,0.25);
    win2 = tukeywin(nlon,0.25);
    win = (win1 * win2');
    K = 0.5*ones(3);
    for imem = meml
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

    %------------------------
    %==========Main algorithm=========
    hmem_FA = zeros(nlon,nlat,nmem,'double');
    for imem = meml  % reference member
        fprintf('reference mem is:');
        disp(imem);
        qxyT_mem = zeros(nlat,nlon,2,nmem,'double');

        %------derive the displacement vectors-------
        %hmemT(:,:,1),hmemT(:,:,2) ->->-> move 1 to 2
        %hmemT(:,:,i),hmemT(:,:,:) ->->-> move i to 1~20 member
        thisHmemT = hmemT(:,:,imem);
        for i = meml  % reference mem to different members
            if imem ~= i
                thatHmemT = hmemT(:,:,i);
                amess = sprintf(['Aligning: ', int2str(imem), ' to ' int2str(i)]);
                disp(amess);
                [qxyT_mem(:,:,1,i),qxyT_mem(:,:,2,i)] = ...
                    FA2DImNoHLIB(thisHmemT,thatHmemT,0.1,iter,vmode,lscale);
            end
        end

        %------calculate the mean of all DVs----
        qxyT = zeros(nlat,nlon,2,'double');
        qxy = zeros(nlon,nlat,2,'double');

        for i = meml
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
[plons,plats] = meshgrid(lons, lats);

% Visual window
%xwin = [min(min(lons)), max(max(lons))];
ywin = [min(min(lats)), max(max(lats))];
xwin = [-77.8, -77.3];
ywin = [-1.25, 0.25];

% Set a contour level
cl = 5.0E-4 * [1.0, 1.0];

% Plot original data into left axes
ax(1) = subplot(221); cc = 1;
for p = meml
    if cc == 8
        break;
    end
	[Mc, cn] = contour(ax(1), plons, plats, hmem0(:,:,p)', cl, '-');
    cn.LineWidth = 2;
    cn.LineColor = ax(1).ColorOrder(cc,:);
    grid on;
    hold on;
    cc = cc + 1;
end
hold off;
alpha(0.1)
xlim(xwin); ylim(ywin);
ax(1).FontSize = 12;
%xlabel('Deg. West');
ylabel('Deg. North');
title('Initial Ensemble');
	
% Plot aligned data into right axes
ax(2) = subplot(222); cc = 1;
for p = meml
    if cc == 8
        break;
    end
	[Mc, cn] = contour(ax(2), plons, plats, hmem_FA(:,:,p)', cl, '-');
    cn.LineWidth = 2;
    cn.LineColor = ax(2).ColorOrder(cc,:);
    grid on;
    hold on;
    cc = cc + 1;
end
hold off;
alpha(0.1)
xlim(xwin); ylim(ywin);
ax(2).FontSize = 12;
xlabel('Deg. West');
ylabel('Deg. North');
title('Aligned Ensemble');

% Plot means
lim1 = 0.0;
lim2 = 0.0035;
colormap('jet')
ax(3) = subplot(223);
AM = mean(hmem0, 3);
[cn, hn] = contourf(ax(3), plons, plats, AM', 11, 'k');
hn.LineStyle = 'none';
xlim(xwin); ylim(ywin);
caxis([lim1 lim2])
colorbar;
ax(3).FontSize = 12;
title('Arithmetic Mean');
ax(4) = subplot(224);
FM = mean(hmem_FA, 3);
[cn, hn] = contourf(ax(4), plons, plats, FM', 11, 'k');
hn.LineStyle = 'none';
xlim(xwin); ylim(ywin);
caxis([lim1 lim2])
colorbar;
ax(4).FontSize = 12;
title('Field-aligned Mean')
	
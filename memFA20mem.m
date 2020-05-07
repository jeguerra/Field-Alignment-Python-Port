%close all;
clear all;
t1 = clock
%Field alignment
%input:
%---two fields---
%A and B
%output:
%--C--> qx, qy : vector field
%--D align B to A
%input files------------------
%dir = ...    filename = ...
    date=1001;
    fcst='00';
    fiinput = ['/test_data/reventador9.nc'];
    fileID1 = fopen(fiinput);
%----parameters----
    precision = 'float';
    nlon = 162;
    nlat = 354;
    nmem = 13;
%----read data----
    hmem0 = fread(fileID1,nlon*nlat*nmem,precision);   % read the file Yall(360*181*20)
    hmem = reshape(hmem0,[nlon,nlat,nmem]);
%---------transpose the data for FA tool------
    hmemT = zeros(nlat,nlon,nmem,'single');
    for i = 1:1:nmem;
	hmemT(:,:,i) = transpose(hmem(:,:,i));  %matrix transposed before input to subroutine
    end;
   %============================== 
   fprintf('Start Calculation');
   %imem=1;    %++++++++++++++++++++where you need to change+++++++
   %fprintf('This is the member',imem)
%----parameters for field alignment
    iter = 2^10;
    lscale = 1;   %change from 1~10,   1 is the best 
%------------------------
%==========Main algorithm=========
for imem = 1:1:20;  % reference member
    	fprintf('reference mem is:');
        disp(imem);
        qxyT_mem = zeros(nlat,nlon,2,nmem,'single');
%------derive the displacement vectors-------
%hmemT(:,:,1),hmemT(:,:,2) ->->-> move 1 to 2
%hmemT(:,:,i),hmemT(:,:,:) ->->-> move i to 1~20 member
    for i = 1:1:nmem;   % reference mem to different members
    	fprintf('imem is:');
        disp(i);
    	[qxyT_mem(:,:,1,i),qxyT_mem(:,:,2,i)] = FA2DImNoHLIB(hmemT(:,:,imem),hmemT(:,:,i),128,.1,iter,2,lscale);
    end
%------calculate the mean of all DVs----
    qxyT = zeros(nlat,nlon,2,'single');
    qxy = zeros(nlon,nlat,2,'single');
    for i = 1:1:nmem;
    	qxyT(:,:,1)=qxyT(:,:,1)+qxyT_mem(:,:,1,i);
    	qxyT(:,:,2)=qxyT(:,:,2)+qxyT_mem(:,:,2,i);
    end;
    qxyT(:,:,1)=qxyT(:,:,1)./nmem;
    qxyT(:,:,2)=qxyT(:,:,2)./nmem;
    qxy(:,:,1) = transpose(qxyT(:,:,1));    %output after transopose:
    qxy(:,:,2) = transpose(qxyT(:,:,2));    %output after transopose:
    fprintf('Derive the displacement vector successfully!')
%-----align the member i----
    hmemT_FA = zeros(nlat,nlon,'single');
    hmem_FA = zeros(nlon,nlat,'single');
     hmemT_FA = advect(hmemT(:,:,imem),qxyT(:,:,1),qxyT(:,:,2));
     hmem_FA = transpose(hmemT_FA);    %output after transopose:
     fprintf('align member imem successfully!')
%-----output data----
    filename1 = ['/scratch3/BMC/det/Jie.Feng/DTC/memAfterFA_date' int2str(date) '_leadt' fcst 'mem' int2str(imem) '.dat'];
    filename2 = ['/scratch3/BMC/det/Jie.Feng/DTC/memDV_date' int2str(date) '_leadt' fcst 'mem' int2str(imem) '.dat'];
    fileID5 = fopen(filename1,'w');
    fileID6 = fopen(filename2,'w');
	fwrite(fileID5,hmem_FA,precision);
	fwrite(fileID6,qxy,precision);
    fclose(fileID5);
    fclose(fileID6);
end
t2 = clock
%e = etime(t2,t1)

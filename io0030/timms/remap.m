%define paths
%dataPath='/work/ollie/orichter/MisomipPlus/fesom_data/';

oldOceFile='/work/ollie/orichter/MisomipPlus/io0030/fesomdata/io0030.1043.oce.nc';
newOceFile='/work/ollie/orichter/MisomipPlus/io0030/fesomdata/io0030.1043.oce.ini.nc';

oldIceFile='/work/ollie/orichter/MisomipPlus/io0030/fesomdata/io0030.1043.ice.nc';
newIceFile='/work/ollie/orichter/MisomipPlus/io0030/fesomdata/io0030.1043.ice.ini.nc';

oldMeshPath='/work/ollie/orichter/MisomipPlus/io0030/fesommesh/1043.08/';
newMeshPath='/work/ollie/orichter/MisomipPlus/io0030/fesommesh/1043.09/';


% Loading old mesh:
global x2dOld y2dOld x3dOld y3dOld z3dOld nod32Old
disp('Loading old mesh ...')

fid=fopen([oldMeshPath,'nod2d.out']);
n2dOld=fscanf(fid,'%g',1);
nodes2=fscanf(fid,'%g', [4,n2dOld]);
fclose(fid);
x2dOld=transpose(nodes2(2,:));
y2dOld=transpose(nodes2(3,:));
% icoord=nodes(4,:);
clear nodes2

fid=fopen([oldMeshPath,'nod3d.out'],'r');
n3d=fscanf(fid,'%g',1);
nodes3=fscanf(fid,'%g',[5 n3d]);
fclose(fid);
x3dOld = transpose(nodes3(2,:));
y3dOld = transpose(nodes3(3,:));
z3dOld = transpose(nodes3(4,:));
%clear nodes3

fid=fopen([oldMeshPath,'shelf.out']);
shelfOld=fscanf(fid,'%g');
fclose(fid);

fid=fopen([oldMeshPath,'aux3d.out'],'r');
nl=fscanf(fid,'%g',1);
nod32Old=fscanf(fid,'%g', [nl n2dOld]);
fclose(fid);


% Loading new mesh:
disp('Loading new mesh ...')
global x2dNew y2dNew x3dNew y3dNew z3dNew idxOld2New idxOld2New3d n2dNew nod32New mask2d mask3d

fid=fopen([newMeshPath,'nod2d.out']);
n2dNew=fscanf(fid,'%g',1);
nodes2=fscanf(fid,'%g', [4,n2dNew]);
fclose(fid);
x2dNew=transpose(nodes2(2,:));
y2dNew=transpose(nodes2(3,:));
% icoord=nodes(4,:);
clear nodes2

fid=fopen([newMeshPath,'nod3d.out'],'r');
n3dNew=fscanf(fid,'%g',1);
nodes3=fscanf(fid,'%g',[5 n3dNew]);
fclose(fid);
x3dNew = nodes3(2,:);
y3dNew = nodes3(3,:);
z3dNew = nodes3(4,:);
clear nodes3

fid=fopen([newMeshPath,'aux3d.out'],'r');
nl=fscanf(fid,'%g',1);
nod32New=fscanf(fid,'%g', [nl n2dNew]);
nod23New=fscanf(fid,'%g', n3dNew);
fclose(fid);

fid=fopen([newMeshPath,'shelf.out']);
shelfNew=fscanf(fid,'%g');
fclose(fid);


%find 2d indices that map from old grid to new grid and create 2d and 3d
%mask (0=old cell, 1=new cell from new water column aka GL retreat, 2= new
%cell from uplift)

[idxOld2New,D]=knnsearch([x2dOld y2dOld],[x2dNew y2dNew]);
mask2d=zeros(n2dNew,1); mask2d(D>0)=1;

lev_diff = sum(nod32New>0,1) - sum(nod32Old(:,idxOld2New)>0,1);
nod32Old2New = ones(size(nod32New)).*-999;
mask32 = zeros(size(nod32New));
mask32(:,mask2d==1)=1;

for i=1:n2dNew
    ld = lev_diff(i);
    wctOld = nod32Old(:,idxOld2New(i));
    if ld == 0
        nod32Old2New(:,i) = wctOld;
           
    elseif ld > 0 %retreat
        nod32Old2New(:,i) = [ones(ld,1).*wctOld(1); wctOld(1:end-ld)] ; 
        mask32(1:ld,i) = 2;

    elseif ld < 0 %advance
        nod32Old2New(1:end+ld,i) = [wctOld(-ld+1:end)];
    end

    if sum(nod32Old2New(:,i)>0,1) ~= sum(nod32New(:,i)>0,1)
        disp('break at i = ')
        break
    end
end
idxOld2New3d = zeros(n3dNew,1);
idxOld2New3d(nod32New(nod32New>0)) = nod32Old2New(nod32New>0);

mask3d = NaN(size(idxOld2New3d));
mask3d(nod32New(nod32New>0)) = mask32(nod32New>0);

if (any(idxOld2New3d==0)) || (any(isnan(mask3d)))
    disp('mask is fucked')
end

disp(['number of new watercolumns from GL retreat: ',num2str(sum(mask2d==1))])
disp(['number of new cells from GL retreat: ',num2str(sum(mask3d==1))])
disp(['number of new cells from Uplift: ',num2str(sum(mask3d==2))])
%global lev_diff
%lev_diff = sum(nod32New>0,1) - sum(nod32Old(:,idxOld2New)>0,1);
%disp(['maximum retreat in levels: ',num2str(max(lev_diff))])
%disp(['maximum advance in levels: ',num2str(min(lev_diff))])

%setting up ncids
disp 'processing oce file'
global ncid_in ncid_out
ncid_in = netcdf.open(oldOceFile);
ncid_out = netcdf.create(newOceFile,'WRITE');

%setting up dimensions
dimidt = netcdf.defDim(ncid_out,'T',1);
dimidn2 = netcdf.defDim(ncid_out,'nodes_2d',n2dNew);
dimidn3 = netcdf.defDim(ncid_out,'nodes_3d',n3dNew);

time_ID=netcdf.defVar(ncid_out,'time','double',dimidt);
netcdf.putAtt(ncid_out,time_ID,'long_name','model time');
netcdf.putAtt(ncid_out,time_ID,'units','s');

ssh_ID = netcdf.defVar(ncid_out,'ssh','double',[dimidn2 dimidt]);
netcdf.putAtt(ncid_out,ssh_ID,'description','sea surface elevation');
netcdf.putAtt(ncid_out,ssh_ID,'units','m');

u_ID = netcdf.defVar(ncid_out,'u','double',[dimidn3 dimidt]);
netcdf.putAtt(ncid_out,u_ID,'description','zonal velocity');
netcdf.putAtt(ncid_out,u_ID,'units','m/s');

v_ID = netcdf.defVar(ncid_out,'v','double',[dimidn3 dimidt]);
netcdf.putAtt(ncid_out,v_ID,'description','meridional velocity');
netcdf.putAtt(ncid_out,v_ID,'units','m/s');

w_ID = netcdf.defVar(ncid_out,'w','double',[dimidn3 dimidt]);
netcdf.putAtt(ncid_out,w_ID,'description','vertical velocity');
netcdf.putAtt(ncid_out,w_ID,'units','m/s');

wpot_ID = netcdf.defVar(ncid_out,'wpot','double',[dimidn3 dimidt]);
netcdf.putAtt(ncid_out,wpot_ID,'description','vertical velocity potential');
netcdf.putAtt(ncid_out,wpot_ID,'units','m.m/s');

temp_ID = netcdf.defVar(ncid_out,'temp','double',[dimidn3 dimidt]);
netcdf.putAtt(ncid_out,temp_ID,'description','potential temperature');
netcdf.putAtt(ncid_out,temp_ID,'units','degC');

salt_ID = netcdf.defVar(ncid_out,'salt','double',[dimidn3 dimidt]);
netcdf.putAtt(ncid_out,salt_ID,'description','salinity');
netcdf.putAtt(ncid_out,salt_ID,'units','psu');
netcdf.endDef(ncid_out);

%writting in time
time = netcdf.getVar(ncid_in,0); time = time(end);
netcdf.putVar(ncid_out,time_ID,time);

%writing other fields
netcdf.putVar(ncid_out,ssh_ID,getDataOut2D(2,'NN'));
netcdf.putVar(ncid_out,u_ID,getDataOut3D(3,'0','0'));
netcdf.putVar(ncid_out,v_ID,getDataOut3D(4,'0','0'));
netcdf.putVar(ncid_out,w_ID,getDataOut3D(5,'0','0'));
netcdf.putVar(ncid_out,wpot_ID,getDataOut3D(6,'0','0'));
netcdf.putVar(ncid_out,temp_ID,getDataOut3D(7,'NN','NN'));
netcdf.putVar(ncid_out,salt_ID,getDataOut3D(8,'NN','NN'));
%here would need to be the function call with all 4 arrays:
% ...get_data_out(3,array1,array2,weights1,weights2)

netcdf.close(ncid_in);
netcdf.close(ncid_out);

disp 'processing ice file'
ncid_in = netcdf.open(oldIceFile);
ncid_out = netcdf.create(newIceFile,'WRITE');

%setting up dimensions
dimidt = netcdf.defDim(ncid_out,'T',1);
dimidn2 = netcdf.defDim(ncid_out,'nodes_2d',n2dNew);

time_ID=netcdf.defVar(ncid_out,'time','double',dimidt);
netcdf.putAtt(ncid_out,time_ID,'long_name','model time');
netcdf.putAtt(ncid_out,time_ID,'units','s');

area_ID = netcdf.defVar(ncid_out,'area','double',[dimidn2 dimidt]);
netcdf.putAtt(ncid_out,area_ID,'description','ice concentration [0-1]');

hice_ID = netcdf.defVar(ncid_out,'hice','double',[dimidn2 dimidt]);
netcdf.putAtt(ncid_out,hice_ID,'description','effective ice thickness');
netcdf.putAtt(ncid_out,hice_ID,'units','m');

hsnow_ID = netcdf.defVar(ncid_out,'hsnow','double',[dimidn2 dimidt]);
netcdf.putAtt(ncid_out,hsnow_ID,'description','effective snow thickness');
netcdf.putAtt(ncid_out,hsnow_ID,'units','m');

uice_ID = netcdf.defVar(ncid_out,'uice','double',[dimidn2 dimidt]);
netcdf.putAtt(ncid_out,uice_ID,'description','zonal velocity');
netcdf.putAtt(ncid_out,uice_ID,'units','m/s');

vice_ID = netcdf.defVar(ncid_out,'vice','double',[dimidn2 dimidt]);
netcdf.putAtt(ncid_out,vice_ID,'description','meridional velocity');
netcdf.putAtt(ncid_out,vice_ID,'units','m/s');

netcdf.endDef(ncid_out);

%writting in time
time = netcdf.getVar(ncid_in,0); time = time(end);
netcdf.putVar(ncid_out,time_ID,time);

%writing other fields
netcdf.putVar(ncid_out,area_ID,getDataOut2D(2,'0'));
netcdf.putVar(ncid_out,hice_ID,getDataOut2D(3,'0'));
netcdf.putVar(ncid_out,hsnow_ID,getDataOut2D(4,'0'));
netcdf.putVar(ncid_out,uice_ID,getDataOut2D(5,'0'));
netcdf.putVar(ncid_out,vice_ID,getDataOut2D(6,'0'));

netcdf.close(ncid_in);
netcdf.close(ncid_out);

exit;
function data_out = getDataOut2D(var_ind,wcFlag)
    global ncid_in mask2d idxOld2New
    data_in = netcdf.getVar(ncid_in,var_ind); data_in = data_in(:,end);
    data_out = data_in(idxOld2New);
    if wcFlag == '0'
        data_out(mask2d==1)=0;
    end
end

function data_out = getDataOut3D(var_ind,wcFlag,upliftFlag)
    global ncid_in mask3d idxOld2New3d z3dNew  
    data_in = netcdf.getVar(ncid_in,var_ind); data_in = data_in(:,end);
    data_out = data_in(idxOld2New3d);
    if wcFlag == '0'
        data_out(mask3d==1)=0;
    end
    if upliftFlag == '0'
        data_out(mask3d==2)=0;
    elseif upliftFlag == 's'
        %z=transpose(z3dNew(mask3d==2));
        %data_out(mask3d==2)=z.*(-1/720);
        data_out(mask3d==2)=5;
    elseif upliftFlag=='Tf'
        z=transpose(z3dNew(mask3d==2));
        S_in = netcdf.getVar(ncid_in,8); S_in = S_in(:,end);
        S_out = S_in(idxOld2New3d);
        S=S_out(mask3d==2);
        data_out(mask3d==2)=-0.0573.*S+0.0832+7.53e-4.*z;
    end
end
%{  
below only old functions, e.g. for barotropic velocity correction
function data_out = get_data_out(var_ind,idxOld2New)
    global ncid_in
    data_in = netcdf.getVar(ncid_in,var_ind); data_in = data_in(:,end);
    data_out = data_in(idxOld2New);
end

function data_out = get_data_out2(var_ind,idxOld2New)
    global ncid_in nod32Old nod32New n2dNew
    data_in = netcdf.getVar(ncid_in,var_ind); data_in = data_in(:,end);
    data_out = data_in(idxOld2New);
    
    cnt = 0;
    for i=1:n2dNew
        idx = nod32New(:,i);
        idx = idx(idx>0);
        column = data_out(idx);
        if any(diff(column)==0)
            cnt= cnt+1;
            column(find(diff(column)==0)+1)=NaN;
            data_out(idx) = fillmissing(column,'linear','EndValues','nearest');
        end
    end
    disp([num2str(cnt),' water columns have been adjusted to linear.'])
end

    
function data_out = get_velo_out(var_ind)
    global ncid_in z3dOld z3dNew nod32Old nod32New n2dNew idxOld2New idxOld2New3d
    data_in = netcdf.getVar(ncid_in,var_ind); data_in = data_in(:,end);
    data_out = data_in(idxOld2New3d);
    
    lev_diff = sum(nod32New>0,1) - sum(nod32Old(:,idxOld2New)>0,1);
    for i=1:n2dNew
    if lev_diff(i) < 0
       idx =   nod32New(:,i);
       idx = idx(idx>0);
       h_new = z3dNew(idx);
       h_old = z3dOld(idxOld2New3d(idx));
       dh_new = abs(h_new(end) - h_new(1));
       dh_old = abs(h_old(end) - h_old(1));
       u_new = mean(data_out(idx));
       u_old = mean(data_in(idxOld2New3d(idx))); 
       
       du = (u_old*dh_old)/dh_new - u_new;
       
       h_diff = abs(diff(h_new));
       h_diff = [h_diff h_diff(end)];
       weights = transpose(h_diff./(dh_new+h_diff(end)));
     
       data_out(idx)=du*weights+data_out(idx);      
    end
    end
    end

function data_out = get_data_out(var_ind,dim,ncid_in)
    global ncid_in x2d y2d x3d y3d z3d x2dNew y2dNew x3dNew y3dNew z3dNew
    data_in = netcdf.getVar(ncid_in,var_ind); data_in = data_in(:,end);
    if dim == '2d'
        interp = scatteredInterpolant(x2d,y2d,data_in,'nearest','nearest');
        data_out = interp(x2dNew,y2dNew);
    elseif dim == '3d'
        weights = -(max(x2d)*111000)/min(z3d);
        interp = scatteredInterpolant(x3d,y3d,z3d./111000.*weights,data_in,'nearest','nearest');
        data_out = interp(x3dNew,y3dNew,z3dNew);
    end
    end
%}


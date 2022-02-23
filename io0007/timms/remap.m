%function remap(oldMeshPath,newMeshPath,oldOceFile,newOceFile,oldIceFile,newIceFile)
%define paths
%dataPath='/work/ollie/orichter/MisomipPlus/fesom_data/';

oldOceFile='/work/ollie/orichter/MisomipPlus/io0007/fesomdata/io0007.1098.oce.nc';
newOceFile='/work/ollie/orichter/MisomipPlus/io0007/fesomdata/io0007.1098.oce.ini.nc';

oldIceFile='/work/ollie/orichter/MisomipPlus/io0007/fesomdata/io0007.1098.ice.nc';
newIceFile='/work/ollie/orichter/MisomipPlus/io0007/fesomdata/io0007.1098.ice.ini.nc';

oldMeshPath='/work/ollie/orichter/MisomipPlus/io0007/fesommesh/1098/';
newMeshPath='/work/ollie/orichter/MisomipPlus/io0007/fesommesh/1099/';


% Loading old mesh:
global x2d y2d x3d y3d z3d
disp('Loading old mesh ...')

fid=fopen([oldMeshPath,'nod2d.out']);
n2d=fscanf(fid,'%g',1);
nodes2=fscanf(fid,'%g', [4,n2d]);
fclose(fid);
x2d=transpose(nodes2(2,:));
y2d=transpose(nodes2(3,:));
% icoord=nodes(4,:);
clear nodes2

fid=fopen([oldMeshPath,'nod3d.out'],'r');
n3d=fscanf(fid,'%g',1);
nodes3=fscanf(fid,'%g',[5 n3d]);
fclose(fid);
x3d = transpose(nodes3(2,:));
y3d = transpose(nodes3(3,:));
z3d = transpose(nodes3(4,:));
%clear nodes3

% Loading new mesh:
disp('Loading new mesh ...')
global x2dNew y2dNew x3dNew y3dNew z3dNew

fid=fopen([newMeshPath,'nod2d.out']);
n2d=fscanf(fid,'%g',1);
nodes2=fscanf(fid,'%g', [4,n2d]);
fclose(fid);
x2dNew=nodes2(2,:);
y2dNew=nodes2(3,:);
% icoord=nodes(4,:);
%clear nodes2

fid=fopen([newMeshPath,'nod3d.out'],'r');
n3d=fscanf(fid,'%g',1);
nodes3=fscanf(fid,'%g',[5 n3d]);
fclose(fid);
x3dNew = nodes3(2,:);
y3dNew = nodes3(3,:);
z3dNew = nodes3(4,:);
clear nodes3

%setting up ncids
disp 'processing oce file'
global ncid_in ncid_out
ncid_in = netcdf.open(oldOceFile);
ncid_out = netcdf.create(newOceFile,'WRITE');

%setting up dimensions
dimidt = netcdf.defDim(ncid_out,'T',1);
dimidn2 = netcdf.defDim(ncid_out,'nodes_2d',n2d);
dimidn3 = netcdf.defDim(ncid_out,'nodes_3d',n3d);

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
netcdf.putVar(ncid_out,ssh_ID,get_data_out(2,'2d'));
netcdf.putVar(ncid_out,u_ID,get_data_out(3,'3d'));
netcdf.putVar(ncid_out,v_ID,get_data_out(4,'3d'));
netcdf.putVar(ncid_out,w_ID,get_data_out(5,'3d'));
netcdf.putVar(ncid_out,wpot_ID,get_data_out(6,'3d'));
netcdf.putVar(ncid_out,temp_ID,get_data_out(7,'3d'));
netcdf.putVar(ncid_out,salt_ID,get_data_out(8,'3d'));

netcdf.close(ncid_in);
netcdf.close(ncid_out);


disp 'processing ice file'
ncid_in = netcdf.open(oldIceFile);
ncid_out = netcdf.create(newIceFile,'WRITE');

%setting up dimensions
dimidt = netcdf.defDim(ncid_out,'T',1);
dimidn2 = netcdf.defDim(ncid_out,'nodes_2d',n2d);

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
netcdf.putVar(ncid_out,area_ID,get_data_out(2,'2d'));
netcdf.putVar(ncid_out,hice_ID,get_data_out(3,'2d'));
netcdf.putVar(ncid_out,hsnow_ID,get_data_out(4,'2d'));
netcdf.putVar(ncid_out,uice_ID,get_data_out(5,'2d'));
netcdf.putVar(ncid_out,vice_ID,get_data_out(6,'2d'));

netcdf.close(ncid_in);
netcdf.close(ncid_out);

exit;
function data_out = get_data_out(var_ind,dim,ncid_in)
    global ncid_in x2d y2d x3d y3d z3d x2dNew y2dNew x3dNew y3dNew z3dNew
    data_in = netcdf.getVar(ncid_in,var_ind); data_in = data_in(:,end);
    if dim == '2d'
        interp = scatteredInterpolant(x2d,y2d,data_in,'nearest','nearest');
        data_out = interp(x2dNew,y2dNew);
    elseif dim == '3d'
        interp = scatteredInterpolant(x3d,y3d,z3d./111000,data_in,'nearest','nearest');
        data_out = interp(x3dNew,y3dNew,z3dNew./111000);
    end
end

%end

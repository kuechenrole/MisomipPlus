% make cavity
% modify the 3D mesh to remove the elements where ice shelf locates
%
% Qiang, 05.02.2017
%--------------------------------------------------------------------

%close all

%--------------------------------------------------------------------
% user specification starts here:

% paths
meshpath=meshOutPath;%'/work/ollie/orichter/MisomipPlus/fesom_mesh/010/';

%meshpath='/work/ollie/cwekerle/MG/mesh_jigsaw_79NG/mesh_79NG_coarse4_final/'


% Loading mesh:
disp('Loading mesh ...')

fid=fopen([meshpath,'nod2d.out']);
n2d=fscanf(fid,'%g',1);
nodes2=fscanf(fid,'%g', [4,n2d]);
fclose(fid);
% xcoord=nodes(2,:);
% ycoord=nodes(3,:);
% icoord=nodes(4,:);
%clear nodes2

fid=fopen([meshpath, 'elem2d.out']);
e2d=fscanf(fid,'%g', 1);
elem2=fscanf(fid,'%g', [3,e2d]);
fclose(fid);

fid=fopen([meshpath,'nod3d.out'],'r');
n3d=fscanf(fid,'%g',1);
nodes3=fscanf(fid,'%g',[5 n3d]);
fclose(fid);

fid=fopen([meshpath,'elem3d.out'],'r');
e3d=fscanf(fid,'%g',1);
elem3=fscanf(fid,'%g', [4,e3d]);
fclose(fid);

fid=fopen([meshpath,'aux3d.out'],'r');
nl=fscanf(fid,'%g',1);
nod32=fscanf(fid,'%g', [nl n2d]);
nod23=fscanf(fid,'%g', n3d);
elem23=fscanf(fid,'%g',e3d);
fclose(fid);

fid=fopen([meshpath,'depth.out'],'r');
depth=fscanf(fid,'%g',n2d);
fclose(fid);

fid=fopen([meshpath,'shelf.out'],'r');
shelf=fscanf(fid,'%g',n2d);
fclose(fid);

%fid=fopen([meshpath,'cavity_flag_nod2d.out'],'r');
%cflag=fscanf(fid,'%g',n2d);
%fclose(fid);


%
%---------------------------------------------------------
disp('Modifying the mesh ...')

% model grid layers
dep_lev=unique(nodes3(4,:));
dep_lev=sort(dep_lev,'descend');
numlay=length(dep_lev);
% layer center
dep_mid=(dep_lev(1:end-1)+dep_lev(2:end))/2;

% determine the bottom level according to ice shelf bathymetry
% there could be ice shelf nodes extending into the water, but this will not
% cause damage.
sb_lev=ones(1,n2d);
for i=1:n2d
    if (shelf(i)>=0) continue, end
    a=find(dep_mid<=shelf(i));
    sb_lev(i)=a(1);
    if(a(1)==1) 
        disp('Error, ice shelf got zero thickness')
        return
    elseif(isempty(a))
        disp('Error, ice shelf reaches the ocean bottom')
        return       
    end
end

% determine surface level in the cavity region
ss_lev=sb_lev;
for i=1:e2d
    aux=min(sb_lev(elem2(:,i)));
    for j=1:3
        ss_lev(elem2(j,i))=min([ss_lev(elem2(j,i)), aux]);
    end
end

disp('3D nodes ...')
% remove 3D nodes
% find which nodes should be removed
% and correct index
indnod=zeros(1,n3d);
for i=1:n2d
     if (shelf(i)>=0), continue, end
     n=nod32(1,i);
     %nodes3(5,n)=11;
     %if(ss_lev(i)==sb_lev(i)), nodes3(5,n)=10; end
     nodes3(4,n)=dep_lev(ss_lev(i));
     
     for k=2:1:ss_lev(i)
         n=nod32(k,i);
         indnod(n)=1;
     end
     for k=ss_lev(i)+1:1:sb_lev(i)-1
        n=nod32(k,i);
        nodes3(5,n)=21;
        %nodes3(5,n)=20;
     end    
end

% new nodes index
newind=cumsum(indnod);
newind=[1:n3d]-newind;

% remove those nodes (where indnod==1)
% the shelf surface nodes are put into the first 1:n2d nodes
nodes3(:,indnod==1)=[];
nod23(indnod==1)=[];
n3d=size(nodes3,2);
nodes3(1,:)=[1:n3d];

disp('nod32 ...')
% new nod32
nod32n=zeros(numlay,n2d);
nod32n(:)=-999;
checklay=numlay;
for i=1:n2d
    nod32n(1,i)=i;
    lay=1;
    for k=ss_lev(i)+1:numlay
        n=nod32(k,i);
        if(n<0); break; end  
        lay=lay+1;      
        nod32n(lay,i)=newind(n);
    end
    checklay=min(checklay,lay);
end
disp('minimal number of layers:')
disp(num2str(checklay))


% correct indnod and newind
% the corrected array is required to correct elem3d: the surface nodes
% were removed, not the shelf surface nodes
for i=1:n2d
    if(ss_lev(i)>1)
        indnod(i)=1;
        indnod(nod32(ss_lev(i),i))=0;
        newind(nod32(ss_lev(i),i))=i;
    end
end

disp('3D elements ...')
% remove 3D elements
indelem=zeros(1,e3d);
indnod(indnod==1)=NaN;
for i=1:e3d
    if(isnan(sum(indnod(elem3(:,i)))))
        indelem(i)=1;
    else
        elem3(:,i)=newind(elem3(:,i));
    end
end
elem3(:,indelem==1)=[];
e3d=size(elem3,2);
elem23(indelem==1)=[];
    
% "surface nodes" should be the first 1:n2d nodes



%
%----------------------------------------------------------
%save
disp('Saving result ...')

%copyfile([meshpath,'cavity_flag_nod2d.out'],[meshpath,'3d_cavity_mesh_new/cavity_flag_nod2d.out']);
%copyfile([meshpath,'elem2d.out'],[meshpath,'3d_cavity_mesh_new/elem2d.out']);
%copyfile([meshpath,'nod2d.out'],[meshpath,'3d_cavity_mesh_new/nod2d.out']);
%copyfile([meshpath,'depth2d.out'],[meshpath,'3d_cavity_mesh_new/depth2d.out']);

fid = fopen([meshpath,'nod3d.out'],'w');
fprintf(fid,'%9i \n',n3d);
fprintf(fid,'%9i %9.4f %9.4f %9.2f %3i\n',nodes3);
fclose(fid);

fid = fopen([meshpath,'aux3d.out'],'w');
fprintf(fid,'%5i \n',numlay);
fprintf(fid,'%9i \n',nod32n);
fprintf(fid,'%8i \n',nod23);
fprintf(fid,'%8i \n',elem23);
fclose(fid);

fid = fopen([meshpath,'elem3d.out'],'w');
fprintf(fid,'%10i \n',e3d);
fprintf(fid,'%10i %10i %10i %10i \n',elem3);
fclose(fid);

%make cavity flag (1 = in cavity; 0 = open ocean)
cavity_flag = ones(size(shelf));
cavity_flag(nodes3(4,nod32(1,:))==0)=0;


fid = fopen([meshpath,'cavity_flag_nod2d.out'],'w');
fprintf(fid,'%u \n',cavity_flag);
fclose(fid);

disp('New mesh saved!')


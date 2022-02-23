% reorder the nodes
% run this routine after rem_unref_nodes.f90

%  RT took this routine from /csys/nobackup2_CLIDYN/dsidoren/make_mesh
%     on 16.01.2014, modified the file names, and added a procedure for the
%     ice shelt thickness file.


meshdir=meshOutPath;
output_dir=meshOutPath;


fidnod=fopen([meshdir,'nod2d.out']);
n2d=fscanf(fidnod,'%g',1);
nodes=fscanf(fidnod,'%g',[4,n2d]);
nodes=nodes';
fclose(fidnod);

fidtri=fopen([meshdir,'elem2d.out']);
e2d=fscanf(fidtri,'%g',1);
elem=fscanf(fidtri,'%g',[3,e2d]);
elem=elem';
fclose(fidtri);


fid=fopen([meshdir, 'depth.out']);
dep=fscanf(fid,'%g',n2d);
fclose(fid);
    
fid=fopen([meshdir, 'shelf.out']);
shelf=fscanf(fid,'%g',n2d);
fclose(fid);

%-------------------------------------------   


polys=elem;
polys=[elem(:,3),polys];

[m,n]=size(polys);
c = ones(4*m,1);
c(:) = polys;

xy=nodes(:,2:3);
i1 = polys(:,1);
j1 = polys(:,2);
i2 = polys(:,2);
j2 = polys(:,3);
i3 = polys(:,3);
j3 = polys(:,4);
i4 = polys(:,4);
j4 = polys(:,1);
S = sparse([i1;i2;i3],[j1;j2;j3],1);


Sr=symrcm(S);
%spy(S(Sr,Sr))
size(Sr)
R = randperm(length(S));
SR = symrcm(S(R,R));
Sr = R(SR);
% spy(S(Sr,Sr))

l = zeros(length(nodes),1);
for i = 1:length(nodes)
    v = Sr(i);
    l(v) = i;
end
%Change values inside these arrays to new node numbers
elem2=elem;
elem2(:,1) = l(elem(:,1));
elem2(:,2) = l(elem(:,2));
elem2(:,3) = l(elem(:,3));
[ju,l]=sort(elem2(:,1));
elem2=elem2(l,:);

%Change the order of these arrays
nodes2=nodes;
nodes2(:,2:4) = nodes(Sr,2:4);

nodes2=nodes2';
elem2=elem2';
nodes=nodes';
elem=elem';

dep=dep(Sr);
shelf=shelf(Sr);

fid = fopen([output_dir,'nod2d.out'],'w');
fprintf(fid,'%8i \n',n2d);
fprintf(fid,'%8i %9.4f %9.4f %8i\n',nodes2);
fclose(fid);

fid=fopen([output_dir,'elem2d.out'],'w');
fprintf(fid,'%8i \n', e2d);
fprintf(fid,'%8i %8i %8i\n',elem2);
fclose(fid);

fid=fopen([output_dir,'depth.out'], 'w');
fprintf(fid,'%7d \n',dep);
fclose(fid);

fid=fopen([output_dir,'shelf.out'], 'w');
fprintf(fid,'%7d \n',shelf);
fclose(fid);


% visualize:
%figure(9)
%plot(nodes2(2,:), nodes2(3,:),'+b');
%ai=find(nodes2(4,:)==1);


disp('reorder.m succesfully completed')

%clear all

meshdir='./jigsaw/out/';
filename='isomip-MESH.msh';

MeshDir='./fesom2d/';

%Ua_path = '/work/ollie/orichter/MisomipPlus/ua_data/ResultsFiles/0000100-Nodes8013-Ele15780-Tri3-kH1000-MismipPlus-ice1rr_t.mat';


[mesh] = loadmsh([meshdir,filename]);
x=mesh.point.coord(:,1);
y=mesh.point.coord(:,2);
n2d=length(x);
indmesh=zeros(size(x));
elem=mesh.tria3.index(:,1:3);
el2d=size(elem,1);


load(Ua_path);
xUa = MUA.coordinates(:,1);
yUa = MUA.coordinates(:,2);
interp = scatteredInterpolant(xUa,yUa,F.GF.node,'nearest','nearest');
GM = interp(x*111000,y*111000);

rem_nod=GM>0.5;
%rem_nod=GM>0.1;
n2d_new = n2d-sum(rem_nod);
x_new = x(~rem_nod);
y_new = y(~rem_nod);

rem_elem=zeros(el2d,1);
for ii=1:el2d
    nod=elem(ii,:);
    sum_rem_nod = sum(rem_nod(nod));
    
    if sum_rem_nod>0
        rem_elem(ii)=1;
    end
    if (sum_rem_nod==2) | (sum_rem_nod==1)
        indmesh(nod)=1;
    end   
end

elem_new=elem;
elem_new(rem_elem==1,:)=[];
el2d_new=el2d-sum(rem_elem);

cnt=0;
indnod_new=nan(n2d,1);
% for every old n2d node, we save the new node number 
% if node is removed, put nan
for ii=1:n2d
    if rem_nod(ii)==0
        cnt=cnt+1;
        indnod_new(ii)=cnt;
    end
end

%%% correct the node number in elem2d
for ii=1:el2d_new
   elem_new(ii,:)=indnod_new(elem_new(ii,:));
end


x_max=max(x);
y_min=min(y);
y_max=max(y);
ind=find(x==x_max | y==y_min | y==y_max);
indmesh(ind)=1;

indmesh_new=indmesh(~rem_nod);


fid=fopen([MeshDir,'nod2d.out'],'w');
fprintf(fid,'%d\n',n2d_new);
for n=1:n2d_new
     fprintf(fid,'%10d %20.8f %20.8f %3d\n',n,x_new(n),y_new(n),indmesh_new(n));  
end
fclose(fid);
  

fid=fopen([MeshDir,'elem2d.out'],'w');
fprintf(fid,'%8i \n', el2d_new);
for ii=1:el2d_new
  fprintf(fid,'%8i %8i %8i \n',elem_new(ii,:));
end
fclose(fid);

%exit;

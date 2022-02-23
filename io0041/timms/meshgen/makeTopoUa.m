% load mesh data

%clear all
%close all

id='br'; % 'ir', 'br', ''

meshInPath='./fesom2d/';

%Ua_path = '/work/ollie/orichter/MisomipPlus/ua_data/ResultsFiles/0000100-Nodes8013-Ele15780-Tri3-kH1000-MismipPlus-ice1rr_t.mat';
%load(Ua_path);

%meshOutPath = '/work/ollie/orichter/MisomipPlus/fesom_mesh/010/';

%if contains(id,'ir')
%    meshOutPath=[meshOutPath_tmp,'ice_removed/'];
%elseif contains(id,'br')
%    meshOutPath=[meshOutPath_tmp,'bed_removed/'];
%else
%    meshOutPath=[meshOutPath_tmp,'nothing_removed/'];
%end


% two cases are set up to ensure 40 m min water column thickness: _ir =
% remove ice; _br = remove bedrock; set id = '' if no manipulations


disp('load mesh...')

fid=fopen([meshInPath,'nod2d.out'],'r');
n2d=fscanf(fid,'%g',1);
nodes=fscanf(fid, '%g', [4,n2d]);
fclose(fid);
xcoord=nodes(2,:);
ycoord=nodes(3,:);
nodind=nodes(4,:);

fid=fopen([meshInPath,'elem2d.out']);
el2d=fscanf(fid,'%g',1);
elem=fscanf(fid,'%g',[3 el2d]);
fclose(fid);


%convert lat lon to meters
x_out = xcoord*111000;
y_out = ycoord*111000;

%%Load isomip geometry
%geom_path = ['/home/csys/orichter/IsomipPlus/geometry/Ocean1_input_geom_v1.01.nc'];

x_in = MUA.coordinates(:,1);
y_in = MUA.coordinates(:,2);
shelf_in = F.b;
depth_in = F.B;


% alternatively we could use the official function for the bedrock, but it
% doesn't 

shelf = griddata(x_in,y_in,shelf_in,x_out,y_out,'linear');
shelf(x_out>640000)=0.0;
%depth = griddata(x_in,y_in,depth_in,x_out,y_out,'linear');
depth = MismBed(x_out,y_out);


%ISOMIP+ calving criterion: no ice, where ice thickness is less than 100 m
%(ice draft is less than 90 m); taken out to do to smoothing for sigma;
%CAUTION: not applicable in coupled simulations! However, fesom doesn't like ice bewteen 0 and half the uppermost layer thickness. Hence, for coupled simulations we set ice that is thinner than 10 m to 11 m. The ice model has to come back with exact zero befor fesom ice gets zero. 
%shelf(shelf>-90)=0;
shelf(shelf<0 & shelf>-11)=-11;


% now ensure 40m minimum water column thickness
wct=-(depth-shelf);

if contains(id,'ir')
    wct(wct<=40)=40;
    shelf = depth+wct;
elseif contains(id,'br')
    wct(wct<=40)=40;
    depth = shelf-wct;
end

%make cavity flag (1 = in cavity; 0 = open ocean)
cavity_flag = ones(size(shelf));
cavity_flag(shelf==0)=0;


% find lonely elements: Elements where all nodes where sum(nodes_new(4,nod))==3 
% (connected only with 1 edge to the rest of the elements)
ind_rem_elem=zeros(el2d,1);
for ii=1:el2d
    nod=elem(:,ii);

    if sum(nodes(4,nod))==3
        ind_rem_elem(ii)=1;
    end   
end

elem_new=elem;
elem_new(:,ind_rem_elem==1)=[];
el2d_new=el2d-sum(ind_rem_elem);

% set boundary index to 1 for nodes of removed elements
nodes_new =nodes;
ind=find(ind_rem_elem==1)
for ii=1:sum(ind_rem_elem)
    nod=elem(:,ind(ii));
    nodes_new(4,nod)=1;
end

% remove 2d nodes
% dermine nodes that we keep
keep_nod=zeros(n2d,1);
for ii=1:el2d_new
    nod=elem_new(:,ii);
    keep_nod(nod)=1;
end
ind_rem_nodes=find(keep_nod==0);


cnt=0;
indnod_new=nan(n2d,1);
% for every old n2d node, we save the new node number 
% if node is removed, put nan
for ii=1:n2d
    if keep_nod(ii)==1
        cnt=cnt+1;
        indnod_new(ii)=cnt;
    end
end
nodes_new(:,ind_rem_nodes)=[];
n2d_new=n2d-length(ind_rem_nodes);
nodes_new(1,:)=[1:n2d_new];


%%% correct the node number in elem2d
for ii=1:el2d_new
   elem_new(:,ii)=indnod_new(elem_new(:,ii));
end

% remove nodes from depth.out, shelf.out, wct.out and cavity_flag_nod2d.out
depth_new=depth;
cavity_flag_new=cavity_flag;
shelf_new=shelf;

depth_new(ind_rem_nodes)=[];
cavity_flag_new(ind_rem_nodes)=[];
shelf_new(ind_rem_nodes)=[];

wct_new = -(depth_new-shelf_new);




%%% write out the data with the id for topography files
fid = fopen([meshOutPath,'nod2d.out'],'w');
fprintf(fid,'%9i \n',n2d_new);
fprintf(fid,'%9i %9.4f %9.4f %3i\n',nodes_new);
fclose(fid);

fid = fopen([meshOutPath,'elem2d.out'],'w');
fprintf(fid,'%10i \n',el2d_new);
fprintf(fid,'%10i %10i %10i \n',elem_new);
fclose(fid);

fid = fopen([meshOutPath,'depth.out'],'w');
fprintf(fid,'%9.4f \n',depth_new);
fclose(fid);

%fid = fopen([meshOutPath,'cavity_flag_nod2d.out'],'w');
%fprintf(fid,'%u \n',cavity_flag_new);
%fclose(fid);

fid = fopen([meshOutPath,'shelf.out'],'w');
fprintf(fid,'%9.4f \n',shelf_new);
fclose(fid);

fid = fopen([meshOutPath,'wct.out'],'w');
fprintf(fid,'%9.4f \n',wct_new);
fclose(fid);


%exit;

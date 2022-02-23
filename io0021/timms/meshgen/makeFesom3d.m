% ==========================
% 3D mesh generation: Find 
% subdivision of prisms into 
% tetrahedra
% ==========================
% Algorithm: Introduce edges as ordered pairs of nodes
% edge=(n1,n2) with n1<n2.
% Since in each element elem = (n1, n2, n3) there is a node with max
% number, a node with min number and the remaining node with an intermediate 
% number, two edges will be directed to nmax. This means that edges are
% never cyclic on triangles. This is the prequisite to the
% realizability of prism cut. The rest is just the rule: cut vertical faces
% according to edge direction: from bottom to top.
% sergey.danilov@awi.de, July 2016.

%clear all
%close all


% =============
% Load 2d mesh:  
% =============
cyclic=0;
cyclic_length=360;
% meshdir='./';
% output_dir='./';
meshpath=meshOutPath;%'/work/ollie/orichter/MisomipPlus/fesom_mesh/010/';

meshdir=meshpath;
output_dir=meshpath;
do_prelim_output=0;

% the levels:
%zbar=-[0 10 20 30 40 100:100:1000, 1250:250:6000]; % Put correct levels here
%zbar=-[0 10 20 30 40 50 60 70 80 90 100 115 135 160 190 230 280 340 410 490 ...
%       580 680 790 910 1040 1180 1330 1500 1700 1920 2150 2400 2650 2900 3150 ...
%       3400 3650 3900 4150 4400 4650 4900 5150 5400 5650 5900 6450]
   
%    zbar=-[0. 10. 20. 30. 40. 50. 60. 70. 80. 90. 100. 115. 135. 160. ...
%        190.	230. 280. 340. 410. 490. 580. 680. 790. 910. 1040. 1180. ...
%        1330. 1500. 1700. 1920. 2150. 2400. 2650. 2900. 3150. 3400. ...
%        3650. 3900. 4150. 4400. 4650. 4900. 5150. 5400. 5650. 5900. 6150.]
   
% 70 layer
%   zbar=[0  -5 -10 -15 -20 -25 -30 -35 -40 -45 -50 -55 -60 -65 -70 -75 -80 -85 -90 ...
%   -95 -100 -115 -130 -145 -160 -175 -190 -205 -220 -235 -250 -265 -280 ...
%   -295 -310 -330 -350 -370 -390 -410 -430 -470 -520 -580 -680 -790 -910 ...
%   -1040 -1180 -1330 -1500 -1700 -1920 -2150 -2400 -2650 -2900 -3150 -3400 ...
%-3650 -3900 -4150 -4400 -4650 -4900 -5150 -5400 -5650 -6000 -6250]
zbar = [0:-20:-720];

% zbar=[0:-5:-100,...
%     -115:-15:-310,...
%    -330:-20:-600,... 
%    -615:-25:-915 ...
%    -960 -1040 -1180 -1330 -1500 -1700 -1920 -2150 -2400 -2650 -2900 -3150 -3400 ...
% -3650 -3900 -4150 -4400 -4650 -4900 -5150 -5400 -5650 -6000 -6250]
% figure
% plot(zbar,'*g')
% set(gca,'fontsize',14)

% 
% figure
% plot(zbar,'*g')
% set(gca,'fontsize',14)
% saveas2('plots/vertical_levels.png',200)

nl=length(zbar); % Their number
Z=0.5*(zbar(1:nl-1)+zbar(2:nl));   
fid=fopen([meshdir,'elem2d.out']);
e2d=fscanf(fid,'%g',1);
tri=fscanf(fid,'%g',[3,e2d]);
fclose(fid);

fid=fopen([meshdir,'nod2d.out']);
n2d=fscanf(fid,'%g',1);
nodes=fscanf(fid, '%g',[4 n2d]);
xcoord=nodes(2,:);
ycoord=nodes(3,:);
nodind=nodes(4,:);

fid=fopen([meshdir,'depth.out']);
depth=fscanf(fid,'%g',[1,n2d]);
fclose(fid);

% ============
% For visualization
% ============
    xxc=xcoord(tri);
    yyc=ycoord(tri);
    
 if cyclic,
    if 1>2, % either min or max
    xmin=min(xxc);
    x1=xxc(1,:)-xmin;
    x2=xxc(2,:)-xmin; 
    x3=xxc(3,:)-xmin;
    
    ai=find(x1>cyclic_length/2);
    xxc(1,ai)=xxc(1,ai)-cyclic_length;
    ai=find(x2>cyclic_length/2);
    xxc(2,ai)=xxc(2,ai)-cyclic_length;
    ai=find(x3>cyclic_length/2);
    xxc(3,ai)=xxc(3,ai)-cyclic_length;
    else
    xmax=max(xxc);
    x1=xxc(1,:)-xmax;
    x2=xxc(2,:)-xmax; 
    x3=xxc(3,:)-xmax;
    
    ai=find(x1<-cyclic_length/2);
    xxc(1,ai)=xxc(1,ai)+cyclic_length;
    ai=find(x2<-cyclic_length/2);
    xxc(2,ai)=xxc(2,ai)+cyclic_length;
    ai=find(x3<-cyclic_length/2);
    xxc(3,ai)=xxc(3,ai)+cyclic_length;
    end;
 end;  

tri_ordered=zeros(size(tri));

% =================
% Order nodes in triangles
% =================
for elem=1: e2d,
  elnodes =tri(:,elem); 
  if elnodes(1)>elnodes(2),
     n1=elnodes(2); 
     elnodes(2)=elnodes(1); elnodes(1)=n1;
  end;
  if elnodes(2)>elnodes(3),
     n1=elnodes(3); 
     elnodes(3)=elnodes(2); elnodes(2)=n1; 
  end; % the largest is at the third position
  if elnodes(1)>elnodes(2),
     n1=elnodes(2); 
     elnodes(2)=elnodes(1); elnodes(1)=n1;
  end;
 tri_ordered(:,elem)=elnodes;
end;
% Assume that 
% face diagonals meet at the top of the prism at nmax (tri_ordered(3,elem)),
% face diagonals meet at the bottom under nmin (tri_ordered(1,elem)), and  
% at nn (tri_ordered(2,elem) one diagonal is at the top 
% and the other one is at the bottom. The rule to form tetrahedra
% is: t1=[nmax top, all bottom]; t2=[all top, nmin]; t3=[nn and nmax top,
% nn and nmin at the bottom];
    
% ==========
% Find levels
% ==========

elevels=zeros([1,e2d]);    % the number of levels on elements
nlevels=zeros([1,n2d]);    % the number of levels on nodes

for elem=1: e2d,
    nodes=tri(:,elem);
    %dmean=min(depth(nodes));
    dmean=max(depth(nodes)); % Ralph
    exit_flag=0;
          for nz=1:nl-1,
                 if Z(nz)<dmean, 
                    exit_flag=1;
                    elevels(elem)=nz;
                    break
                 end;
          end;
          if exit_flag==0 & dmean<0, elevels(elem)=nl; end;
          if dmean>=0, elevels(elem)=3; end;
          if elevels(elem)<=2, elevels(elem)=3; end;
end;
% figure(2)
% patch(xxc,yyc,elevels)
 %caxis([0 10])
% axis('equal')
% colorbar
% title('Levels on elements')
% trim_bottom
% We can remove prisms that go into the land (there is no use of them
% for velocity is zero).
% script trim_bottom does it, but requires to assemble many additional 
% arrays (all are 2d)


% =============
% The number of levels at nodes is the largest of the 
% numver of levels at elements that contain this node
% =============
for n=1:e2d,
    for j=1:3,
       node=tri(j,n);
       if nlevels(node)<elevels(n),
       nlevels(node)=elevels(n);
       end;
    end;
end;
% figure(3)
%  pp=patch(xxc,yyc,nlevels(tri));
%  set(pp,'EdgeColor','none');
%  axis('equal')
  %caxis([0 10])
%  title('Levels at nodes')
%  colorbar

%  saveas2([meshdir,'levels_at_node.png'],100)

% ====================  
% Save tri_ordered and level arrays:
% ====================
if do_prelim_output==1,
fid=fopen([output_dir,'elem2d_ord.out'],'w');
        fprintf(fid,'%10i \n', e2d);
        fprintf(fid,'%10i %10i %10i\n',tri_ordered);
        fclose(fid);

fid=fopen([output_dir,'elem_levels.out'],'w');
        fprintf(fid,'%4i\n',elevels);
        fclose(fid);

fid=fopen([output_dir,'node_levels.out'],'w');
        fprintf(fid,'%4i\n',nlevels);
        fclose(fid);
end;
% ====================
% In principle this information is already sufficient 
% to build the rest
% ====================
%
% node 3d below node 2d. This is the main array, it sets numbers
%
nod3d_below=-999*ones([nl,n2d]);   % If this array is too large
                                   % one can write it directly to the file.
                                   % However, we should manage with n2d
                                   % about 5M
nod3d_below(1,:)=1:n2d;
count1=n2d;
for n=1:n2d,
    nod3d_below(2:nlevels(n),n)=count1+[1:nlevels(n)-1];
    count1=count1+nlevels(n)-1;
end;

% ================
% Writing aux3d
% ================
disp('Writing aux3d')
fidaux=fopen([output_dir,'aux3d.out'],'w');
    fprintf(fidaux,'%10i\n', nl);
    fprintf(fidaux,'%10i\n',nod3d_below);
        
%
% nod2d corresponds to nod 3d
%
fprintf(fidaux,'%10i\n',[1:n2d]);
for n=1:n2d,
aux1=n*ones([1,nlevels(n)-1]);
fprintf(fidaux,'%10i\n',aux1);
end; 
%
% elem2d corresp to elem3d
for elem=1:e2d,
    fprintf(fidaux,'%10i\n', elem*ones([1,3*(elevels(elem)-1)]));
end;

fclose(fidaux);      

% ===============
% Writing nod3d:
% =============== 

disp('Writing nod3d')    
fid3d=fopen([output_dir,'nod3d.out'],'w');
n3d=sum(nlevels);
fprintf(fid3d, '%10i\n',n3d);
for n=1:n2d,
    ind=10; if nodind(n)==1, ind=11; end;
    aux2=[n, xcoord(n), ycoord(n), zbar(1), ind];
    fprintf(fid3d,'%10i %8.4f %8.4f %8.4f %4i\n',aux2);
end;
for n=1:n2d,
    aux1=zeros([5,nlevels(n)-1]);
    aux1(1,:) = nod3d_below(2:nlevels(n),n)';
    aux1(2,:) = xcoord(n)*ones([1,nlevels(n)-1]);
    aux1(3,:) = ycoord(n)*ones([1,nlevels(n)-1]);
    aux1(4,:) = zbar(2:nlevels(n));
    if nodind(n)==0, 
    aux1(5,:) = 20*ones([1,nlevels(n)-1]);
    aux1(5,nlevels(n)-1)=30;
    end;
    if nodind(n)==1, 
    aux1(5,:) = 21*ones([1,nlevels(n)-1]);
    aux1(5,nlevels(n)-1)=31;
    end;
fprintf(fid3d,'%10i %8.4f %8.4f %8.4f %4i\n',aux1);
end;
% Check the rules for indices and output format for coordinates
fclose(fid3d);
% ============
% Writting elem3d
% ============
disp('Writing elem3d')
fid3d=fopen([output_dir,'elem3d.out'],'w');
e3d=(sum(elevels)-e2d)*3;
fprintf(fid3d, '%9i\n',e3d);
for elem=1:e2d,      
   for nz=1:elevels(elem)-1,
       t1=[nod3d_below(nz,tri_ordered(3,elem)), nod3d_below(nz+1,tri_ordered(:,elem))];
       t2=[nod3d_below(nz,tri_ordered(:,elem)), nod3d_below(nz+1,tri_ordered(1,elem))];
       t3=[nod3d_below(nz,tri_ordered(2:3,elem)), nod3d_below(nz+1,tri_ordered(1:2,elem))];
       fprintf(fid3d, '%10i %10i %10i %10i\n',[t1,t2,t3]);    
   end;
end;
fclose(fid3d);

disp('writing finished!!!!')
% ====================
% We can consider three possibilities:
% (i) Generate files in the old fashion. 
% (ii) Since files are large and we do not really need ascii, we can write
% directly to netcdf
% (iii) Having elevels and nlevels the rest can be generated in the code.
% (ii) and (iii) imply some additional but elementary work.
% The largest allocated array is nod2d*nl in size. For n2d=5e+5 
% and nl=100 it takes 0.4G memory. The rest is negligible unless bottom 
% trimming is done. Yet even in this case the memory demand should not be 
% exceeding twice the amount above. It should therefore be possible to carry
% the generation in matlab, without translating it to fortran.


      
      
      
      

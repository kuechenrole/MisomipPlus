% ensure_minLayers_ovl.m in place of ensure_pos_wct.m for z-layer models
% ensure minimum number of layers and vertical overlap of adjoining 2d elems
% WARNING: takes forever for large grids!
% adapted from ensure_pos_wct.m by Verena Haid, Nov 2020
% reduced to contain only the ensure-overlap part by Ole Richter, May 2021

%clear all
%close all

meshpath=meshOutPath;
outpath=meshOutPath;

nodfile          =[meshpath,'nod2d.out'];
elemfile         =[meshpath,'elem2d.out'];
depfile_in       =[meshpath,'depth.out'];
shelffile_in     =[meshpath,'shelf.out'];

depfile_out      =[outpath,'depth.out'];
shelffile_out    =[outpath,'shelf.out'];

ovlmin        =  1.      % minimum number of 3d elem to overlap vertically


% make sure these are the same depth levels as for 3d step:
% zbar=-[   0.0000    5.0000   10.0000   20.0000   30.0000   40.0000   50.0000   60.0000 ...
%          70.0000   80.0000   90.0000  100.0000  111.7500  124.5000  137.7500  150.5000 ...
%         162.7500  175.5000  189.0000  203.5000  219.0000  235.5000  253.0000  271.5000 ...
%         291.0000  311.0000  331.0000  351.0000  371.0000  391.0000  411.2500  431.7500 ...
%         452.2500  472.7500  493.2500  513.7500  534.2500  554.7500  575.2500  595.7500 ...
%         616.2500  636.7500  657.2500  678.0000  699.0000  720.2500  741.7500  763.5000 ...
%         785.5000  807.7500  830.2500  853.0000  876.0000  899.5000  924.3000  950.3000 ...
%         977.3000 1005.3000 1034.3000 1064.3000 1095.8000 1129.3000 1164.8000 1202.3000 ...
%        1241.3000 1281.3000 1322.8000 1365.3000 1408.8000 1454.3000 1501.8000 1551.3000 ...
%        1602.8000 1657.3000 1715.6000 1778.1000 1845.6000 1920.6000 2005.6000 2100.6000 ...
%        2205.6000 2320.6000 2445.6000 2580.6000 2725.6000 2880.6000 3045.6000 3220.6000 ...
%        3405.6000 3600.6000 3806.9000 4025.6000 4256.9000 4500.6000 4756.9000 5025.6000 ...
%        5313.1000 5625.6000 5963.1000 6325.6000 6713.1000 ] ;
%

zbar=[0:-20:-720];

%Ralph Timmermann 28.02.11
%--------------------------------------------------------------------
% no step                no step                    no step

fid=fopen(nodfile);
n2d=fscanf(fid,'%g',1);
nodes=fscanf(fid,'%g', [4,n2d]);
fclose(fid);
xcoord=nodes(2,:);
ycoord=nodes(3,:);
icoord=nodes(4,:);

fid=fopen(elemfile);
e2d=fscanf(fid,'%g', 1);
elem=fscanf(fid,'%g', [3,e2d]);
fclose(fid);


fid=fopen(depfile_in);
dep=fscanf(fid,'%g',n2d);
fclose(fid);

fid=fopen(shelffile_in);
shelf=fscanf(fid,'%g',n2d);
fclose(fid);

cavity_flag=shelf*0.;
cavity_flag(find(shelf<0))=1;
%fid=fopen([meshpath,'cavity_flag_nod2d.out'],'r');
%cavity_flag=fscanf(fid,'%g',n2d);
%fclose(fid);

dep_mid=(zbar(1:end-1)+zbar(2:end))/2;
for ei=1:e2d
    % find uppermost level
    ai=find(dep_mid<=min(shelf(elem(:,ei))));
    ul(ei)=ai(1);
    %minbl=ul(ei)+levmin-1; % minimum number for bottom layer
    % find bottom level
    if max(dep(elem(:,ei)))<=dep_mid(end)
        bl(ei)=length(dep_mid);
        %   disp(['hit bottom: ',num2str([max(dep(elem(:,ei)))])])
    else
        au=find(dep_mid<max(dep(elem(:,ei))));
        bl(ei)=au(1)-1;
    end
    %bli(ei)=bl(ei); % for debugging
    %if bl(ei)<minbl
        %for i=1:3
            %if dep(elem(i,ei))>dep_mid(minbl)-1 % if clause for counting
            %      disp([num2str(i) ' ' num2str(ei) ' ' num2str(elem(i,ei)) ' ' num2str(ul(ei)) ' ' num2str(bli(ei)) ' '...
            %           num2str(minbl) ' ' num2str(shelf(elem(:,ei))') ' ' num2str(dep(elem(i,ei))) ' ' num2str(dep_mid(minbl)-1)])
            %    dep(elem(i,ei))=dep_mid(minbl)-1;
            %    cnt=cnt+1;
            %end
            % dep(elem(i,ei))=min(dep(elem(i,ei)),dep_mid(minbl)-1);
        %end
    %    bl(ei)=minbl;
    %end
end

disp('check for overlap')     %   VH
cnt=0;
changednodes=[];
maxct=0;
checked_flag=0*bl;
for ei=1:e2d
    if sum(cavity_flag(elem(:,ei)))==3% & sum(xcoord(elem(:,ei))>=4.2)==3
        if mod(ei,20)==0 disp([num2str(ei),'/',num2str(e2d)]), end
        cnt2=0;
        % find neighbor elems
        for ei2=1:e2d
            if sum(cavity_flag(elem(:,ei2)))==3 & checked_flag(ei2)==0
                
                % schnelle version dieser Abfrage?
                if sum([ismember(elem(1,ei),elem(:,ei2)),ismember(elem(2,ei),elem(:,ei2)),ismember(elem(3,ei),elem(:,ei2))])==2
                    % [ei, ei2]
                    
                    cnt2=cnt2+1;
                    % check overlap
                    % disp([ul(ei) bl(ei) ul(ei2) bl(ei2)])
                    if bl(ei)<ul(ei2)+ovlmin-1
                        for i=1:3
                            %  disp([i shelf(elem(i,ei)) dep(elem(i,ei)) dep_mid(ul(ei2)) dep_mid(ul(ei2)+ovlmin-1)-1 ])
                            if dep(elem(i,ei))>dep_mid(ul(ei2)+ovlmin-1)-1 % if clause for counting
                                dep(elem(i,ei))=dep_mid(ul(ei2)+ovlmin-1)-1;
                                cnt=cnt+1;
                                %   disp('changed')
                            end
                            % dep(elem(i,ei))=min(dep(elem(i,ei)),dep_mid(ul(ei2)+ovlmin-1)-1);
                        end
                    elseif bl(ei2)<ul(ei)+ovlmin-1
                        for i=1:3
                            %  disp([i shelf(elem(i,ei)) dep(elem(i,ei2)) dep_mid(ul(ei2)) dep_mid(ul(ei)+ovlmin-1)-1 ])
                            if dep(elem(i,ei2))>dep_mid(ul(ei)+ovlmin-1)-1 % if clause for counting
                                dep(elem(i,ei2))=dep_mid(ul(ei)+ovlmin-1)-1;
                                cnt=cnt+1;
                                %    disp('changed')
                            end
                            % dep(elem(i,ei2))=min(dep(elem(i,ei2)),dep_mid(ul(ei)+ovlmin-1)-1);
                        end
                    end
                end
            end
            maxct=max(maxct,cnt2);
        end
    end
    checked_flag(ei)=1;
end
disp([' nodes deepened in ',num2str(cnt),' operations to enforce ',num2str(ovlmin),' overlapping layers.'])
%toc

disp('write data')
fid=fopen(depfile_out,'w');
fprintf(fid,'%8i\n', dep);
fclose(fid);

disp('creation of cavity_flag_nod2d_postwct.out succesfully completed')

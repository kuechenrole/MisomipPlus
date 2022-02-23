 function makeJigsaw

    
    %path(path,'~/IsomipPlus/jigsaw/mesh-util/')
    path(path,'/work/ollie/orichter/MisomipPlus/jigsaw/mesh-util/')
 
    close all ;
  
    opts.geom_file = ...
        ['jigsaw/geo/isomip-GEOM.msh'];
    
    opts.jcfg_file = ...
        ['jigsaw/out/isomip.jig'];
    
    opts.mesh_file = ...
        ['jigsaw/out/isomip-MESH.msh'];
    
    opts.hfun_file = ...
        ['jigsaw/out/isomip-HFUN.msh'];
    
%------------------------------------ define JIGSAW geometry

    geom.mshID = 'EUCLIDEAN-MESH';
    
    %M = load('polygon_coords.mat');
    %geom.point.coord = M.coords;
    
    %M = load('polygon_index.mat');
    %geom.edge2.index = M.index;
    
    lon_max = 800/111;
    lat_max = 80/111;

    geom.point.coord = [    % list of xy "node" coordinates
        0, 0, 0             % outer square
        lon_max, 0, 0
        lon_max, lat_max, 0
        0, lat_max, 0 ] ;

    
    geom.edge2.index = [    % list of "edges" between nodes
        1, 2, 0             % outer square 
        2, 3, 0
        3, 4, 0
        4, 1, 0 ] ;

        
    savemsh(opts.geom_file,geom) ;
    
%------------------------------------ compute HFUN over GEOM
    
    XPOS = linspace(0,max(geom.point.coord(:,1)),20);
    YPOS = linspace(0,max(geom.point.coord(:,2)),20);
    
    HFUN = 0.02 * ones(length(YPOS),length(XPOS)) ;
    
    hfun.mshID = 'EUCLIDEAN-GRID';
    hfun.point.coord{1} = XPOS ;
    hfun.point.coord{2} = YPOS ;
    hfun.value = HFUN ;
    
    savemsh (opts.hfun_file,hfun);
    
%------------------------------------ build mesh via JIGSAW! 
    
    opts.hfun_scal = 'absolute';
    opts.hfun_hmax = +inf ;
    opts.hfun_hmin = +0.0 ;
    
    opts.mesh_dims = +2 ;
   %opts.optm_iter = +0 ;
   
   opts.geom_feat = true ; % smooth the coastline
   opts.mesh_top1 = true ; % include small channels 
    
    mesh = jigsaw  (opts) ;
    
    %plotplanar(geom,mesh,hfun) ;
    savemsh (opts.mesh_file,mesh);

   % exit
end



function plotsphere(mesh,hfun)
%PLOT-SPHERE draw JIGSAW output for sphere problems.

    topo = loadmsh('jigsaw/geo/topo.msh');
    
    xpos = topo.point.coord{1};
    ypos = topo.point.coord{2};
    zlev = reshape( ...
    topo.value,length(ypos),length(xpos));

    tlev = ...
        findalt(mesh,xpos,ypos,zlev) ;

    figure;
    surf(hfun.point.coord{1}*180/pi, ...
         hfun.point.coord{2}*180/pi, ...
         hfun.value) ;
    view(2); axis image; hold on ;
    shading interp;
    title('JIGSAW HFUN data') ;
    
    twet = tlev <= +0. ;
    tdry = tlev >  +0. ;
    
    figure;
    patch ('faces',mesh.tria3.index(:,1:3), ...
        'vertices',mesh.point.coord(:,1:3), ...
        'facevertexcdata',tlev(:,:), ...
        'facecolor','flat', ...
        'edgecolor',[.2,.2,.2]) ;
    hold on; axis image off;
    set(gca,'clipping','off') ;
    caxis([min(zlev(:))*4./3., +0.]);
    colormap('hot');
    brighten(+0.75);
    title('JIGSAW TRIA mesh') ;
    
    figure;
    patch ('faces',mesh.tria3.index(twet,1:3), ...
        'vertices',mesh.point.coord(:,1:3), ...
        'facevertexcdata',tlev(twet,:), ...
        'facecolor','flat', ...
        'edgecolor',[.2,.2,.2]) ;
    hold on; axis image off;
    patch ('faces',mesh.tria3.index(tdry,1:3), ...
        'vertices',mesh.point.coord(:,1:3), ...
        'facecolor','w', ...
        'edgecolor','none');
    set(gca,'clipping','off') ;
    caxis([min(zlev(:))*4./3., +0.]);
    colormap('hot');
    brighten(+0.75);
    title('JIGSAW TRIA mesh') ;
    
    drawscr(mesh.point.coord (:,1:3), ...
            [], ...
            mesh.tria3.index (:,1:3)) ;
            
    drawnow ;        
    set(figure(1),'units','normalized', ...
        'position',[.35,.55,.30,.35]) ;
    set(figure(2),'units','normalized', ...
        'position',[.05,.55,.30,.35]) ;
    set(figure(3),'units','normalized', ...
        'position',[.05,.10,.30,.35]) ;
    set(figure(4),'units','normalized', ...
        'position',[.35,.10,.30,.35]) ;
    drawnow ;

end

function plotplanar(geom,mesh,hfun)
%PLOT-PLANAR draw JIGSAW output for planar problems.
  
    figure;
    patch ('faces',geom.edge2.index(:,1:2), ...
        'vertices',geom.point.coord(:,1:2), ...
        'facecolor','w', ...
        'edgecolor',[.1,.1,.1], ...
        'linewidth',1.5) ;
    hold on; axis image;
    title('JIGSAW GEOM data') ;

    figure;
    surf(hfun.point.coord{1}, ...
         hfun.point.coord{2}, ...
         hfun.value) ;
    view(2); axis image; hold on ;
    shading interp;
    title('JIGSAW HFUN data') ;

    figure;
    patch ('faces',mesh.tria3.index(:,1:3), ...
        'vertices',mesh.point.coord(:,1:2), ...
        'facecolor','w', ...
        'edgecolor',[.2,.2,.2]) ;
    hold on; axis image;
    patch ('faces',mesh.edge2.index(:,1:2), ...
        'vertices',mesh.point.coord(:,1:2), ...
        'facecolor','w', ...
        'edgecolor',[.1,.1,.1], ...
        'linewidth',1.5) ;
    patch ('faces',geom.edge2.index(:,1:2), ...
        'vertices',geom.point.coord(:,1:2), ...
        'facecolor','w', ...
        'edgecolor',[.1,.1,.8], ...
        'linewidth',1.5) ;
    title('JIGSAW TRIA mesh') ;

    drawscr(mesh.point.coord (:,1:2), ...
            mesh.edge2.index (:,1:2), ...
            mesh.tria3.index (:,1:3)) ;
    
    drawnow ;        
    set(figure(1),'units','normalized', ...
        'position',[.05,.55,.30,.35]) ;
    set(figure(2),'units','normalized', ...
        'position',[.35,.55,.30,.35]) ;
    set(figure(3),'units','normalized', ...
        'position',[.35,.10,.30,.35]) ;
    set(figure(4),'units','normalized', ...
        'position',[.05,.10,.30,.35]) ;
    drawnow ;

end

function [zlev] = findalt(mesh,alon,alat,topo)
%FINDALT calc. an "altitude" for each tria-cell in the mesh.

    xrad = mesh.point.coord(:,1) .^ 2 ...
         + mesh.point.coord(:,2) .^ 2 ...
         + mesh.point.coord(:,3) .^ 2 ;
    xrad = max(sqrt(xrad),eps) ;
    
    xlat = asin (mesh.point.coord(:,3)./xrad);
    xlon = atan2(mesh.point.coord(:,2), ...
                 mesh.point.coord(:,1)) ;
                 
    xlat = xlat * 180 / pi;
    xlon = xlon * 180 / pi;
    
    xlev = interp2 (alon,alat,topo,xlon,xlat);
    
    zlev = xlev (mesh.tria3.index(:,1)) ...
         + xlev (mesh.tria3.index(:,2)) ...
         + xlev (mesh.tria3.index(:,3)) ;
    zlev = zlev / +3.;

end



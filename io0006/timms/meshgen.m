meshOutPath='/work/ollie/orichter/MisomipPlus/io0006/fesommesh/1036.10/';
Ua_path='/work/ollie/orichter/MisomipPlus/io0006/uadata/1036.10-Nodes8573-Ele16897-Tri3-kH1000-MismipPlus-io0006.mat';
goodfile_path='/work/ollie/orichter/MisomipPlus/io0006/fesommesh/meshgen.goodfile.1036.10';

cd meshgen;

disp('jigsaw2fesom');
jigsaw2fesomUa;

disp('makeTopoUa');
makeTopoUa;

disp('ensure min overlap');
ensure_minLayers_ovl;

disp('reorder');
reorder;

disp('makeFesom3d');
makeFesom3d;

pause(2);
disp('makeCavity');
makeCavity;

cd ..;

fid = fopen(goodfile_path,'w');
fid = fclose(fid);

exit;

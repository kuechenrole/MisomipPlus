meshOutPath='/work/ollie/orichter/MisomipPlus/io0009/fesommesh/1007.70/';
Ua_path='/work/ollie/orichter/MisomipPlus/io0009/uadata/1007.70-Nodes8363-Ele16480-Tri3-kH1000-MismipPlus-io0009.mat';
goodfile_path='/work/ollie/orichter/MisomipPlus/io0009/fesommesh/meshgen.goodfile.1007.70';

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

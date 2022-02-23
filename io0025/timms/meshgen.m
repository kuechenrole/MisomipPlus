meshOutPath='/work/ollie/orichter/MisomipPlus/io0025/fesommesh/1042.30/';
Ua_path='/work/ollie/orichter/MisomipPlus/io0025/uadata/1042.30-Nodes8294-Ele16342-Tri3-kH1000-MismipPlus-io0025.mat';
goodfile_path='/work/ollie/orichter/MisomipPlus/io0025/fesommesh/meshgen.goodfile.1042.30';

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

meshOutPath='/work/ollie/orichter/MisomipPlus/io0011/fesommesh/1036.20/';
Ua_path='/work/ollie/orichter/MisomipPlus/io0011/uadata/1036.20-Nodes8632-Ele17012-Tri3-kH1000-MismipPlus-io0011.mat';
goodfile_path='/work/ollie/orichter/MisomipPlus/io0011/fesommesh/meshgen.goodfile.1036.20';
UaSourcePath='/home/ollie/orichter/MisomipPlus/io0011/ua/UaSource';

addpath(UaSourcePath);
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

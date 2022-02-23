meshOutPath='/work/ollie/orichter/MisomipPlus/io0042/fesommesh/1042.10/';
Ua_path='/work/ollie/orichter/MisomipPlus/io0042/uadata/1042.10-Nodes8169-Ele16096-Tri3-kH1000-MismipPlus-io0042.mat';
goodfile_path='/work/ollie/orichter/MisomipPlus/io0042/fesommesh/meshgen.goodfile.1042.10';

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

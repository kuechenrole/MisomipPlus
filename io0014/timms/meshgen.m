meshOutPath='/work/ollie/orichter/MisomipPlus/io0014/fesommesh/1071.10/';
Ua_path='/work/ollie/orichter/MisomipPlus/io0014/uadata/1071.10-Nodes9563-Ele18868-Tri3-kH1000-MismipPlus-io0014.mat';
goodfile_path='/work/ollie/orichter/MisomipPlus/io0014/fesommesh/meshgen.goodfile.1071.10';

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

meshOutPath='/work/ollie/orichter/MisomipPlus/io0003/fesommesh/1099/';
Ua_path='/work/ollie/orichter/MisomipPlus/io0003/uadata/0109900-Nodes8163-Ele16069-Tri3-kH1000-MismipPlus-io0003.mat';
goodfile_path='/work/ollie/orichter/MisomipPlus/io0003/fesommesh/meshgen.goodfile.1099';

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

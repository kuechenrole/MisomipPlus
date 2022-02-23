meshOutPath='/work/ollie/orichter/MisomipPlus/io0007/fesommesh/1099/';
Ua_path='/work/ollie/orichter/MisomipPlus/io0007/uadata/0109900-Nodes9136-Ele18014-Tri3-kH1000-MismipPlus-io0007.mat';
goodfile_path='/work/ollie/orichter/MisomipPlus/io0007/fesommesh/meshgen.goodfile.1099';

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

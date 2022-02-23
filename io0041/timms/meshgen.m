meshOutPath='/work/ollie/orichter/MisomipPlus/io0041/fesommesh/1199.11/';
Ua_path='/work/ollie/orichter/MisomipPlus/io0041/uadata/1199.11-Nodes15293-Ele30331-Tri3-kH1000-MismipPlus-io0041.mat';
goodfile_path='/work/ollie/orichter/MisomipPlus/io0041/fesommesh/meshgen.goodfile.1199.11';

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

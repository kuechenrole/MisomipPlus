meshOutPath='/work/ollie/orichter/MisomipPlus/io0040/fesommesh/1199.11/';
Ua_path='/work/ollie/orichter/MisomipPlus/io0040/uadata/1199.11-Nodes15572-Ele30884-Tri3-kH1000-MismipPlus-io0040.mat';
goodfile_path='/work/ollie/orichter/MisomipPlus/io0040/fesommesh/meshgen.goodfile.1199.11';

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

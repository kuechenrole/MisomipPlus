meshOutPath='/work/ollie/orichter/MisomipPlus/io0032/fesommesh/1034.11/';
Ua_path='/work/ollie/orichter/MisomipPlus/io0032/uadata/1034.11-Nodes6696-Ele13142-Tri3-kH1000-MismipPlus-io0032.mat';
goodfile_path='/work/ollie/orichter/MisomipPlus/io0032/fesommesh/meshgen.goodfile.1034.11';

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

meshOutPath='/work/ollie/orichter/MisomipPlus/io0015/fesommesh/1065.40/';
Ua_path='/work/ollie/orichter/MisomipPlus/io0015/uadata/1065.40-Nodes9572-Ele18889-Tri3-kH1000-MismipPlus-io0015.mat';
goodfile_path='/work/ollie/orichter/MisomipPlus/io0015/fesommesh/meshgen.goodfile.1065.40';

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

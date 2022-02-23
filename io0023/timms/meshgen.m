meshOutPath='/work/ollie/orichter/MisomipPlus/io0023/fesommesh/1054.60/';
Ua_path='/work/ollie/orichter/MisomipPlus/io0023/uadata/1054.60-Nodes9142-Ele18032-Tri3-kH1000-MismipPlus-io0023.mat';
goodfile_path='/work/ollie/orichter/MisomipPlus/io0023/fesommesh/meshgen.goodfile.1054.60';

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

meshOutPath='/work/ollie/orichter/MisomipPlus/io0012/fesommesh/1199.90/';
Ua_path='/work/ollie/orichter/MisomipPlus/io0012/uadata/1199.90-Nodes9551-Ele18844-Tri3-kH1000-MismipPlus-io0012.mat';
goodfile_path='/work/ollie/orichter/MisomipPlus/io0012/fesommesh/meshgen.goodfile.1199.90';

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

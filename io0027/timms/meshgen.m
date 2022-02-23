meshOutPath='/work/ollie/orichter/MisomipPlus/io0027/fesommesh/1034.90/';
Ua_path='/work/ollie/orichter/MisomipPlus/io0027/uadata/1034.90-Nodes8739-Ele17225-Tri3-kH500-MismipPlus-io0027.mat';
goodfile_path='/work/ollie/orichter/MisomipPlus/io0027/fesommesh/meshgen.goodfile.1034.90';

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

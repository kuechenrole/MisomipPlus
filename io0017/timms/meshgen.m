meshOutPath='/work/ollie/orichter/MisomipPlus/io0017/fesommesh/1068.60/';
Ua_path='/work/ollie/orichter/MisomipPlus/io0017/uadata/1068.60-Nodes8709-Ele17169-Tri3-kH1000-MismipPlus-io0017.mat';
goodfile_path='/work/ollie/orichter/MisomipPlus/io0017/fesommesh/meshgen.goodfile.1068.60';

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

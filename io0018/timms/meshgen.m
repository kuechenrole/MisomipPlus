meshOutPath='/work/ollie/orichter/MisomipPlus/io0018/fesommesh/1051.80/';
Ua_path='/work/ollie/orichter/MisomipPlus/io0018/uadata/1051.80-Nodes9749-Ele19246-Tri3-kH1000-MismipPlus-io0018.mat';
goodfile_path='/work/ollie/orichter/MisomipPlus/io0018/fesommesh/meshgen.goodfile.1051.80';

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

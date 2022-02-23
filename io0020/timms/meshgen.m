meshOutPath='/work/ollie/orichter/MisomipPlus/io0020/fesommesh/1034.20/';
Ua_path='/work/ollie/orichter/MisomipPlus/io0020/uadata/1034.20-Nodes2027-Ele3844-Tri3-kH1000-MismipPlus-io0020.mat';
goodfile_path='/work/ollie/orichter/MisomipPlus/io0020/fesommesh/meshgen.goodfile.1034.20';

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

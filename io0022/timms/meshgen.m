meshOutPath='/work/ollie/orichter/MisomipPlus/io0022/fesommesh/1039.60/';
Ua_path='/work/ollie/orichter/MisomipPlus/io0022/uadata/1039.60-Nodes8944-Ele17636-Tri3-kH1000-MismipPlus-io0022.mat';
goodfile_path='/work/ollie/orichter/MisomipPlus/io0022/fesommesh/meshgen.goodfile.1039.60';

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

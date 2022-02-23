meshOutPath='/work/ollie/orichter/MisomipPlus/io0016/fesommesh/1090.60/';
Ua_path='/work/ollie/orichter/MisomipPlus/io0016/uadata/1090.60-Nodes9745-Ele19235-Tri3-kH1000-MismipPlus-io0016.mat';
goodfile_path='/work/ollie/orichter/MisomipPlus/io0016/fesommesh/meshgen.goodfile.1090.60';

cd meshgen;

disp('jigsaw2fesom');
jigsaw2fesomUa;

disp('makeTopoUa');
makeTopoUa;

%disp('ensure min overlap');
%ensure_minLayers_ovl;

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

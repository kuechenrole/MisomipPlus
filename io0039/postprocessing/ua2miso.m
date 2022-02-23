%execute with e.g.: matlab.sh -s -S"-wprod-0304" -M"-nojvm -r ua2miso('uadata_path','postprocessing_path/IceOcean1r_COM_ice_UaFesom_SSATsai.nc',1,100);exit"

function ua2miso(resultsFolder,ncfilename,start,stop)

filePattern = fullfile(resultsFolder, '*.mat'); % Change to whatever pattern you need.
theFiles = dir(filePattern);
theFiles = theFiles(start:stop);

nx = 800; dx=1e3;
ny = 80; dy=1e3;

x = linspace(dx/2,dx*(nx-1/2),nx);
y = linspace(dy/2,dy*(ny-1/2),ny);

if exist(ncfilename,'file')==2
    delete(ncfilename);
end

        
variables.name={'time','iceVolume','iceVAF','groundedArea','xGL','yGL',...
                'iceThicknessGL','uBaseGL','vBaseGL','uSurfaceGL','vSurfaceGL','uMeanGL','vMeanGL',...
                'x','y','iceThickness','upperSurface','lowerSurface','basalMassBalance',...
                'groundedMask','floatingMask','basalTractionMagnitude','uBase','vBase','uSurface','vSurface',...
                'uMean','vMean'};
dimensions.name={'nx','ny','nPointGL','nTime'};
dimensions.dim=[length(x);length(y);Inf;length(theFiles)];
variables.dim= [0 0 0 1; 0 0 0 1; 0 0 0 1; 0 0 0 1; 0 0 1 1; 0 0 1 1; ...
                0 0 1 1; 0 0 1 1; 0 0 1 1; 0 0 1 1; 0 0 1 1; 0 0 1 1; 0 0 1 1;...
                1 0 0 0; 0 1 0 0; 1 1 0 1; 1 1 0 1; 1 1 0 1; 1 1 0 1;...
                1 1 0 1; 1 1 0 1; 1 1 0 1; 1 1 0 1; 1 1 0 1; 1 1 0 1; 1 1 0 1;...
                1 1 0 1; 1 1 0 1];

variables.units={'s','m^3','m^3','m^2','m','m',...
                         'm','m/s','m/s','m/s','m/s','m/s','m/s',...
                         'm','m','m','m','m','m/s',...
                         'unitless','unitless','Pa','m/s','m/s','m/s','m/s',...
                         'm/s','m/s'};
%variables.description={'ice thickness','upper surface elevation','lower surfae elevation','basal mass balance of ice',...
%    'fraction of grounded ice in a given cell','fraction of floating ice in a given cell','magnitude of tangential basal traction field',...
%    'x component of basal velocity','y component of basal velocity',''};
     
for ii=1:length(variables.name)
    I = find(variables.dim(ii,:));
    dim = {dimensions.name{I(1)},dimensions.dim(I(1))};
    for jj=2:length(I)
        dim = {dim{:}, dimensions.name{I(jj)},dimensions.dim(I(jj))};
    end
    nccreate(ncfilename,variables.name{ii},...
         'Dimensions', dim,'FillValue', NaN,...
         'Format','netcdf4','Datatype','single');
    ncwriteatt(ncfilename,variables.name{ii},'units',variables.units{ii});
    %ncwriteatt(ncfilename,variables.name{ii},'description',variables.description{ii});
end
ncwriteatt(ncfilename,'/','_FillValue',NaN);

%read Ua data for all variables
for ii=1:length(variables.name)
    tstart = tic;
    data = readUaVariables(variables.name{ii},resultsFolder,ncfilename,start,stop);
    %'nx','ny','nPointGL','nTime'
    ncwrite(ncfilename,variables.name{ii},data);
    telapsed = toc(tstart);
    disp([variables.name{ii},', elapsed time: ',num2str(telapsed,3),'s']);
end
end


function data = readUaVariables(varName,resultsFolder,ncfilename,start,stop)

filePattern = fullfile(resultsFolder, '*.mat'); % Change to whatever pattern you need.
theFiles = dir(filePattern);
theFiles = theFiles(start:stop);

files_subset=[];
for k=1:length(theFiles)
    baseFileName = theFiles(k).name;
    files_subset{k} = fullfile(theFiles(k).folder, baseFileName);
end


nx = 800; dx=1e3;
ny = 80; dy=1e3;

x = linspace(dx/2,dx*(nx-1/2),nx);
y = linspace(dy/2,dy*(ny-1/2),ny);

data=[];

switch varName
    case 'time' % in seconds, sampled monthly
        %data=[startmonth:endmonth]*365/12*24*60*60;
        for jj=1:length(files_subset)
            fileNameParts= strsplit(theFiles(jj).name,'-');
            timeStr = strsplit(fileNameParts{1},'.');
            year = str2num(timeStr{1});
            month = str2num(timeStr{2});
            dim = [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
            data(jj)=((year-1000)*365 + sum(dim(1:month+1)))*24*3600;
        end
        
    case 'iceVolume' % in m3
        for jj=1:length(files_subset)
            load(string(files_subset(jj)));
            [~,iceVolume,~]=CalcVAF(CtrlVar,MUA,F.h,F.B,F.S,F.rho,F.rhow,F.GF);
            data(jj) = iceVolume.Total;
        end
        
    case 'iceVAF' % in m3
        for jj=1:length(files_subset)
            load(string(files_subset(jj)));
            [iceVAF,~,~]=CalcVAF(CtrlVar,MUA,F.h,F.B,F.S,F.rho,F.rhow,F.GF);
            data(jj) = iceVAF.Total;
        end
        
    case 'groundedArea' % in m2
        for jj=1:length(files_subset)
            load(string(files_subset(jj)));
            [~,~,groundedArea] = CalcVAF(CtrlVar,MUA,F.h,F.B,F.S,F.rho,F.rhow,F.GF);
            data(jj) = groundedArea.Total;
        end
        
    case 'xGL'
        for jj=1:length(files_subset)
            load(string(files_subset(jj)));
            GLgeo=GLgeometry(MUA.connectivity,MUA.coordinates,F.GF,CtrlVar);
            [xGL,~] = ArrangeGroundingLinePos(CtrlVar,GLgeo,1);
            [m,n] = size(data);
            I = max([m,length(xGL)]);
            datanew = NaN*ones(I,jj);
            datanew(1:m,1:n) = data;
            datanew(1:length(xGL),jj) = xGL(:);
            data = datanew;
        end
        
     case 'yGL'
        for jj=1:length(files_subset)
            load(string(files_subset(jj)));
            GLgeo=GLgeometry(MUA.connectivity,MUA.coordinates,F.GF,CtrlVar);
            [~,yGL] = ArrangeGroundingLinePos(CtrlVar,GLgeo,1);
            [m,n] = size(data);
            I = max([m,length(yGL)]);
            datanew = NaN*ones(I,jj);
            datanew(1:m,1:n) = data;
            datanew(1:length(yGL),jj) = yGL(:);
            data = datanew;
        end   
        
    case 'iceThicknessGL'
        for jj=1:length(files_subset)
            load(string(files_subset(jj)));
            GLgeo=GLgeometry(MUA.connectivity,MUA.coordinates,F.GF,CtrlVar);
            [xGL,yGL] = ArrangeGroundingLinePos(CtrlVar,GLgeo,1);
            Fh = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),F.h);
            hGL = Fh(xGL,yGL);
            [m,n] = size(data);
            I = max([m,length(xGL)]);
            datanew = NaN*ones(I,jj);
            datanew(1:m,1:n) = data;
            datanew(1:length(xGL),jj) = hGL(:);
            data = datanew;
        end
        
    case {'uBaseGL', 'uSurfaceGL', 'uMeanGL'} % in m/s
        for jj=1:length(files_subset)
            load(string(files_subset(jj)));
            GLgeo=GLgeometry(MUA.connectivity,MUA.coordinates,F.GF,CtrlVar);
            [xGL,yGL] = ArrangeGroundingLinePos(CtrlVar,GLgeo,1);
            Fub = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),F.ub);
            ubGL = Fub(xGL,yGL)/(365*24*60*60);
            [m,n] = size(data);
            I = max([m,length(xGL)]);
            datanew = NaN*ones(I,jj);
            datanew(1:m,1:n) = data;
            datanew(1:length(xGL),jj) = ubGL(:);
            data = datanew;
        end
        
    case {'vBaseGL', 'vSurfaceGL', 'vMeanGL'} % in m/s
        for jj=1:length(files_subset)
            load(string(files_subset(jj)));
            GLgeo=GLgeometry(MUA.connectivity,MUA.coordinates,F.GF,CtrlVar);
            [xGL,yGL] = ArrangeGroundingLinePos(CtrlVar,GLgeo,1);
            Fvb = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),F.vb);
            vbGL = Fvb(xGL,yGL)/(365*24*60*60);
            [m,n] = size(data);
            I = max([m,length(xGL)]);
            datanew = NaN*ones(I,jj);
            datanew(1:m,1:n) = data;
            datanew(1:length(xGL),jj) = vbGL(:);
            data = datanew;
        end
        
    case 'x'
        data=x;
        
    case 'y'
        data=y;
        
    case 'iceThickness'
        for jj=1:length(files_subset)
            load(string(files_subset(jj)));
            Fh = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),F.h);
            [X,Y]=ndgrid(x,y);
            h = Fh(X,Y); h(641:end,:)=NaN;
            [m,n,p] = size(data);
            I = max([m,length(x)]);
            J = max([n,length(y)]);
            datanew = NaN*ones(I,J,jj);
            datanew(1:m,1:n,1:p) = data;
            datanew(1:length(x),1:length(y),jj) = h;
            data = datanew;
        end
        
    case 'upperSurface'
        for jj=1:length(files_subset)
            load(string(files_subset(jj)));
            Fs = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),F.s);
            [X,Y]=ndgrid(x,y);
            s = Fs(X,Y); s(641:end,:)=NaN;
            [m,n,p] = size(data);
            I = max([m,length(x)]);
            J = max([n,length(y)]);
            datanew = NaN*ones(I,J,jj);
            datanew(1:m,1:n,1:p) = data;
            datanew(1:length(x),1:length(y),jj) = s;
            data = datanew;
        end
        
    case 'lowerSurface'
        for jj=1:length(files_subset)
            load(string(files_subset(jj)));
            Fb = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),F.b);
            [X,Y]=ndgrid(x,y);
            b = Fb(X,Y); b(641:end,:)=NaN;
            [m,n,p] = size(data);
            I = max([m,length(x)]);
            J = max([n,length(y)]);
            datanew = NaN*ones(I,J,jj);
            datanew(1:m,1:n,1:p) = data;
            datanew(1:length(x),1:length(y),jj) = b;
            data = datanew;
        end
        
    case 'basalMassBalance' % m/s
        for jj=1:length(files_subset)
            load(string(files_subset(jj)));
            Fab = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),F.ab);
            [X,Y]=ndgrid(x,y);
            ab = -Fab(X,Y)/(365.24*24*60*60); ab(641:end,:)=NaN;
            ab(ab==0)=NaN;
            [m,n,p] = size(data);
            I = max([m,length(x)]);
            J = max([n,length(y)]);
            datanew = NaN*ones(I,J,jj);
            datanew(1:m,1:n,1:p) = data;
            datanew(1:length(x),1:length(y),jj) = ab;
            data = datanew;
        end
        
    case 'groundedMask'
        for jj=1:length(files_subset)
            load(string(files_subset(jj)));
            Fgf = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),double(F.GF.node));
%            fun = @(x,y) Fgf(x,y);
            [X,Y]=ndgrid(x,y);
            % run over all cells and find those that have either 0s or 1s
            % at all corners
            gf_ll = Fgf(X(:)-dx/2,Y(:)-dy/2);
            gf_ul = Fgf(X(:)-dx/2,Y(:)+dy/2);
            gf_tr = Fgf(X(:)+dx/2,Y(:)+dy/2);
            gf_lr = Fgf(X(:)+dx/2,Y(:)-dy/2);
            gf = (gf_ll+gf_ul+gf_tr+gf_lr)/4;
%             % run integral over those cells that have gf>0 and <1
%             I = find(gf>0 & gf<1);
%             for kk=1:length(I)
%                     gf(I(kk)) = integral2(fun,X(I(kk))-dx/2,X(I(kk))+dx/2,Y(I(kk))-dy/2,Y(I(kk))+dy/2);
%             end
%             % check that for other cells, there are no nodes inside cell
%             % that have values <1 or >0
%             I = find(gf==0 | gf==1);
%             count=0;
%             for kk=1:length(I)
%                squareX = [X(I(kk))-dx/2,X(I(kk))-dx/2,X(I(kk))+dx/2,X(I(kk))+dx/2];
%                squareY = [Y(I(kk))-dy/2,Y(I(kk))+dy/2,Y(I(kk))+dy/2,Y(I(kk))-dy/2];
%                J = find(inpoly([MUA.coordinates(:,1),MUA.coordinates(:,2)],[squareX(:),squareY(:)]));
%                K = find(GF.node(J)~=0 | GF.node(J)~=1);
%                if ~isempty(K)
%                    gf(I(kk)) = integral2(fun,X(I(kk))-dx/2,X(I(kk))+dx/2,Y(I(kk))-dy/2,Y(I(kk))+dy/2);
%                    count = count+1;            
%                    disp(count);
%                end
%             end
            gf = reshape(gf,size(X)); gf(641:end,:)=NaN;
            [m,n,p] = size(data);
            I = max([m,length(x)]);
            J = max([n,length(y)]);
            datanew = NaN*ones(I,J,jj);
            datanew(1:m,1:n,1:p) = data;
            datanew(1:length(x),1:length(y),jj) = gf;
            data = datanew;
        end
        
    case 'floatingMask'
        data = 1-ncread(ncfilename,'groundedMask');
 
        
    case 'basalTractionMagnitude'
        for jj=1:length(files_subset)
            load(string(files_subset(jj)));
            [~,~,tb]=CalcBasalTraction(CtrlVar,0,MUA,F);
            Ftb = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),tb);
            [X,Y]=ndgrid(x,y);
            tb = Ftb(X,Y);  tb(641:end,:)=NaN;
            [m,n,p] = size(data);
            I = max([m,length(x)]);
            J = max([n,length(y)]);
            datanew = NaN*ones(I,J,jj);
            datanew(1:m,1:n,1:p) = data;
            datanew(1:length(x),1:length(y),jj) = tb;
            data = datanew;
        end
        
    case {'uBase', 'uSurface', 'uMean'}
        for jj=1:length(files_subset)
            load(string(files_subset(jj)));
            Fub = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),F.ub);
            [X,Y]=ndgrid(x,y);
            ub = Fub(X,Y)/(365*24*60*60);  ub(641:end,:)=NaN;
            [m,n,p] = size(data);
            I = max([m,length(x)]);
            J = max([n,length(y)]);
            datanew = NaN*ones(I,J,jj);
            datanew(1:m,1:n,1:p) = data;
            datanew(1:length(x),1:length(y),jj) = ub;
            data = datanew;
        end
        
    case {'vBase', 'vSurface', 'vMean'}
        for jj=1:length(files_subset)
            load(string(files_subset(jj)));
            Fvb = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),F.vb);
            [X,Y]=ndgrid(x,y);
            vb = Fvb(X,Y)/(365*24*60*60); vb(641:end,:)=NaN;
            [m,n,p] = size(data);
            I = max([m,length(x)]);
            J = max([n,length(y)]);
            datanew = NaN*ones(I,J,jj);
            datanew(1:m,1:n,1:p) = data;
            datanew(1:length(x),1:length(y),jj) = vb;
            data = datanew;
        end
end

end

%'iceVAF','groundedArea','xGL','yGL',...
%   'iceThicknessGL','uBaseGL','vBaseGL','uSurfaceGL','vSurfaceGL','uMeanGL','vMeanGL',...
%   'x','y','iceThickness','upperSurface','lowerSurface','basalMassBalance',...
%   'groundedMask','floatingMask','basalTractionMagnitude','uBase','vBase','uSurface','vSurface',...
%   'uMean','vMean'

function [ball] = pwrbal1(pp,pw,ee)
%PWRBAL1 compute the ortho-balls associated with a 1-simplex
%triangulation embedded in R^2.
%   [BB] = PWRBAL1(PP,PW,TT) returns the set of power balls
%   associated with the edges in [PP,PW,EE], such that BB = 
%   [XC,YC,RC.^2]. PW is a vector of vertex weights.

%   Darren Engwirda : 2017 --
%   Email           : de2363@columbia.edu
%   Last updated    : 04/12/2017

%---------------------------------------------- basic checks    
    if ( ~isnumeric(pp) || ...
         ~isnumeric(pw) || ...
         ~isnumeric(ee) )
        error('pwrbal1:incorrectInputClass' , ...
            'Incorrect input class.');
    end
    
%---------------------------------------------- basic checks
    if (ndims(pp) ~= +2 || ...
        ndims(pw) ~= +2 || ...
        ndims(ee) ~= +2 )
        error('pwrbal1:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end
    
    if (size(pp,2)~= +2 || ...
            size(pp,1)~= size(pw,1) || ...
                size(ee,2) < +2 )
        error('pwrbal1:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end
    
%------------------------------------------------ lin offset
    xx12 = pp(ee(:,1),:) - ...
           pp(ee(:,2),:) ;
       
    ww12 = pw(ee(:,1),1) - ...
           pw(ee(:,2),1) ;
       
    dp12 = sum(xx12.*xx12,2) ;
    
    tpwr = +.5 * (ww12+dp12)./dp12 ;

    ball = zeros(size(ee,1),3) ;
    ball(:,1:2) = ...
        pp(ee(:,1),:) - tpwr.*xx12 ;
    
    vsq1 = ...
        pp(ee(:,1),:) - ball(:,1:2);
    vsq2 = ...
        pp(ee(:,2),:) - ball(:,1:2);
 
%------------------------------------------------ mean radii
    rsq1 = sum(vsq1 .^ 2,2) ;
    rsq2 = sum(vsq2 .^ 2,2) ;
    
    rsq1 = rsq1-pw(ee(:,1)) ;
    rsq2 = rsq2-pw(ee(:,2)) ;
    
    ball(:,3) = (rsq1 + rsq2) / 2. ;
    
end




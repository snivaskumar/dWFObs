function [xe,Ce] = dfuze(z,Z,x,hr,n,Zkk1,xkk1,x_est,typeCZ,typeWeight);
% function [xe,Ce] = fuze(z,Z,x,hr,n,x_est,typeCZ);
% [xe,Ce] = fuze(z,Z,x,hr,n,x_est,typeCZ);
% typeCZ  1 if Z = Co-Variance,  
%       2 if Z = Information,

H = [];
k = 0;
nn = 0;
xfus = x{1};
for i = 2:hr
    xfus = union( xfus,x{i} );
end
n = length(xfus);

for i = 1:hr
    [x12,x1,x2]     = intersect( xfus,x{i} );
    x12             = sort(x12);
    l               = length(x{i});
    nn              = nn + l;
    H(k+1:nn,:)     = zeros(l,n);
    H(k+1:nn,x1)  = eye(l,l);
    k               = k + l;
end

Ze  = [];
nn  = 0;
k   = 0;  

if strcmp(typeCZ,'C')
    C = Z;
elseif strcmp(typeCZ,'Z')
    for i = 1:hr
        C{i} = inv(Z{i});
    end
end

Zp = zeros(n,n);
for i = 1:hr
    [x12,x1,x2]         = intersect( xfus,x{i} );
    x12                 = sort(x12);
    l                   = length(x{i});
    nn                  = nn + l;
    
    if strcmp(typeWeight,'CONSTANT')
%         w = l;
        w = 1/hr;
        Ze(k+1:nn,k+1:nn)   = pinv( (1/w).*C{i} );
    elseif strcmp(typeWeight,'OPTIMAL')
        Ze(k+1:nn,k+1:nn)   = C{i};
%         f = @(w) trace((1/w).*C{i});
%         w = fminbnd(f,0,1,optimset('Display','off'));
%         Ze(k+1:nn,k+1:nn)   = inv( (1/w).*C{i} );
    end
    if strcmp(typeCZ,'C')
        xx(k+1:nn,:)        = z{i};
    elseif strcmp(typeCZ,'Z')
        xx(k+1:nn,:)        = C{i}*z{i};
    end

    k                   = k + l;
end
if strcmp(typeWeight,'OPTIMAL')
    f   = @(w) trace((1/w).*Ze);
    w   = fminbnd(f,0,1,optimset('Display','off'));
    Ze  = inv( (1/w).*Ze );
end

H = sparse(H);

% size(Zp)
% size(H'*Ze*H)
Ce = pinv(H'*Ze*H);

% size(Ze)
% size(H)
% size(xx)
% size(Ce)
xe = Ce*H'*Ze*xx;
end
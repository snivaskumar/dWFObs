% % % function [zf, Zf] = d_fuze2(z1,z2,Z1,Z2,x1,x2,type);
% % % xlkk = cell(2,1);
% % % Slkk = cell(2,1);
% % % Slkk{1} = inv(Z1);
% % % Slkk{2} = inv(Z2);
% % % xlkk{1} = Z1*z1;
% % % xlkk{2} = Z2*z2;
% % % 
% % % x = cell(2,1);
% % % x{1} = x1;
% % % x{2} = x2;
% % % z = cell(2,1);
% % % Z = cell(2,1);
% % % for i = 1:2
% % %     for j = 1:2
% % %         if i~=j
% % %             a = xlkk{i};
% % %             b = xlkk{j};
% % %             P1 = Slkk{i};
% % %             P2 = Slkk{j};
% % % 
% % %             l1 = length(x{i});
% % %             [cc1,ia1]=setdiff(x{i},x{j});
% % %             Tmp1 = zeros(l1,l1);
% % %             Tmp1(ia1,ia1) = eye(length(ia1),length(ia1));
% % %             Tmp1 = Tmp1(ia1,[1:l1]);
% % % 
% % %             l2 = length(x{j});
% % %             [ha,cc2,ia2] = intersect(x{i},x{j});
% % %             Tmp2 = zeros(l2,l2);
% % %             Tmp2(ia2,ia2) = eye(length(ia2),length(ia2));
% % %             Tmp2 = Tmp2(ia2,[1:l2]);
% % % 
% % %             xa = x{i};
% % %             xb = [Tmp2*x{j}; Tmp1*x{i}];
% % % 
% % %             [xb,ix] = sort(xb);
% % %             xb;
% % % 
% % %             a       = a;
% % %             btmp    = [Tmp2*b; Tmp1*a];
% % %             b       = btmp(ix);
% % %             xb;
% % %             
% % %             lc1 = length(cc1);
% % %             lc2 = length(cc2);
% % % 
% % %             PP = zeros(l1,l1);
% % %             PP(cc2,cc2) = Tmp2*P2*Tmp2';
% % % 
% % %             [Sa,Da] = eig(P1);
% % %             [Sb,Db] = eig(PP);
% % % 
% % %             [~,I] = sort(diag(Db),'descend');
% % %             Sb = Sb(I,I);
% % %             Db = Db(I,I);
% % % 
% % %             for ii = 1:length(x{i})
% % %                 if Db(ii,ii) == 0
% % %                     Gb(ii,ii) = 1;
% % %                 else
% % %                     Gb(ii,ii) = Db(ii,ii);
% % %                 end
% % %             end
% % %             Xb = Sa*(Da^0.5)*Sb*Gb*pinv(Sb)*(Da^0.5)*pinv(Sa);
% % %             Xa = P1;
% % % 
% % %             x_a{i} = a;
% % %             x_b{i} = b;
% % %             X_a{i} = Xa;
% % %             X_b{i} = Xb;
% % %         end
% % %     end
% % % end

% function [zf, Zf] = d_fuze2(z1,z2,Z1,Z2,x1,x2,type);
% % [zf{i}, Zf{i}, xf] = fuze2(z1,z2,Z1,Z2,x1,x2,type);
% % type    1 for CI, 2 for EI, 3 for ICI
% 
% % x1 = [1;2;3;4];
% % x2 = [2;4];
% % z1 = rand(4,1);
% % Z1 = rand(4,4);
% % z2 = rand(2,1);
% % Z2 = rand(2,2);
% 
% 
% l1 = length(x1);
% l2 = length(x2);
% 
% Z1 = inv(Z1);
% Z2 = inv(Z2);
% z1 = Z1*z1;
% z2 = Z2*z2;
% 
% x1a = z1;
% P1a = Z1;
% x1b = zeros(l1,1);
% P1b = zeros(l1,l1);
% tmpx1 = [];
% tmpx2 = [];
% for i = 1:l1
%     if sum(x1(i) == x2)
%         tmp = [1:l1]';
%         tmp = tmp( x1(i) == x2 );
%         tmpx2 = [tmpx2;tmp];
%         tmpx1 = [tmpx1;i];
%     end
% end
% x1b(tmpx1)       = z2(tmpx2);
% P1b(tmpx1,tmpx1)  = Z2(tmpx2,tmpx2);
% 
% x2b = z2;
% P2b = Z2;
% x2a = zeros(l2,1);
% P2a = zeros(l2,l2);
% tmpx1 = [];
% tmpx2 = [];
% for i = 1:l2
%     if sum(x2(i) == x1)
%         tmp = [1:l1]';
%         tmp = tmp( x2(i) == x1 );
%         tmpx2 = [tmpx2;tmp];
%         tmpx1 = [tmpx1;i];
%     end
% end
% x2a(tmpx1)       = z1(tmpx2);
% P2a(tmpx1,tmpx1)  = Z1(tmpx2,tmpx2);
% 
% % z1'
% % z2'
% % x1a'
% % x1b'
% % 
% % z1'
% % z2'
% % x2a'
% % x2b'
% 
% x_a{1} = x1a;
% x_b{1} = x1b;
% X_a{1} = P1a;
% X_b{1} = P1b;
% x_a{2} = x2a;
% x_b{2} = x2b;
% X_a{2} = P2a;
% X_b{2} = P2b;
% for i = 1:2
%     if strcmp(type,'CIN')
%         Zf{i}      = X_a{i} + X_b{i};
%         zf{i}      = x_a{i} + x_b{i};
%         0;
%     elseif strcmp(type,'CI')
%         % CI
%         ZA      = X_a{i};
%         ZB      = X_b{i};
%         f       = @(w) trace( pinv(w*ZA + (1-w)*ZB) ); % arg (min -f) = arg (max f)
%     omega   = fminbnd(f,0,1,optimset('Display','off'));
% %         omega   = 0.2;
%         Zf{i}      = omega*ZA + (1-omega)*ZB;
%         zf{i}      = omega*x_a{i} + (1-omega)*x_b{i};
%         1;
%     elseif strcmp(type,'EI')
%         % EI
%         [Si,Di] = eig(X_a{i});
%         [Sj,Dj] = eig(pinv(Di^0.5)*pinv(Si)*X_b{i}*Si*pinv(Di^0.5));
% 
%         Dij = zeros(size(Dj));
%         for ii = 1: length(Dij)
%             Dij(ii,ii) = min(1,Dj(ii,ii));
%         end
%         Xij = Si*(Di^0.5)*Sj*Dij*pinv(Sj)*(Di^0.5)*pinv(Si);
% 
%     %         w1 = (X_b{i} - 0.5*Xij)*pinv(X_a{i} + X_b{i} - Xij);
%     %         w2 = (X_a{i} - 0.5*Xij)*pinv(X_a{i} + X_b{i} - Xij);
% 
%         w1 = (X_b{i})*pinv(X_a{i} + X_b{i});
%         w2 = (X_a{i})*pinv(X_a{i} + X_b{i});
%         xij = w1*x_a{i} + w2*x_b{i};
% 
%         zf{i}  = x_a{i} + x_b{i} - xij;
%         Zf{i}  = X_a{i} + X_b{i} - Xij;
%         2;
%     elseif strcmp(type,'ICI')
%         % ICI
%         ZA      = X_a{i};
%         ZB      = X_b{i};
% % %         ff      = @(w) trace(pinv(ZA + ZB - ZA*pinv(w*ZB + (1-w)*ZA)*ZB)); % min -f = max f
%         ff      = @(w) trace(-(ZA + ZB - ZA*inv(w*ZB + (1-w)*ZA)*ZB)); % min -f = max f
%         omega   = fminbnd(ff,0,1,optimset('Display','off'));
% %         omega   = 0.5;
%         Xij     = ZA*inv(omega*ZB + (1-omega)*ZA)*ZB;
%         Zf{i}      = X_a{i} + X_b{i} - Xij;
% 
%         K = ( X_a{i} - (omega)*Xij )*pinv(Zf{i});
%         L = ( X_b{i} - (1 - omega)*Xij )*pinv(Zf{i});
%         zf{i} = K*x_a{i} + L*x_b{i};
%         3;
%     end   
%     Zf{i} = inv(Zf{i});
%     zf{i} = Zf{i}*zf{i};
% end

function [zf, Zf, xf] = d_fuze2(z1,z2,Z1,Z2,x1,x2,type,fusion_weight,constant,iteration);
% [zf, Zf, xf] = fuze2(z1,z2,Z1,Z2,x1,x2,type);
% type    1 for CI, 2 for EI, 3 for ICI

Z1 = pinv(Z1);
Z2 = pinv(Z2);
z1 = Z1*z1;
z2 = Z2*z2;

x_tmp                           = union( x1,x2 );
[cc{1},ia{1}]                   = setdiff( x_tmp,x1 );
T_tmp{1}                        = zeros(max(x_tmp),max(x_tmp));
T_tmp{1}(x_tmp,x_tmp)           = eye(length(x_tmp),length(x_tmp));
T_tmp{1}(cc{1},cc{1})           = zeros(length(cc{1}),length(cc{1}));
% T_tmp{i}                        = T_tmp{i}(x_tmp,setdiff(x_tmp,ia{i}));
T_tmp{1}                        = T_tmp{1}(x_tmp,x1);
T_tmp{1}                        = sparse(T_tmp{1});
x_a                             = T_tmp{1}*z1;
x1 = T_tmp{1}*x1;

% clear cc ia
[cc{2},ia{2}]                   = setdiff( x_tmp,x2 );
T_tmp{2}                        = zeros(max(x_tmp),max(x_tmp));
T_tmp{2}(x_tmp,x_tmp)           = eye(length(x_tmp),length(x_tmp));
T_tmp{2}(cc{2},cc{2})           = zeros(length(cc{2}),length(cc{2}));  
% T_tmp{i+1}                      = T_tmp{i+1}(x_tmp,setdiff(x_tmp,ia{i+1}));
T_tmp{2}                        = T_tmp{2}(x_tmp,x2);
T_tmp{2}                        = sparse(T_tmp{2});
x_b                             = T_tmp{2}*z2;
x2 = T_tmp{2}*x2;

X_a = T_tmp{1}*Z1*T_tmp{1}';
X_b = T_tmp{2}*Z2*T_tmp{2}';

% X_a = pinv(X_a);
% X_b = pinv(X_b);
% x_a = X_a*x_a;
% x_b = X_b*x_b;

%%%%%%
% % Extra
% x_a(ia{1}) = x_b(ia{1});
% x_b(ia{2}) = x_a(ia{2});
% 
% % x1(ia{1}) = x2(ia{1});
% % x2(ia{2}) = x1(ia{2});
% % x1;
% % x2;
% 
% X_a(ia{1},ia{1}) = X_b(ia{1},ia{1});
% X_b(ia{2},ia{2}) = X_a(ia{2},ia{2});
%%%%%%
% x1a = x_a;
% X1a = X_a;
% x1b = x_b;
% X1b = X_b;
% x2a = x_a;
% X2a = X_a;
% x2b = x_b;
% X2b = X_b;
x1a = x_a(x1~=0);
X1a = X_a(x1~=0,x1~=0);
x1b = x_b(x1~=0);
X1b = X_b(x1~=0,x1~=0);
x2a = x_a(x2~=0);
X2a = X_a(x2~=0,x2~=0);
x2b = x_b(x2~=0);
X2b = X_b(x2~=0,x2~=0);
x_a = cell(2,1);
x_b = cell(2,1);
X_a = cell(2,1);
X_b = cell(2,1);
x_a{1} = x1a;
x_b{1} = x1b;
X_a{1} = X1a;
X_b{1} = X1b;
x_a{2} = x2a;
x_b{2} = x2b;
X_a{2} = X2a;
X_b{2} = X2b;
for i = 1:2
    if strcmp(type,'CIN')
        Zf{i}      = X_a{i} + X_b{i};
        zf{i}      = x_a{i} + x_b{i};
        0;
    elseif strcmp(type,'CI')||strcmp(type,'CI2')
        % CI
        ZA      = X_a{i};
        ZB      = X_b{i};
        f       = @(w) trace( pinv(w*ZA + (1-w)*ZB) ); % arg (min -f) = arg (max f)
        if strcmp(fusion_weight,'OPTIMAL')
            options = optimset('MaxIter',iteration,'Display','off');
            omega   = fminbnd(f,0,1,options);
        else
            omega   = constant;
        end
        Zf{i}      = omega*ZA + (1-omega)*ZB;
        zf{i}      = omega*x_a{i} + (1-omega)*x_b{i};
        1;
    elseif strcmp(type,'EI')
        % EI
        [Si,Di] = eig(X_a{i});
        [Sj,Dj] = eig(pinv(Di^0.5)*pinv(Si)*X_b{i}*Si*pinv(Di^0.5));

        Dij = zeros(size(Dj));
        for ii = 1: length(Dij)
            Dij(ii,ii) = min(1,Dj(ii,ii));
        end
        Xij = Si*(Di^0.5)*Sj*Dij*pinv(Sj)*(Di^0.5)*pinv(Si);

    %         w1 = (X_b{i} - 0.5*Xij)*pinv(X_a{i} + X_b{i} - Xij);
    %         w2 = (X_a{i} - 0.5*Xij)*pinv(X_a{i} + X_b{i} - Xij);

        w1 = (X_b{i})*pinv(X_a{i} + X_b{i});
        w2 = (X_a{i})*pinv(X_a{i} + X_b{i});
        xij = w1*x_a{i} + w2*x_b{i};

        zf{i}  = x_a{i} + x_b{i} - xij;
        Zf{i}  = X_a{i} + X_b{i} - Xij;
        2;
    elseif strcmp(type,'ICI')
        % ICI
        ZA      = X_a{i};
        ZB      = X_b{i};
        ff      = @(w) trace( pinv(ZA + ZB - ZA*pinv(w*ZB + (1-w)*ZA)*ZB) ); % min -f = max f
%         ff      = @(w) trace(-(ZA + ZB - ZA*inv(w*ZB + (1-w)*ZA)*ZB)); % min -f = max f
        if strcmp(fusion_weight,'OPTIMAL')
            options = optimset('MaxIter',iteration,'Display','off');
            omega   = fminbnd(ff,0,1,options);
        else
            omega   = constant;
        end
        Xij     = ZA*pinv(omega*ZB + (1-omega)*ZA)*ZB;
        Zf{i}      = X_a{i} + X_b{i} - Xij;

        K = ( X_a{i} - (omega)*Xij )*pinv(Zf{i});
        L = ( X_b{i} - (1 - omega)*Xij )*pinv(Zf{i});
        zf{i} = K*x_a{i} + L*x_b{i};
        3;
    end   
    Zf{i} = pinv(Zf{i});
    zf{i} = Zf{i}*zf{i};
end
% xx = zf{1};
% zf{1} = xx(x1~=0);
% zf{2} = xx(x2~=0);
% pp = Zf{1};
% Zf{1} = pp(x1~=0,x1~=0);
% Zf{2} = pp(x2~=0,x2~=0);

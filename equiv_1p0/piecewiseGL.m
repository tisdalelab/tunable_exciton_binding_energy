function curve = piecewiseGL(wL,wG,Ec,choose,x)
lz = @(wL,Ec,x) (wL.^2)./((x-Ec).^2 + wL.^2);    %generic lorentzian
gss = @(wG,Ec,x) exp(-((x-Ec)./wG).^2);   %generic gaussian
p = zeros(size(x,1),1);     % input matrix x 
for k = 1:size(x,1)
    if choose == -1
            if x(k) > Ec
                p(k) = lz(wL,Ec,x(k));
            elseif x(k) <= Ec
                p(k) = gss(wG,Ec,x(k));
            end
                          
   elseif choose == 1
            if x(k) < Ec
                p(k) = gss(wG,Ec,x(k));
            elseif x(k) >= Ec
                p(k) = lz(wL,Ec,x(k));
                
            end

    end
end

%curve = p';

% output p is 401x1
if size(x,1) == 1 % 401, i.e. number of rows
    curve = p';
elseif size(x,2) == 1 % 1, i.e. number of columns
    curve = p;
end

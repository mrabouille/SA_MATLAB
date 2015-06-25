function p_n = legendrecoef(n)
%LEGENDRECOEF compute the Legendre polynomial coefficient
%return a n+1-by-n+1 matrix
%
% The result p_n is a i-by-j matrix, each element is the {j-1}-order 
% coefficient of P_{i} polynome
%
% P_{i} = p_n(i,1)+ p_n(i,2)*X + ... + p_n(i,j)*X^(j-1) + ...


morm = true;
sparceResult = true;

if sparceResult
    p_n = sparse([],[],[],n+1,n+1,ceil((n+1)/2)*(n+2)/2);
else
    p_n = zeros(n+1,n+1);
end

if n>=0
    p_n(1,1) = 1;
end

if n>=1
    p_n(2,2) = 1;
end

if n>1

    for k=3:n+1     
        for a=k:-2:2
            p_n(k,a) = (2*k-3)*p_n(k-1,a-1) - (k-2)*p_n(k-2,a);
        end        
        p_n(k,1) = p_n(k,1) - (k-2)*p_n(k-2,1);
        p_n(k,:) = p_n(k,:)/(k-1);
    end

    if morm
        for k=3:n+1
            p_n(k,:) = p_n(k,:)*sqrt(2*k-1);
        end
    end
end

% full(p_n)

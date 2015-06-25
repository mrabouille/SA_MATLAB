function h_n = hermitecoef(n,choise)
%HERMITECOEF compute the Hermite polynomial coefficient 
%return a n+1-by-n+1 matrix
%
% The result h_n is a i-by-j matrix, each element is the {j-1}-order
% coefficient of H_{i} polynome
%
% H_{i} = h_n(i,1)+ h_n(i,2)*X + ... + h_n(i,j)*X^(j-1) + ...


morm = true;
sparceResult = true;

if nargin==2 && ~choise
    % for the 'physicists' Hermite polynomials
    c=2;
else
    % for 'probabilists' Hermite polynomials
    c=1;
end

if sparceResult
    h_n = sparse([],[],[],n+1,n+1,ceil((n+1)/2)*(n+2)/2);
else
    h_n = zeros(n+1,n+1);
end
    
if n>=0
    h_n(1,1) = 1;
end

if n>=1
    h_n(2,2) = c;
end

if n>1

    for k=3:n+1
        for a=k-1:-2:1
            h_n(k,a+1) = c*( h_n(k-1,a) - (k-2)*h_n(k-2,a+1) );
        end
        h_n(k,1) = -c*(k-2)*h_n(k-2,1);
    end

    if morm
        for k=3:n+1
            h_n(k,:) = h_n(k,:)/sqrt(factorial(k-1));
        end
    end
%     for k=3:n+1
%         for a=k-1:-2:1
%             h_n(k,a+1) = c*( h_n(k-1,a) - (k-2)*h_n(k-2,a+1)/sqrt(k-2) )/sqrt(k-1);
%         end
%         h_n(k,1) = -c*(k-2)*h_n(k-2,1)/(sqrt(k-1)*sqrt(k-2));
%     end

end

% full(h_n)

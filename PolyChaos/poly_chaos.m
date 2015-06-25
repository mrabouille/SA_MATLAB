function [model,effets,indices,facteurs,famille]=poly_chaos(input,output,ordre,type)
%poly_chaos Calcule une aproximation polynomiale du modèle évalué
%
%   [model,indices,facteurs,famille]=poly_chaos(input,output,ordre,type)
%
%       input       : Entrées du modèle
%       output      : Sortie du modèle
%       ordre       : Ordre maximum admissible pour le polynôme
%       type        : Type de polynome à utiliser
%
%       model       : Méta-modèle trouvé
%       effets      : Represente les effets du 1er et second ordre
%       indices     : Indices des polynomes composant la famille créée
%       facteurs    : Facteurs de la famille
%       famille     : Famille de polynome créée
%  
% 


%% test
if nargin==0
    if false
        
    % ISHIGAMI fonction
        % Crestaux et al. (2007) and Marrel et al. (2009) use: a = 7 and b = 0.1.
        % Sobol' & Levitan (1999) use: a = 7 and b = 0.05. 
        a = 7;	b = 0.05;
        % def: X=[x1,x2,x3] -> xi ~ U(-pi,pi)
        ninput = 3;
        
        fonc = @(X) sin(X(:,1)) + a*(sin(X(:,2))).^2 + b*(X(:,3)).^4.*sin(X(:,1))

        E = a/2; % experance    
        Vx1 = 1/2*(1+b*pi^4/5)^2;
        Vx2 = a^2/8;
        Vx13 = b^2*pi^8*8/225;
        V = Vx1 + Vx2 + Vx13;

        exact = [ Vx1/V,    0,      0;
                  0,        Vx2/V,  0;
                  Vx13/V,   0,      0]
              
        rng shuffle
        input = -pi + 2*pi.*rand(5000,ninput);
        output = fonc(input);
   
        test_input = -pi + 2*pi.*rand(5000,ninput);
        test_output = fonc(test_input);
        
        ordre = 8;
        type = 2;
        
    else
    % SOBOL fonction
      %  a=[1 2 5 10 33 99];
        a=[5 10 20];
        ninput = length(a);
        fonc = '@(X) 1';
        for i=1:ninput
            fonc = [fonc, sprintf('.*( abs(4*X(:,%1$d)-2) + %2$d )/( 1+ %2$d )',i,a(i) )];
        end
        fonc = eval(fonc);
 
        
        Vx = 1./(3*(1+a).^2);
        V = prod(Vx+1)-1;
        
        exact = Vx/V;     
        
        rng shuffle
        input = rand(5000,ninput);
        output = fonc(input);
        
        ordre = 8;
        type = 2;

        test_input = rand(5000,ninput);
        test_output = fonc(test_input);
    end
    
end   

numVar=size(input,2);


%% initialisation
switch type
    case {1, 'He'}
        % Hermite: gaussienne sur R
        % probabilists polynomials ! Normé !
        base = {'1'
                'X'       
                '(X^2-1)/sqrt(2)'
                '(X^3-3*X)/sqrt(6)'
                '(X^4-6*X^2+3)/sqrt(24)'
                '(X^5-10*X^3+15*X)/sqrt(120)'
                '(X^6-15*X^4+45*X^2-15)/sqrt(720)'
                '(X^7-21*X^5+105*X^3-105*X)/sqrt(5040)'
                '(X^8-28*X^6+210*X^4-420*X^2+105)/sqrt(40320)'
                '(X^9-36*X^7+378*X^5-1260*X^3+945*X)/sqrt(362880)'
                '(X^10-45*X^8+630*X^6-3150*X^4+4725*X^2-945)/sqrt(3628800)'};
            
     %   full(HermitePoly(n))
        
   %     maxi = max(max(input)); mini = min( min(input));
        polyInput = input;

  

    case {2, 'P'}
        % Legendre: uniforme sur [-1;1]
        %{
        'P(0,x) = 1'
        'P(1,x) = x'
        'P(n,x)=(2*n-1)/n * x * P(n-1,x) - (n-1)/n * P(n-2,x)'
        %}
        base = {'1'
                'X*sqrt(2*1+1)'
                '(3*X^2-1)/2*sqrt(2*2+1)'
                '(5*X^3-3*X)/2*sqrt(2*3+1)'
                '(35*X^4-30*X^2+3)/8*sqrt(2*4+1)'
                '(63*X^5-70*X^3+15*X)/8*sqrt(2*5+1)'
                '(231*X^6-315*X^4+105*X^2-5)/16*sqrt(2*6+1)'
                '(429*X^7-693*X^5+315*X^3-35*X)/16*sqrt(2*7+1)'
                '(6435*X^8-12012*X^6+6930*X^4-1260*X^2+35)/128*sqrt(2*8+1)'
                '(12155*X^9-25740*X^7+18018*X^5-4620*X^3+315*X)/128*sqrt(2*9+1)'
                '(46189*X^10-109395*X^8+90090*X^6-30030*X^4+3465*X^2-63)/256*sqrt(2*10+1)'};
       
     %   full(LegendrePoly(n))
        
        maxi = max(max(input)); mini = min( min(input));
        polyInput = (input-mini)/(maxi-mini)*2-1;
        
        test_input = (test_input-mini)/(maxi-mini)*2-1;
 
end

if ordre>(size(base,1)-1)
    error('Le degrés du polynôme est trop élevé pour la famille choisie')
end


%% Constrction de la famille

% création des indices
[indices,numpoly]=indices_poly(ordre,numVar);

% nétoyage (pas d'interactions d'ordre supérieur à 3)
% indices(sum(indices>0,2)>2,:)=[];
% numpoly = size(indices,1);

% vérification de la résolubilité du système
% ==== A FAIRE ====
%size(indices,1)

% génération de la famille
A=base(indices+1);


for i=1:numVar
    A(:,i)=regexprep(A(:,i), '[xX]', ['X(' num2str(i) ')']);
end

%symbolic expression -> permet de faire des operations 
syms X
famille(X) = prod(sym(A),2);


%produit des warning >> A CORRIGER !
warning off
    equa = matlabFunction(famille);  %function -> permet d'être evaluée
warning on

[~, msgid] = lastwarn;
if ~strncmp(msgid, 'symbolic:mupadmex:MuPADTextWarning', length(msgid))
    warning(lastwarn)
end


%% Résolution
% recherche de la solution par projection

% evaluation de chaque polynome de la base sur les points d'evaluation
P=zeros(size(polyInput,1),numpoly);
for i=1:size(polyInput,1)
   P(i,:)= equa(polyInput(i,:));
end


if false
    % recherche de la solution par inversion
    facteurs=P\output;

else
    % recherche de la solution par projection

    % Initialisation avec coeffs=0 et Residu = Y-P*coeff = Y
    facteurs = zeros(numpoly,1);
    R = output;

    % Mesure de la moyenne et definition du 1er facteur
    Moy = mean(R);  R=R-Moy;  facteurs(1)=Moy;

    %Vari = var(R);  R = R/sqrt(Vari);
    
% Bradley Efron, Trevor Hastie, Iain Johnstone, and Robert Tibshirani. 2004.
% Least Angle Regression. The Annals of Statistics, 32:2, 407–99. doi:<10.1214/009053604000000067>
    [beta, A, mu, C, c, gamma] = lars(P(:,2:end), R, 'lar', Inf, false);

    facteurs(2:end) = beta(end,:)';
end


effets=zeros(numVar);
tot=sum(facteurs(2:end).^2);
for i=1:size(facteurs,1)
    elements = find(indices(i,:)>0);
    if isempty(elements)
        %effets(1,end) = facteurs(i);
        %fprintf('mean:%f\n',facteurs(i));
    elseif size(elements,2)==1
        effets(elements,elements) = effets(elements,elements)+ facteurs(i)^2/tot;
    elseif size(elements,2)==2
        effets(max(elements),min(elements)) = effets(max(elements),min(elements))+ facteurs(i)^2/tot;
    elseif size(elements,2)==3
        effets(1,end) = effets(1,end)+ facteurs(i)^2/tot;
    end
    
end

% Création du meta modèle
B=char(sum(famille.*facteurs));
B=regexprep(B, 'X\(([0-9])\)', 'X(:,$1)');
B=regexprep(B, '([*/^])', '.$1');
model=inline(B, 'X');

%Effets
exact
effets

%Err
R2= mean( (test_output-model(test_input)).^2 )










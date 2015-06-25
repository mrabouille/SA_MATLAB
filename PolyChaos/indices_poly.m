function [mo,id]=indices_poly(deg_max,nb_poly,poly_cur, debutLigne, id, mo)
%INDICES_POLY génère une indexation lexicographique par degré total croissant
% 
%   mo = indices_poly(3,4)
%       deg_max    : degré maximal du polynome final
%       nb_poly    : nombre de polynôme à indexer
%       ...opt     : suite des paramètres, ne sert que pour la recurence
% 
%       mo         : matrice des ordres
%       id         : nombre de polynômes produit
% 
% 

if nargin<=2
    id=0;
    mo=zeros(factorial(deg_max+nb_poly)/(factorial(deg_max)*factorial(nb_poly)),nb_poly, 'uint8');
    
    % Augmentation progressive de l'ordre total des lignes (nombre de degres disponibles)
    for deg=0:deg_max
        [mo,id]=indices_poly(deg,nb_poly,1,[],id,mo);
    end
    
    % +1 à toute les valeurs pour être compatible avec les indices MatLab
    %mo=mo+1;   
    return
end

if poly_cur==nb_poly   
    % dernier polynome de la ligne atteint
    % la ligne est completee par le degre restant defini
    id=id+1;
    mo(id,:)=[debutLigne deg_max];
elseif deg_max==0   
    % degre max atteint 
    % la ligne est completee par des zeros
    id=id+1;
    mo(id,:)=[debutLigne zeros(1,nb_poly-poly_cur+1)];
else
    
    % reduit progressivement l'ordre du polynome courant (poly_cur) et
    % passe au polynome suivant ... qui sera lui meme deduit
    % progressivement avec de nombre de degres restant disponible.
    for cur_deg=deg_max:-1:0
        [mo,id]=indices_poly(deg_max-cur_deg,nb_poly,poly_cur+1, [debutLigne cur_deg], id, mo);
    end
end
end

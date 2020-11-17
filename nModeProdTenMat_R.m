function P = nModeProdTenMat_R(modes, T, A)
%MATxTEN multiplica un tensor por una matriz (proyecta un tensor sobre una 
%base). T es el tensor, A la matriz y dim la dimension en la que se opera.
%El tama�o de T ser� N1xN2x...Nn. El tama�o de A ser� de MxNi donde i es
%el modo o dimensi�n en la cual hay que multiplicar; n�mero argumento de   
%entrada mode. El resultado de la proyecci�n P ser� un tensor del mismo 
%tama�o que T en todas sus diemnsiones salvo en la i-�sima, donde tendr�
%tama�o M en vez de Ni
% 1.) https://doi.org/10.1145/1186785.1186794
% 2.) https://math.stackexchange.com/questions/2005612/tensor-mode-product
% 3.) ISBN 978-3-642-28027-6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Devuleve el producto P de T por A para esa dimensi�n
% modes: dimensiones sobre las que aplicar el producto 1,[1 2 3 ...], [1 5]
% T un tensor
% A Matriz si s�lo se multiplica una dimensi�n o cell de matrices, {A B C}
% Ej: nModeProdTenMat_R(3, T,  A) � nModeProdTenMat_R([1 3], T, {U1 [] U2})
if ~iscell(A)
    aux = cell([ 1 ndims(T)]); 
    for iii = modes
        aux{iii} = A;
    end
    A = aux;
end

ndim = ndims(T);
dimorder = 1:1:ndim;

P = T;

for iii = modes
    %% El tama�o original
    sz = size(P);
    %% y el tama�o del producto
    szP        = sz;
    szP(iii) = size(A{iii},1);
    %% Almacenar las permutaciones de las dimeniones
    shiftdiml = circshift(dimorder,1-iii);
    shiftdiml(2:end)= circshift(shiftdiml(2:end),1);
    for j = ndim:-1:1
        shiftdimr(j) =  find(shiftdiml == j);
    end
    %% Rotar dimensiones y matrizar
    P = permute(P,shiftdiml);
    P = reshape(P, [ sz(iii) numel(P)/sz(iii) ]);
    %% Multiplicar
    P = A{iii} * P;
    %% Volver a ensamblar en forma de tensor
    P = reshape(P, szP(shiftdiml) );
    P = permute(P, shiftdimr );
end
end
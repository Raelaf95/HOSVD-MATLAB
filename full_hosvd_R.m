function T = full_hosvd_R(S, U, DIMS)
%full_hosvd_R reconstruye el tensor original a partir de su descomposición
% HOSVD T = SUM(  U^(i) * S ).
% 1.)   https://doi.org/10.1016/j.expthermflusci.2016.06.017
% 2.)   DOI: 10.1137/S0895479896305696
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% S almacena el nucleo de la descomposición.
% U es un cell de matrices. Cada Matriz almacena los modos 1-D del tensor
% S y U tal como se han generado en la función hosvd_R.m
% Con el tercer argumento, DIMS, indicas en qué dimensiones multiplicar el
% tensor S por la matriz. Ejemplo DIMS = [1 2 3 4... mdims(S)] es lo mismo 
% que full_hosvd_R(S, U). DIMS = [ 3 4 ] sólo multiplica en la 3a y 4a dim.
% T = full_hosvd_R(S, U, DIMS)
T = S;
ndim = length(U);

if nargin < 3
    DIMS = 1:1:ndim;
end

% sz = zeros([1 ndim]); szTrc = sz;
for i = ndim:-1:1
    %% En caso de descomposiciones truncadas, el tamaño original
    sz(i) = size(U{i},1);
    %% y el tamaño truncado size(Struncado) = szTrc
    szTrc(i) = size(U{i},2);
end

dimorder = 1:1:ndim;

siz = szTrc;
for i = DIMS %1:1:ndim
    %% Controlar como se expande el tensor
    siz(i) = sz(i);
    %% Almacenar las permutaciones de las dimeniones
    shiftdiml = circshift(dimorder,1-i);    
    shiftdiml(2:end)= circshift(shiftdiml(2:end),1);
    for j = ndim:-1:1
        shiftdimr(j) =  find(shiftdiml == j);
    end
    %% Rotar dimensiones y matrizar
    T = permute(T,shiftdiml);
    T = reshape(T, [ szTrc(i) numel(T)/szTrc(i) ]);
    %% Multiplicar
    T = U{i} * T;
    %% Volver a ensamblar en forma de tensor
    T = reshape(T, siz(shiftdiml) );
    T = permute(T, shiftdimr );
end

end
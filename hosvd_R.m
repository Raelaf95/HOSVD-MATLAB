function [S, U, sv] = hosvd_R(T)
%hosvd_R generaliza la descomposición SVD de una matriz para un tensor de orden n>=2.
% Fuentes:
% 1.)   https://doi.org/10.1016/j.expthermflusci.2016.06.017
% 2.)   DOI: 10.1137/S0895479896305696
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% VERSION 1.0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% S almacena el núcleo de la descomposición.
% U es un cell de matrices. Cada Matriz almacena los modos 1-D del tensor
% sv son los valores singulares de las matrices B

ndim = ndims(T);
sz   = size(T);

S  = T;
U  = cell(1,ndim);
sv = cell(1,ndim);

dimorder = 1:1:ndim;

for i = 1:1:ndim % Para la dimensión i    
    shiftdiml = circshift(dimorder,1-i);    
    shiftdiml(2:end)= circshift(shiftdiml(2:end),1);
    for j = ndim:-1:1
        shiftdimr(j) =  find(shiftdiml == j);
    end
    %% Matricization para construir la matriz B^(i) de [ref. 1]
    B = permute(T,shiftdiml);
    B = reshape( B, [ size(B,1) numel(B)/size(B,1) ] );
    B = B -  .5 * ( min(B,[],2) + max(B,[],2) );
    B2 = B*B';
   
    %% SVD De la matriz B
    [V, D]  = eig(B2);
    [D,idx] = sort(diag(D),'descend');
    U{i}    = V(:, idx);    
    sv{i}   = sqrt(D);
% %     [U{i}, D, ~] = SVD_R(B); 
%     [U{i}, D, ~] = svd(B,'econ');
%     sv{i}   = diag(D);
    %% Para calcular S hay que matrizar el tensor original T
    S = permute(S,shiftdiml);
    S = reshape( S, [ size(S,1) numel(S)/size(S,1) ] );
    %% Multiplicar por la conjugada transpuesta de U (despejar U*S = T)
    S = U{i}' * S;
    %% Y volver a ensamblar en forma de tensor
    S = reshape(S, sz(shiftdiml) );
    S = permute(S, shiftdimr );
end

end
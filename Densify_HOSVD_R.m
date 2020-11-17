function T =  Densify_HOSVD_R(T,NewSize,method,VEC)
% VEC trunca. Densify_HOSVD_R(T,NewSize,method) no trunca
%Desify_HOSVD_R Interpola un tensor devolviendo otro de un tamaño distinto 
%N1xN2...xNn
% method = 'linear' 'spline' 'makima' ó Densify_HOSVD_R(T,NewSize,VEC)
if nargin < 3
    method = 'makima';
end

SzOriginal = size(T);
    
if nargin <= 3 || any(VEC == -1) % Si no hay VEC, no trunco    
    [S, U, ~] = hosvd_R(T);
    for i = ndims(T):-1:1
        
        if NewSize(i)~= SzOriginal(i)
            %% 1D
            x  = linspace( 1, size(U{i},1), size(U{i},1) );
            xq = linspace( 1, size(U{i},1), NewSize(i)   );
            for j=1:1:size(U{i},2)
                U_int{i}(:,j) = interp1(x,U{i}(:,j),xq,method);
            end
            %% 2D
%             [x, y]  = meshgrid(linspace( 1, size(U{i},2), size(U{i},2) ));
%             [xq, yq] =  meshgrid(linspace( 1, size(U{i},1),size(U{i},2) ),linspace( 1, size(U{i},2), NewSize(i) ) );   
%             U_int{i} = interp2(x,y,U{i},xq,yq,method); 
        else
                U_int{i} = U{i};
        end
    end
    
    T = full_hosvd_R(S, U_int);
    
elseif nargin == 4 % Si hay VEC trunco. Si VEC = [] elegiré a dedo    
    [ ~, S, U, ~ ] = trunc_hosvd_R(VEC,T);
    
    for i = ndims(T):-1:1
        
        if NewSize(i)~= SzOriginal(i)
            x  = linspace( 1, size(U{i},1), size(U{i},1) );
            xq = linspace( 1, size(U{i},1), NewSize(i)   );
            
            for j=1:1:size(U{i},2)
                U_int{i}(:,j) = interp1(x,U{i}(:,j),xq,method);
            end
        else 
                U_int{i} = U{i};
        end
    end
    
    T = full_hosvd_R(S, U_int);
    
else 
    warning('Estás utilizando mal Desify_HOSVD_R pollito')
end

end


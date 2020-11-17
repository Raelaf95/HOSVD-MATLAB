function P = Tensor3DTimes3MatricesFORS(T,U1,U2,U3)
%Tensor3DTimes3MatricesFORS es un ejemplo de cómo se puede particularizar
%el algoritmo de la función nModeProdTenMat_R para multiplicar un tensor de
%orden tres, T, por tres matrices U1, U2, U3 mediante bulcles. Si:
% P1 = Tensor3DTimes3MatricesFORS(T,U1,U2,U3), y 
% P2 =  nModeProdTenMat_R([1 2 3],T, { U1 U2 U3 }), entonces,
% P1 == P2 -> TRUE.

for j1 = size(U1,1):-1:1
    for j2 = size(U2,1):-1:1
        for j3 = size(U3,1):-1:1
            P(j1,j2,j3)=0 ;
            for i = size(U1,2):-1:1
                for j = size(U2,2):-1:1
                    for k = size(U3,2):-1:1
                        AUX = T(i,j,k) * U1(j1,i) * U2(j2,j) * U3(j3,k) ;                       
                        P(j1,j2,j3)=P(j1,j2,j3) + AUX ;
                    end
                end
            end
        end
    end
end

end


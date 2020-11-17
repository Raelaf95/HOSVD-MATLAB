function [ T_t, S_t, U_t, SV_t, VEC  ] = trunc_hosvd_R(VEC,A,U,sv)
% trunc_hosvd_R(VEC,A,U,sv) trunca un tensor.
% trunc_hosvd_R(VEC,A) trunca el tensor A en los modos del vector VEC. 
% trunc_hosvd_R(VEC,S,U,sv) trunca el nucleo S y los modos unidimensionales
% U en los modos del vector VEC. Devuleve el tensor denso A = full(S,U)
% Si VEC = [], eliges a manualemente los valores.

T_t=[]; S_t=[]; U_t=[]; SV_t=[];

if VEC == -1
    T_t = A;
    if nargout > 1
        [S_t, U_t, SV_t] = hosvd_R(A);
    end
    return
end

if nargin == 2 % Si doy trunc_hosvd_R(VEC,T). VEC puede ser []
    [S, U, sv] = hosvd_R(A);
    
    if isempty(VEC)
        
        figure()
%         setPlot
        set(gcf,'units','centimeters','Position',[6 6 16 16/.5/(1+sqrt(5))]); %  16/(.5*(1+sqrt(5)))]);
        set(gca,'Units','normalized','Position',[.0925 .11 .89 .89],'Box','on'); %0.1300    0.1100    0.3347    0.8150
        ylabel('SV','Rotation',0);%,'units','normalized','Position',[-.06 .5 .00])
        xlabel('Ordinal');%,'units','normalized','Position',[.5 -.07 .00])
        hold on
        grid on
        for i=ndims(A):-1:1
            N=length(sv{i}); sv{i} = real(sv{i});
            h(i)=scatter(linspace(1,N,N),sv{i},24,RBW(i,ndims(A)),'filled');
            set(gca, 'YScale', 'log')
            LOL{i}=['Dimensión ' num2str(i)];% '\hspace{5mm}'];
        end
        legend(h,LOL)%,'Interpreter','Latex'),'Fontsize',14)
        dcm_obj = datacursormode(gcf);
        set(dcm_obj,'DisplayStyle','datatip','SnapToDataVertex','on',      ...
            'Enable','on');
        for i=ndims(A):-1:1
            pause()
            c_info= getCursorInfo(dcm_obj);
            VEC(ndims(A)-i+1)=c_info.Position(1);
        end
        datacursormode off;
        close(gcf)
        disp(['***** VEC == ' num2str(VEC)] )
    end
    
    for i = ndims(S):-1:1
        Idx{i} = 1:1:VEC(i);
    end    
    S_t = S(Idx{:});
    
    for i = ndims(A):-1:1
        U_t{i}  = U{i}(:,1:VEC(i));
        SV_t{i} = sv{i}(1:VEC(i));
    end
    
elseif nargin == 3 || 4 % Si doy trunc_hosvd_R(VEC,S,U,sv). VEC puede ser []
    
    if isempty(VEC)
        
        figure()
        hold on
        grid on
        for i=ndims(A):-1:1
            N=length(sv{i}); sv{i} = real(sv{i});
            h(i)=scatter(linspace(1,N,N),sv{i},24,RBW(i,ndims(A)),'filled');
            set(gca, 'YScale', 'log')
            LOL{i}=['Dimensión número ' num2str(i)];
        end
        legend(h,LOL)
        dcm_obj = datacursormode(gcf);
        set(dcm_obj,'DisplayStyle','datatip','SnapToDataVertex','on',      ...
            'Enable','on');
        for i=ndims(A):-1:1
            pause()
            c_info= getCursorInfo(dcm_obj);
            VEC(ndims(A)-i+1)=c_info.Position(1);
        end
        datacursormode off;
        close(gcf)
        disp(['***** VEC == ' num2str(VEC)] )
    end
    
    for i = ndims(A):-1:1
        Idx{i} = 1:1:VEC(i);
    end     
    S_t = A(Idx{:});
    
    for i = ndims(A):-1:1
        U_t{i}  = U{i}(:,1:VEC(i));
        if nargin == 4 
            SV_t{i} = sv{i}(1:VEC(i));
        end
    end
end

if nargout == 4 || nargout == 1
    T_t = full_hosvd_R(S_t, U_t);
end

    function vec = RBW(iii,N)        
        RGB_vecR = [1 , 1 , 1 ,0.5, 0 , 0 , 0 , 0 , 0 ,0.5, 1 ];
        RGB_vecG = [0 ,0.5, 1 , 1 , 1 , 1 , 1 ,0.5, 0 , 0 , 0 ];
        RGB_vecB = [0 , 0 , 0 , 0 , 0 ,0.5, 1 , 1 , 1 , 1 , 1 ];
        RGB_vecx = linspace(-1e-8,1+1e-8,11);
        
        vec = [interp1(RGB_vecx,RGB_vecR,(iii-1)/(N-1)),...
            interp1(RGB_vecx,RGB_vecG,(iii-1)/(N-1)),...
            interp1(RGB_vecx,RGB_vecB,(iii-1)/(N-1))];        
    end

end
%% Sensitivity analysis
% Closed multiregional economy
% No residential damage
% Hypothetical case
% 30%-24 COVID control over all regions
% Large flood in C
% 50% export restriction on MANR-C and MANR-B (retaliation)
% without production specialization
% T = 82
% on inventory restoration rate: 0.2 0.4 0.6 0.8 1
% close all
clear
clc
%% Main loop
%%% scenario sets
% Flood scale
% Intensity of flood on pop, agr, mang, manr, con and com
Intense_Fld = [0.6 0.6 0.3 0.3 0.45 0.45]; % percent of damage
% Export restriction
imc = 0.5;
RestrMatrix = [3 2;3 3]; % Restricted sectors
% Inventory restoration rate
tauss = [0.2 0.4 0.6 0.8 1];

%%% load data
R = 4;  % region
S = 5;  % sector
F = 2;  % number of productive factors
T = 82;

MRIOT = xlsread('MRIOT_test.xlsx','Sheet1','C3:AA25'); % weekly flow; 10000 RMB, 2015
CurrencyEX = xlsread('MRIOT_test.xlsx','Sheet2','BH13'); % 1RMB = 0.1413USD
MRIOT = MRIOT.*CurrencyEX./100000; % billion USD,2015
MRIOT(isnan(MRIOT)) = 0;

CapitalMatrix = xlsread('Capital_Matrix.xlsx','Sheet1','C2:F21'); % intermediate input needed to produce one unit of capital in each region
CapitalMatrix = CapitalMatrix(:,sort(repmat(1:R,1,S)));

VAL = zeros(size(tauss,2),S*R);
VAL_T = zeros(size(tauss,2),T);
VA_AT = zeros(size(tauss,2),T);
VA_BT = zeros(size(tauss,2),T);
VA_CT = zeros(size(tauss,2),T);
VA_DT = zeros(size(tauss,2),T);

for ss = 1:size(tauss,2)
    
    ss
    
    zz = MRIOT(1:R*S, 1:R*S);
    ff = MRIOT(1:R*S, R*S+1:end-1);
    va = MRIOT(R*S+1:R*S+F,1:R*S); % capital; labor
    xx = MRIOT(end,1:R*S);
    ex = zeros(S*R,1);
    for i = 1:R
        ex(1+(i-1)*S:i*S)=sum(zz(1+(i-1)*S:i*S,:),2)+sum(ff(1+(i-1)*S:i*S,:),2)...
            -sum(zz(1+(i-1)*S:i*S,1+(i-1)*S:i*S),2)-ff(1+(i-1)*S:i*S,i);
    end
        
    va_region = zeros(1,R);
    for i = 1:R
        va_region(i) = sum(sum(va(:,1+S*(i-1):S*i)));
    end
    va_sum = sum(va_region);
    
    capital_0 = va(1,:)*52*4;
    
    I_sum = repmat(eye(S),1,R);
    E_C_Capital = I_sum * CapitalMatrix;
    CapitalMatrix_dis = CapitalMatrix ./ repmat(E_C_Capital,R,1);
    CapitalMatrix_dis(isnan(CapitalMatrix_dis))=0;
    
    si = RestrMatrix(:,1); % Export restrictions on specific goods
    sr = RestrMatrix(:,2);
    
    %%% preprocess data
    zz_c = I_sum * zz;  % sum of intermediate use by sector
    z_dis = zz ./ repmat(zz_c,R,1);  % distribution of intermediate use by sector among regions
    z_dis(isnan(z_dis)) = 0;
    
    ff_c = I_sum * ff;  % sum of final use by sector
    f_dis = ff ./ repmat(ff_c,R,1);
    f_dis(isnan(f_dis)) = 0;
    
    E_CZ = zz_c ./ xx;   % intermediate use by one unit of output
    E_CZ(isnan(E_CZ)) = 0;
    E_VA = va ./ xx;     % value-added use by one unit of output
    E_VA(isnan(E_VA)) = 0;
    
    mat_trivial = E_CZ<0.001;
    
    %%%% Overproduction module
    OverProd = ones(F,R*S);
    OverProdSign = zeros(F,R*S);
    OverProdStep = repmat(0.25/52,F,R*S);
    OverProdUpBd = repmat(1.25,F,R*S);
    
    %%%% stock & order
    StockP = 4;  % weeks of production sustained by stock
    StockMat = StockP*zz_c;
    StockObj = StockMat+zz_c;
    OrderMat = [zz,ff];
    
    %%% Initialization
    Rcod = [1,2,4];
    Rfld = [];
    Rcfd = 3;
    Roth = [];
    
    Sagr = 1;  % 0.5(c);1(f);
    Sman = [2,3,4];  % 1; construction has the same labour-impact coefficient with manufacture
    Scom = 5;  % 0.1;
    
    Stran = 5;  % Ids of transport sectors = 5
    p_rescue = 0.05; % portion of rescue products in manufacture general
    p_outdoor = 0.2; % portion of outdoor services in commerce
    p_health = 0.1;
    
    Duration_Cod = 24;  % weeks
    Intense_Cod = 0.3; % Intensity of transport control
    Intense_Cod_ex = 0.3; % Intensity of export control to other regions
    Intense_Cod_F = 0.3; % percentage of food and outdoor consumption declines
    
    Intense_Cod_Labor = zeros(1,S*R); % Proportion of Covid-affected labor
    Intense_Cod_Tran = zeros(1,R);
    Intense_Cod_Tran_ex = zeros(1,R);
    Intense_Cod_F_R = zeros(1,R);
    for i = [Rcod,Rcfd]
        Intense_Cod_Labor(1+S*(i-1):S*i) = Intense_Cod;
        Intense_Cod_Tran(i) = Intense_Cod;
        Intense_Cod_Tran_ex(i) = Intense_Cod_ex;
        Intense_Cod_F_R(i) = Intense_Cod_F;
    end
    
    Duration_Fld = 2;  % weeks at maximum level
    Intense_Fld_F = Intense_Fld(1); % percentage of outdoor consumption declines
    
    Intense_Fld_Labor = zeros(1,S*R); % Proportion of flood-affected labor
    Intense_Fld_Capital = zeros(1,S*R);
    Intense_Fld_Tran = zeros(1,R);
    Intense_Fld_F_R = zeros(1,R);
    for i = [Rfld,Rcfd]
        Intense_Fld_Labor(1+S*(i-1):S*i) = Intense_Fld(1);
        Intense_Fld_Capital(1+S*(i-1):S*i) = Intense_Fld(2:end);
        Intense_Fld_Tran(i) = Intense_Fld(1+Stran);
        Intense_Fld_F_R(i) = Intense_Fld_F;
    end
    
    %%%% Constraints
    Cod_Labor = zeros(T, S*R);
    for i = 1:S*R
        if Intense_Cod_Labor(i) ~= 0
            Cod_Labor(:,i) = -[zeros(4,1);...
                repmat(Intense_Cod_Labor(i),Duration_Cod,1);(Intense_Cod_Labor(i)...
                :(-Intense_Cod_Labor(i)/4):0)';zeros(T-Duration_Cod-9,1)];
        end
    end
    
    Cod_Tran = zeros(T, R);
    Cod_Tran_ex = zeros(T, R);
    Cod_F = zeros(T, R);
    for i = [Rcod,Rcfd]
        if Intense_Cod_Tran(i) ~= 0
            Cod_Tran(:,i) = -[zeros(4,1);...
                repmat(Intense_Cod_Tran(i),Duration_Cod,1);(Intense_Cod_Tran(i)...
                :(-Intense_Cod_Tran(i)/4):0)';zeros(T-Duration_Cod-9,1)];
        end
        if Intense_Cod_Tran_ex(i) ~= 0
            Cod_Tran_ex(:,i) = -[zeros(4,1);...
                repmat(Intense_Cod_Tran_ex(i),Duration_Cod,1);(Intense_Cod_Tran_ex(i)...
                :(-Intense_Cod_Tran_ex(i)/4):0)';zeros(T-Duration_Cod-9,1)];
        end
        if Intense_Cod_F_R(i) ~= 0
            Cod_F(:,i) = -[zeros(4,1);...
                repmat(Intense_Cod_F_R(i),Duration_Cod,1);(Intense_Cod_F_R(i)...
                :(-Intense_Cod_F_R(i)/4):0)';zeros(T-Duration_Cod-9,1)];
        end
    end
    
    Fld_Capital = -[zeros(4,S*R);Intense_Fld_Capital;zeros(T-4-1,S*R)];
    
    Fld_Labor = zeros(T,S*R);
    for i = 1:S*R
        if Intense_Fld_Labor(i) ~=0
            Fld_Labor(:,i) = -[zeros(4,1);repmat(Intense_Fld_Labor(i),Duration_Fld,1);...
                (Intense_Fld_Labor(i):(-Intense_Fld_Labor(i)/2):0)';zeros(T-Duration_Fld-7,1)];
        end
    end
    
    Fld_Tran = zeros(T,R);
    Fld_F = zeros(T,R);
    for i = [Rfld,Rcfd]
        if Intense_Fld_Tran(i) ~= 0
            Fld_Tran(:,i) = -[zeros(4,1);repmat(Intense_Fld_Tran(i),Duration_Fld,1);...
                (Intense_Fld_Tran(i):(-Intense_Fld_Tran(i)/2):0)';zeros(T-Duration_Fld-7,1)];
        end
        if Intense_Fld_F_R(i) ~= 0
            Fld_F(:,i) = -[zeros(4,1);repmat(Intense_Fld_F_R(i),Duration_Fld,1);...
                (Intense_Fld_F_R(i):(-Intense_Fld_F_R(i)/2):0)';zeros(T-Duration_Fld-7,1)];
        end
    end
    
    imc_t = zeros(T,2);
    if imc~=0
        imc_t(:,1) = -[zeros(4,1);repmat(imc,Duration_Cod,1);...
            (imc:(-imc/4):0)';zeros(T-Duration_Cod-9,1)];
        imc_t(:,2) = -[zeros(6,1);repmat(imc,Duration_Cod-2,1);...
            (imc:(-imc/4):0)';zeros(T-Duration_Cod-9,1)];
    end
    
    Labor_Cons = zeros(T,S*R);
    for i = 1:S
        if ismember(i,Sagr)
            Labor_Cons(:,i:S:end) = ones(T,R) + min(0.5*Cod_Labor(:,i:S:end), Fld_Labor(:,i:S:end));
        elseif ismember(i,Sman)
            Labor_Cons(:,i:S:end) = ones(T,R) + min(Cod_Labor(:,i:S:end), Fld_Labor(:,i:S:end));
        else
            Labor_Cons(:,i:S:end) = ones(T,R) + 0.1*min(Cod_Labor(:,i:S:end), Fld_Labor(:,i:S:end));
        end
    end
    
    Tran_Cons = zeros(S*R,R,T);
    for i = 1:S
        for j = 1:R
            for k = 1:R
                if ismember(i,Sagr)
                    if j == k
                        Tran_Cons(i+(j-1)*S,k,:) = ones(1,1,T)+ reshape(min([0.5*Cod_Tran(:,j),Fld_Tran(:,j)],[],2),1,1,T);
                    else
                        Tran_Cons(i+(j-1)*S,k,:) = ones(1,1,T)+ reshape(min([Fld_Tran(:,j),0.5*Cod_Tran_ex(:,j)],[],2),1,1,T);
                    end
                elseif ismember(i,Sman)
                    if j == k
                        Tran_Cons(i+(j-1)*S,k,:) = ones(1,1,T)+ reshape(min([Cod_Tran(:,j),Fld_Tran(:,j)],[],2),1,1,T);
                    else
                        Tran_Cons(i+(j-1)*S,k,:) = ones(1,1,T)+ reshape(min([Fld_Tran(:,j),Cod_Tran_ex(:,j)],[],2),1,1,T);
                    end
                else
                    if j == k
                        Tran_Cons(i+(j-1)*S,k,:) = ones(1,1,T)+ reshape(0.1*min([Cod_Tran(:,j),Fld_Tran(:,j)],[],2),1,1,T);
                    else
                        Tran_Cons(i+(j-1)*S,k,:) = ones(1,1,T)+ reshape(0.1*min([Fld_Tran(:,j),Cod_Tran_ex(:,j)],[],2),1,1,T);
                    end
                end
            end
        end
    end
    
    FF_T = repmat(ff,1,1,T);
    for j = Rcod
        for k = 1:R
            if k ~= j
                FF_T(5+(j-1)*S,k,:) = FF_T(5+(j-1)*S,k,:).*reshape(ones(T,1)+p_outdoor*Cod_F(:,j),1,1,T);
                FF_T(5+(k-1)*S,j,:) = FF_T(5+(k-1)*S,j,:).*reshape(ones(T,1)+p_outdoor*Cod_F(:,j)-p_health*0.1*Cod_F(:,j),1,1,T);
            else
                FF_T(5+(j-1)*S,k,:) = FF_T(5+(j-1)*S,k,:).*reshape(ones(T,1)+p_outdoor*Cod_F(:,j)-p_health*0.2*Cod_F(:,j),1,1,T);
            end
        end
    end
    for j = Rcfd
        for k = 1:R
            if k ~= j
                FF_T(5+(j-1)*S,k,:) = FF_T(5+(j-1)*S,k,:).*reshape(ones(T,1)+p_outdoor*min(Cod_F(:,j),Fld_F(:,j)),1,1,T);
                FF_T(5+(k-1)*S,j,:) = FF_T(5+(k-1)*S,j,:).*reshape(ones(T,1)+p_outdoor*min(Cod_F(:,j),Fld_F(:,j))-p_health*0.1*Cod_F(:,j),1,1,T);
                FF_T(2+(k-1)*S,j,:) = FF_T(2+(k-1)*S,j,:).*reshape(ones(T,1)-p_rescue*0.05*Fld_F(:,j),1,1,T);
            else
                FF_T(5+(j-1)*S,k,:) = FF_T(5+(j-1)*S,k,:).*reshape(ones(T,1)+p_outdoor*min(Cod_F(:,j),Fld_F(:,j))-p_health*0.2*Cod_F(:,j),1,1,T);
                FF_T(2+(k-1)*S,j,:) = FF_T(2+(k-1)*S,j,:).*reshape(ones(T,1)-p_rescue*0.1*Fld_F(:,j),1,1,T);
            end
        end
    end
    for j = Rfld
        for k = 1:R
            FF_T(5+(j-1)*S,k,:) = FF_T(5+(j-1)*S,k,:).*reshape(ones(T,1)+p_outdoor*Fld_F(:,j),1,1,T);
            if k ~= j
                FF_T(5+(k-1)*S,j,:) = FF_T(5+(k-1)*S,j,:).*reshape(ones(T,1)+p_outdoor*Fld_F(:,j),1,1,T);
                FF_T(2+(k-1)*S,j,:) = FF_T(2+(k-1)*S,j,:).*reshape(ones(T,1)-p_rescue*0.05*Fld_F(:,j),1,1,T);
            else
                FF_T(2+(k-1)*S,j,:) = FF_T(2+(k-1)*S,j,:).*reshape(ones(T,1)-p_rescue*0.1*Fld_F(:,j),1,1,T);
            end
        end
    end
    
    Im_cons = ones(S*R,R,T);
    for k = 1:T
        for i = 1:size(RestrMatrix,1)
            if ismember(sr(i),[Rfld,Rcfd])
                Im_cons(si(i)+S*(sr(i)-1),:,k) = ones(1,R)+repmat(imc_t(k,1),1,R);
                Im_cons(si(i)+S*(sr(i)-1),sr(i),k) = 1;
            else
                Im_cons(si(i)+S*(sr(i)-1),:,k) = ones(1,R)+repmat(imc_t(k,2),1,R);
                Im_cons(si(i)+S*(sr(i)-1),sr(i),k) = 1;
            end
        end
    end
    
    %%% Modelling
    x_t = zeros(T,S*R);
    x_t_max = zeros(T,S*R);
    va_sum_t = zeros(T,S*R);
    f_t = zeros(S*R,R,T);
    
    x_t(1,:) = xx;
    x_t_max(1,:) = xx;
    va_sum_t(1,:) = sum(va);
    f_t(:,:,1) = ff;
    
    Capital_Cons = ones(T,S*R);
    
    RA = zeros(T,S*R);
    K = zeros(T,S*R);
    K(1,:) = capital_0;
    RD = zeros(S*R);
    
    for t = 2:T
        t
        
        K(t,:) = K(t-1,:)+Fld_Capital(t,:).*K(t-1,:)+RA(t-1,:);
        Capital_Cons(t,:) = min(K(t,:)./capital_0,1);
        
        VA_Cons = [Capital_Cons(t,:);Labor_Cons(t,:)];
        VA_t = VA_Cons.*va.*OverProd;
        
        x_max = Production_max(StockMat,E_CZ,VA_t,E_VA,xx,S,mat_trivial);
        x = Production(StockMat,E_CZ,VA_t,E_VA,OrderMat,RD,xx,S,mat_trivial);
        
        OrderInter = sum(OrderMat(:,1:S*R),2);
        F_DisMat = [6*OrderMat(:,S*R+1:end),RD];
        x_dis_inter = zeros(S*R);
        x_remain = zeros(S*R,1);
        for i = 1:S*R
            if OrderInter(i) <= x(i)
                x_dis_inter(i,:) = OrderMat(i,1:S*R);
                x_remain(i) = x(i)-OrderInter(i);
            else
                x_dis_inter(i,:) = OrderMat(i,1:S*R)./OrderInter(i).*x(i);
                x_remain(i) = 0;
            end
        end
        x_dis_final = F_DisMat(:,1:R)./sum(F_DisMat,2).*x_remain;
        x_dis_final(isnan(x_dis_final)) = 0;
        x_dis_rec = F_DisMat(:,R+1:end)./sum(F_DisMat,2).*x_remain;
        x_dis_rec(isnan(x_dis_rec)) = 0;
        
        x_dis_inter_ex = x_dis_inter;
        x_dis_rec_ex = x_dis_rec;
        x_dis_final_ex = x_dis_final;
        for i = 1:R
            x_dis_inter_ex(1+(i-1)*S:i*S,1+(i-1)*S:i*S)=0;
            x_dis_rec_ex(1+(i-1)*S:i*S,1+(i-1)*S:i*S)=0;
            x_dis_final_ex(1+(i-1)*S:i*S,i)=0;
        end
        x_dis_inter_loc = x_dis_inter-x_dis_inter_ex;
        x_dis_rec_loc = x_dis_rec-x_dis_rec_ex;
        x_dis_final_loc = x_dis_final-x_dis_final_ex;
        
        Im_cons_min = min(Im_cons(:,:,t),[],2);
        Im_cons_min(Im_cons_min==1) = Inf;
        x_dis_inter_ex_adj = x_dis_inter_ex;
        x_dis_final_ex_adj = x_dis_final_ex;
        x_dis_rec_ex_adj = x_dis_rec_ex;
        x_dis_inter_loc_adj = x_dis_inter_loc;
        x_dis_final_loc_adj = x_dis_final_loc;
        x_dis_rec_loc_adj = x_dis_rec_loc;
        for i = 1:S*R
            if sum(x_dis_inter_ex(i,:)) >= Im_cons_min(i)*ex(i)
                x_dis_inter_ex_adj(i,:) = x_dis_inter_ex(i,:)./sum(x_dis_inter_ex(i,:)).*Im_cons_min(i).*ex(i);
                x_dis_final_ex_adj(i,:) = zeros(1,R);
                x_dis_rec_ex_adj(i,:) = zeros(1,S*R);
            elseif sum(x_dis_inter_ex(i,:))+sum(x_dis_final_ex(i,:))+sum(x_dis_rec_ex(i,:)) >= Im_cons_min(i)*ex(i)
                x_dis_inter_ex_adj(i,:) = x_dis_inter_ex(i,:);
                x_dis_final_ex_adj(i,:) = x_dis_final_ex(i,:)./(sum(x_dis_final_ex(i,:))+sum(x_dis_rec_ex(i,:))).*(Im_cons_min(i).*ex(i)-sum(x_dis_inter_ex_adj(i,:)));
                x_dis_rec_ex_adj(i,:) = x_dis_rec_ex(i,:)./(sum(x_dis_final_ex(i,:))+sum(x_dis_rec_ex(i,:))).*(Im_cons_min(i).*ex(i)-sum(x_dis_inter_ex_adj(i,:)));
            end
        end
        x_dis_inter_ex_adj(isnan(x_dis_inter_ex_adj)) = 0;
        x_dis_final_ex_adj(isnan(x_dis_final_ex_adj)) = 0;
        x_dis_rec_ex_adj(isnan(x_dis_rec_ex_adj)) = 0;
        
        for i = 1:S*R
            x_loc = x(i)-(sum(x_dis_inter_ex_adj(i,:))+sum(x_dis_final_ex_adj(i,:))+sum(x_dis_rec_ex_adj(i,:)));
            x_dis_inter_loc_adj(i,:) = x_dis_inter_loc(i,:)./(sum(x_dis_inter_loc(i,:))+sum(x_dis_final_loc(i,:))+sum(x_dis_rec_loc(i,:))).*x_loc;
            x_dis_final_loc_adj(i,:) = x_dis_final_loc(i,:)./(sum(x_dis_inter_loc(i,:))+sum(x_dis_final_loc(i,:))+sum(x_dis_rec_loc(i,:))).*x_loc;
            x_dis_rec_loc_adj(i,:) = x_dis_rec_loc(i,:)./(sum(x_dis_inter_loc(i,:))+sum(x_dis_final_loc(i,:))+sum(x_dis_rec_loc(i,:))).*x_loc;
        end
        x_dis_inter_loc_adj(isnan(x_dis_inter_loc_adj)) = 0;
        x_dis_final_loc_adj(isnan(x_dis_final_loc_adj)) = 0;
        x_dis_rec_loc_adj(isnan(x_dis_rec_loc_adj)) = 0;
        
        x_dis_inter = x_dis_inter_ex_adj+x_dis_inter_loc_adj;
        x_dis_final = x_dis_final_ex_adj+x_dis_final_loc_adj;
        x_dis_rec = x_dis_rec_ex_adj+x_dis_rec_loc_adj;
        
        temp = Tran_Cons(:,:,t);
        temp = temp(:,sort(repmat(1:R,1,S)));
        temp_im = Im_cons(:,:,t);
        temp_im = temp_im(:,sort(repmat(1:R,1,S)));
        temp_CapitalMatrix_dis = CapitalMatrix_dis.*temp.*temp_im.*min(VA_Cons,[],1)';
        temp_CapitalMatrix_dis_sum = I_sum*temp_CapitalMatrix_dis;
        CapitalMatrix_dis_u = temp_CapitalMatrix_dis./repmat(temp_CapitalMatrix_dis_sum,R,1);
        CapitalMatrix_dis_u(isnan(CapitalMatrix_dis_u)) = 0;
        CapitalMatrix_u = CapitalMatrix_dis_u .* repmat(E_C_Capital,R,1);
        
        RA(t,:) = sum(x_dis_rec);
        capital_gap = max(capital_0-K(t,:)-RA(t,:),0);
        RD = CapitalMatrix_u * diag(capital_gap).*temp.*temp_im;
        
        StockUse = E_CZ*diag(x);
        StockAdd = I_sum*x_dis_inter;
        StockMat = StockMat-StockUse+StockAdd;
        
        taus = tauss(ss);
        StockGap = taus.*(StockObj-StockMat-zz_c)+zz_c;
        StockGap(StockGap<0)=0;
        
        temp_z_dis = z_dis.*temp.*temp_im.*min(VA_Cons,[],1)';
        temp_z_dis_sum = I_sum*temp_z_dis;
        z_dis_u = temp_z_dis./repmat(temp_z_dis_sum,R,1);
        z_dis_u(isnan(z_dis_u)) = 0;
        
        OrderMat(:,1:S*R) = z_dis_u .* repmat(StockGap,R,1);
        a = OrderMat(:,1:S*R);
        b = zz.*temp.*temp_im;
        a(b<a) = b(b<a);
        OrderMat(:,1:S*R) = a;
        
        temp = Tran_Cons(:,:,t);
        temp_im = Im_cons(:,:,t);
        temp_f_dis = f_dis.*temp.*temp_im.*min(VA_Cons,[],1)';
        temp_f_dis_sum = I_sum*temp_f_dis;
        f_dis_u = temp_f_dis ./ repmat(temp_f_dis_sum,R,1);
        f_dis_u(isnan(f_dis_u)) = 0;
        
        a = I_sum*FF_T(:,:,t);
        OrderMat(:,S*R+1:end) = repmat(a,R,1).*f_dis_u;
        a = OrderMat(:,S*R+1:end);
        b = ff.*temp.*temp_im;
        a(b<a) = b(b<a);
        OrderMat(:,S*R+1:end) = a;
        a = OrderMat(:,S*R+1:end);
        b = FF_T(:,:,t);
        a(b<a) = b(b<a);
        OrderMat(:,S*R+1:end) = a;
        
        OverProdSign = OverProdSignFun(StockMat,E_CZ,VA_t,E_VA,OrderMat,RD,xx,F,S,R);
        OverProd = OverProd + OverProdSign .* OverProdStep;
        OverProd(OverProd<1) = 1;
        OverProd(OverProd>OverProdUpBd) = OverProdUpBd(OverProd>OverProdUpBd);
        
        x_t(t,:) = x;
        x_t_max(t,:) = x_max;
        va_sum_t(t,:) = x - sum(StockUse);
        % va_sum_t(t,:) = sum(E_VA*diag(x));
        f_t(:,:,t) = x_dis_final;
    end
    
    %%% Result analysis
    ValueAddedLoss = sum(va,1)-va_sum_t; % T*SR
    
    VAL(ss,:) = -sum(ValueAddedLoss,1);
    VAL_T(ss,:) = -sum(ValueAddedLoss,2)./va_sum;
    VA_AT(ss,:) = (sum(va_sum_t(:,1:S),2)-va_region(1))./va_region(1); % Region A
    VA_BT(ss,:) = (sum(va_sum_t(:,S+1:2*S),2)-va_region(2))./va_region(2); % Region B
    VA_CT(ss,:) = (sum(va_sum_t(:,2*S+1:3*S),2)-va_region(3))./va_region(3);
    VA_DT(ss,:) = (sum(va_sum_t(:,3*S+1:end),2)-va_region(4))./va_region(4); % Region D
end
%% Result analysis
VAL_s = -sum(VAL,2); % positive values
% VAL_p = VAL_s./(52.*va_sum);
VAL_A = -sum(VAL(:,1:S),2);   % positive values
VAL_B = -sum(VAL(:,S+1:2*S),2);   % positive values
VAL_C = -sum(VAL(:,2*S+1:3*S),2); % positive values
VAL_D = -sum(VAL(:,3*S+1:end),2); % positive values
%% Figure
dx = 1:T;

linespec = ["-","-","-","-","-","--"];

colorspec = [...
    182, 115, 101;... % brown
    204, 185, 116;... % yellow
    184, 135, 195;... % purple
    100, 181, 205;... % blue
    0, 0, 0;...       % black    
    164, 160, 155]./255; % grey
colorbg = [...
    242, 242, 242;...
    217, 217, 217;...
    191, 191, 191;...
    166, 166, 166;...
    128, 128, 128]./255;
linewidthspec = [1,1,1,1,1,1];

figure(2);
clf;

set(gcf,'Units','centimeter','Position',[1 0.8 15.92 15.92]); % set figure size

f = tiledlayout('flow');
f.TileSpacing = 'compact';
f.Padding = 'compact';

% Region A
nexttile
dy = 100.*[zeros(1,T);VA_AT];
p1 = plot(dx,dy);
for j = 1:size(tauss,2)+1
    if j ==1
        p1(j).Color = colorspec(end,:);
        p1(j).LineStyle = linespec(end);
        p1(j).LineWidth = linewidthspec(end);
    else
        p1(j).Color = colorspec(j-1,:);
        p1(j).LineStyle = linespec(j-1);
        p1(j).LineWidth = linewidthspec(j-1);
    end
end
ylim([-60,10])
xlim([0,T])
set(gca,'YTick',-60:20:0,'XTick',0:20:80,'Layer','top','FontName','Helvetica','Fontsize',8,'LineWidth',0.5,'Color',[colorbg(2,:),0.3]);
box off
ylabel('Change in GDP (%, A)')
xlabel('Weeks')

% Region B
nexttile
dy = 100.*[zeros(1,T);VA_BT];
p1 = plot(dx,dy);
for j = 1:size(tauss,2)+1
    if j ==1
        p1(j).Color = colorspec(end,:);
        p1(j).LineStyle = linespec(end);
        p1(j).LineWidth = linewidthspec(end);
    else
        p1(j).Color = colorspec(j-1,:);
        p1(j).LineStyle = linespec(j-1);
        p1(j).LineWidth = linewidthspec(j-1);
    end
end
ylim([-60,10])
xlim([0,T])
set(gca,'YTick',-60:20:0,'XTick',0:20:80,'Layer','top','FontName','Helvetica','Fontsize',8,'LineWidth',0.5,'Color',[colorbg(2,:),0.3]);
box off
ylabel('Change in GDP (%, B)')
xlabel('Weeks')

% Region C
nexttile
dy = 100.*[zeros(1,T);VA_CT];
p1 = plot(dx,dy);
for j = 1:size(tauss,2)+1
    if j ==1
        p1(j).Color = colorspec(end,:);
        p1(j).LineStyle = linespec(end);
        p1(j).LineWidth = linewidthspec(end);
    else
        p1(j).Color = colorspec(j-1,:);
        p1(j).LineStyle = linespec(j-1);
        p1(j).LineWidth = linewidthspec(j-1);
    end
end
ylim([-60,10])
xlim([0,T])
set(gca,'YTick',-60:20:0,'XTick',0:20:80,'Layer','top','FontName','Helvetica','Fontsize',8,'LineWidth',0.5,'Color',[colorbg(2,:),0.3]);
box off
ylabel('Change in GDP (%, C)')
xlabel('Weeks')

% Region D
nexttile
dy = 100.*[zeros(1,T);VA_DT];
p1 = plot(dx,dy);
for j = 1:size(tauss,2)+1
    if j ==1
        p1(j).Color = colorspec(end,:);
        p1(j).LineStyle = linespec(end);
        p1(j).LineWidth = linewidthspec(end);
    else
        p1(j).Color = colorspec(j-1,:);
        p1(j).LineStyle = linespec(j-1);
        p1(j).LineWidth = linewidthspec(j-1);
    end
end
ylim([-60,10])
xlim([0,T])
set(gca,'YTick',-60:20:0,'XTick',0:20:80,'Layer','top','FontName','Helvetica','Fontsize',8,'LineWidth',0.5,'Color',[colorbg(2,:),0.3]);
box off
ylabel('Change in GDP (%, D)')
xlabel('Weeks')

% World aggregate
nexttile
dy = 100.*[zeros(1,T);VAL_T];
p1 = plot(dx,dy);
for j = 1:size(tauss,2)+1
    if j ==1
        p1(j).Color = colorspec(end,:);
        p1(j).LineStyle = linespec(end);
        p1(j).LineWidth = linewidthspec(end);
    else
        p1(j).Color = colorspec(j-1,:);
        p1(j).LineStyle = linespec(j-1);
        p1(j).LineWidth = linewidthspec(j-1);
    end
end
ylim([-60,10])
xlim([0,T])
set(gca,'YTick',-60:20:0,'XTick',0:20:80,'Layer','top','FontName','Helvetica','Fontsize',8,'LineWidth',0.5,'Color',[colorbg(2,:),0.3]);
box off
ylabel('Change in GDP (%, global)')
xlabel('Weeks')

% legend
lgd = legend('Pre-disaster','InvRate = 0.2','InvRate = 0.4','InvRate = 0.6','InvRate = 0.8','InvRate = 1');
set(lgd,'NumColumns',1,'FontSize',8,'FontName','Helvetica','Box','on','Location','SouthEast');
%% End.


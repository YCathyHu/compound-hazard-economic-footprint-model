%% COVID-flood modelling
% Closed multiregional economy
% No residential damage
% Hypothetical case
% 30%-24 COVID control over all regions
% Multi-scale flooding in region C for 2 weeks
% Specialized production: MANR-C, MANR-B
% Free trade without export restrictions
% T = 82;
% close all
clear
clc
%% Main loop
%%% scenario sets
% Flood scale
Flood_set_int = [...
    0.2 0.2 0.1 0.1 0.15 0.15;...
    0.4 0.4 0.2 0.2 0.3 0.3;...
    0.6 0.6 0.3 0.3 0.45 0.45]; % percent of damage
% specialized production
Scenarios = {[3 2;3 3]};

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

DD = zeros(size(Scenarios,1)*size(Flood_set_int,1),S*R);
FOD = zeros(size(Scenarios,1)*size(Flood_set_int,1),S*R);
OL = zeros(size(Scenarios,1)*size(Flood_set_int,1),S*R);
VAL = zeros(size(Scenarios,1)*size(Flood_set_int,1),S*R);
VA_CT = zeros(size(Scenarios,1)*size(Flood_set_int,1),T);
VA_OT = zeros(size(Scenarios,1)*size(Flood_set_int,1),T);
KK_f = zeros(size(Scenarios,1)*size(Flood_set_int,1),T);
VAL_T = zeros(size(Scenarios,1)*size(Flood_set_int,1),T);
VA_AT = zeros(size(Scenarios,1)*size(Flood_set_int,1),T);
VA_DT = zeros(size(Scenarios,1)*size(Flood_set_int,1),T);

for sc =1:size(Scenarios,1)
    for sf = 1:size(Flood_set_int,1)
        
        [sc sf]
        
        zz = MRIOT(1:R*S, 1:R*S);
        ff = MRIOT(1:R*S, R*S+1:end-1);
        va = MRIOT(R*S+1:R*S+F,1:R*S); % capital; labor
        xx = MRIOT(end,1:R*S);
        
        va_region = zeros(1,R);
        for i = 1:R
            va_region(i) = sum(sum(va(:,1+S*(i-1):S*i)));
        end
        va_sum = sum(va_region);
        
        capital_0 = va(1,:)*52*4;
        
        I_sum = repmat(eye(S),1,R);
        si = Scenarios{sc}(:,1);
        sr = Scenarios{sc}(:,2);
        for i = 1:size(si,1)
            I_sum(si(i),si(i)+(sr(i)-1)*S) = 0;
            I_sum(S+i,si(i)+(sr(i)-1)*S) = 1;  % agr in D is a specialized product
        end
        
        E_C_Capital = I_sum * CapitalMatrix;
        CapitalMatrix_dis = substi(CapitalMatrix,E_C_Capital,I_sum);
        
        %%% preprocess data
        zz_c = I_sum * zz;  % sum of intermediate use by sector
        z_dis = substi(zz,zz_c,I_sum); % distribution of intermediate use by sector among regions
        
        ff_c = I_sum * ff;  % sum of final use by sector
        f_dis = substi(ff,ff_c,I_sum);
        
        E_CZ = zz_c ./ xx;   % intermediate use by one unit of output
        E_CZ(isnan(E_CZ)) = 0;
        E_VA = va ./ xx;     % value-added use by one unit of output
        E_VA(isnan(E_VA)) = 0;
        
        mat_trivial = E_CZ<0.001; %1e-04;
        
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
        Intense_Cod_F = 0.3;
        
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
        Intense_Fld = Flood_set_int(sf,:); % Intensity of flood on pop, agr, mang, manr, con and com
        Intense_Fld_F = Flood_set_int(sf,1); % percentage of outdoor consumption declines
        
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
        
        FirstOrderDamage_Capital = zeros(T,S*R);
        
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
            
            x = Production2(StockMat,E_CZ,VA_t,E_VA,OrderMat,RD,xx,mat_trivial);
            
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
            
            temp = Tran_Cons(:,:,t);
            temp = temp(:,sort(repmat(1:R,1,S)));
            temp_CapitalMatrix_dis = CapitalMatrix_dis.*temp.*min(VA_Cons,[],1)';
            temp_CapitalMatrix_dis_sum = I_sum*temp_CapitalMatrix_dis;
            CapitalMatrix_dis_u = substi(temp_CapitalMatrix_dis,temp_CapitalMatrix_dis_sum,I_sum);
            CapitalMatrix_u = substi_inv(CapitalMatrix_dis_u,E_C_Capital,I_sum);
            
            RA(t,:) = sum(x_dis_rec);
            capital_gap = max(capital_0-K(t,:)-RA(t,:),0);
            RD = CapitalMatrix_u * diag(capital_gap).*temp;
            
            StockUse = E_CZ*diag(x);
            StockAdd = I_sum*x_dis_inter;
            StockMat = StockMat-StockUse+StockAdd;
            StockGap = StockObj - StockMat;
            StockGap(StockGap<0)=0;
            
            temp_z_dis = z_dis.*temp.*min(VA_Cons,[],1)';
            temp_z_dis_sum = I_sum*temp_z_dis;
            z_dis_u = substi(temp_z_dis,temp_z_dis_sum,I_sum);
            
            OrderMat(:,1:S*R) = substi_inv(z_dis_u,StockGap,I_sum);
            a = OrderMat(:,1:S*R);
            b = zz.*temp;
            a(b<a) = b(b<a);
            OrderMat(:,1:S*R) = a;
            
            temp = Tran_Cons(:,:,t);
            temp_f_dis = f_dis.*temp.*min(VA_Cons,[],1)';
            temp_f_dis_sum = I_sum*temp_f_dis;
            f_dis_u = substi(temp_f_dis,temp_f_dis_sum,I_sum);
            
            a = I_sum*FF_T(:,:,t);
            OrderMat(:,S*R+1:end) = substi_inv(f_dis_u,a,I_sum);
            a = OrderMat(:,S*R+1:end);
            b = ff.*temp;
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
            va_sum_t(t,:) = x - sum(StockUse);
            % va_sum_t(t,:) = sum(E_VA*diag(x));
            f_t(:,:,t) = x_dis_final;
            
            FirstOrderDamage_Capital(t,:) = (1-Capital_Cons(t,:)).*sum(va,1);
        end
        
        %%% Result analysis
        DirectDamage = -Fld_Capital.*capital_0;
        FirstOrderDamage_Labor = (1-Labor_Cons).*sum(va,1);
        FirstOrderDamage = max(FirstOrderDamage_Capital,FirstOrderDamage_Labor);
        OutputLoss = xx-x_t;
        ValueAddedLoss = sum(va,1)-va_sum_t;
        
        DD(sf+(sc-1)*size(Flood_set_int,1),:) = -sum(DirectDamage,1);
        FOD(sf+(sc-1)*size(Flood_set_int,1),:) = -sum(FirstOrderDamage,1);
        OL(sf+(sc-1)*size(Flood_set_int,1),:) = -sum(OutputLoss,1);
        VAL(sf+(sc-1)*size(Flood_set_int,1),:) = -sum(ValueAddedLoss,1);
        VAL_T(sf+(sc-1)*size(Flood_set_int,1),:) = -sum(ValueAddedLoss,2)./va_sum;
        VA_CT(sf+(sc-1)*size(Flood_set_int,1),:) = (sum(va_sum_t(:,2*S+1:3*S),2)-va_region(3))./va_region(3);
        %VA_OT(sf+(sc-1)*size(Flood_set_int,1),:) = (sum(va_sum_t(:,[1:2*S,3*S+1:end]),2)-sum(va_region([1,2,4])))./sum(va_region([1,2,4]));
        VA_OT(sf+(sc-1)*size(Flood_set_int,1),:) = (sum(va_sum_t(:,S+1:2*S),2)-va_region(2))./va_region(2); % Region B
        VA_AT(sf+(sc-1)*size(Flood_set_int,1),:) = (sum(va_sum_t(:,1:S),2)-va_region(1))./va_region(1); % Region A
        VA_DT(sf+(sc-1)*size(Flood_set_int,1),:) = (sum(va_sum_t(:,3*S+1:end),2)-va_region(4))./va_region(4); % Region D
        KK_f(sf+(sc-1)*size(Flood_set_int,1),:) = (sum(K(:,2*S+1:3*S),2)-sum(capital_0(2*S+1:3*S)))./sum(capital_0(2*S+1:3*S));
    end
end
%% Result analysis
recovery_time_capital = zeros(size(Scenarios,1)*size(Flood_set_int,1),1);
recovery_time_gdp = zeros(size(Scenarios,1)*size(Flood_set_int,1),1);
for i = 1:size(Scenarios,1)*size(Flood_set_int,1)
    II = find(abs(KK_f(i,:))<1e-03);
    LL = zeros(size(II));
    for j = 1:size(II,2)
        LL(j) = all(abs(KK_f(i,II(j)+1:end))<1e-03);
    end
    if all(LL==0)
        recovery_time_capital(i) = T;
    else
        recovery_time_capital(i) = II(find(LL,1));
    end
end
for i = 1:size(Scenarios,1)*size(Flood_set_int,1)
    II = find(abs(VAL_T(i,:))<1e-03);
    LL = zeros(size(II));
    for j = 1:size(II,2)
        LL(j) = all(abs(VAL_T(i,II(j)+1:end))<1e-03);
    end
    if all(LL==0)
        recovery_time_gdp(i) = T;
    else
        recovery_time_gdp(i) = II(find(LL,1));
    end
end

DD_s = -sum(DD,2);   % positive values
FOD_s = -sum(FOD,2); % positive values
VAL_s = -sum(VAL,2); % positive values

VAL_A = -sum(VAL(:,1:S),2);   % positive values
VAL_B = -sum(VAL(:,S+1:2*S),2);   % positive values
VAL_C = -sum(VAL(:,2*S+1:3*S),2); % positive values
VAL_D = -sum(VAL(:,3*S+1:end),2); % positive values

VAL_C_p = VAL_C./(52.*va_region(3));
VAL_B_p = VAL_B./(52.*va_region(2));
VAL_p = VAL_s./(52.*va_sum);

DD_sorted = resultsort3(DD_s,Scenarios,Flood_set_int);
FOD_sorted = resultsort3(FOD_s,Scenarios,Flood_set_int);
VAL_sorted = resultsort3(VAL_s,Scenarios,Flood_set_int);
rtc_sorted = resultsort3(recovery_time_capital,Scenarios,Flood_set_int);
rtg_sorted = resultsort3(recovery_time_gdp,Scenarios,Flood_set_int);
VAL_A_sorted = resultsort3(VAL_A,Scenarios,Flood_set_int);
VAL_B_sorted = resultsort3(VAL_B,Scenarios,Flood_set_int);
VAL_C_sorted = resultsort3(VAL_C,Scenarios,Flood_set_int);
VAL_D_sorted = resultsort3(VAL_D,Scenarios,Flood_set_int);

TD = DD+VAL;
TD_s = DD_s+VAL_s;
TD_s_p = TD_s./(52.*va_sum);
TD_sorted = DD_sorted+VAL_sorted;
TD_sorted_p = TD_sorted./(52.*va_sum);
%% End.
VA_CT3 = VA_CT;
VA_OT3 = VA_OT;
VA_AT3 = VA_AT;
VA_DT3 = VA_DT;
save FigureData_basefrtrade_sp.mat VA_CT3 VA_OT3 VA_AT3 VA_DT3
function [ A, B, E, G, b, F, c, d, N_conti,N_bi, N_renew, N_para, N_dem ] = model(data_b, data_g, data_l, data_r, data_s, data_gcost,PTDF)

baseMVA = 100;
T = 24;
eff_c  = data_s(:,2);
eff_dc = data_s(:,3);

E_lower = data_s(:,5);
ramp =  data_g(:,11); % 检验  0.5*(P_max-P_min)
N_g = length(data_g(:,1));
N_l = length(data_l(:,1));
N_b = length(data_b(:,1));
N_r = length(data_r(:,1));
N_s = length(data_s(:,1));
%% 计算线性化成本的系数，k和b后面直接用标幺值就可以，得到的是美元。
[ k_gen,b_gen ] = gencost(baseMVA,N_g,data_g,data_gcost);

%% 定义变量 注意顺序！
% 连续变量
P_g = sdpvar(N_g,T,'full'); % 发电功率
Gcost = sdpvar(N_g,T,'full'); % 成本线性化的辅助变量
P_c = sdpvar(N_s,T,'full'); % 储能充电
P_dc = sdpvar(N_s,T,'full'); % 储能放电
E_init = sdpvar(N_s,1,'full'); %储能初始状态
dem = sdpvar(N_b,T,'full'); % 负荷

P_s = sdpvar(N_s,1,'full'); % 参数变量
E_s = sdpvar(N_s,1,'full');% 参数变量
P_r = sdpvar(N_r,T,'full');% 新能源变量
% 布尔量
b_state = binvar(N_g,T,'full');
b_up = binvar(N_g,T,'full');
b_down = binvar(N_g,T,'full');
%% 表示0时段发电机状态，但这个约束不一定要用，可以放开第一时刻的状态及出力。看情况。

b_0 = data_g(:,12);
P_g0 = b_0.*data_g(:,9)*0.5;
%% 最小启停机时间。
minT = data_g(:,13);
%% 负荷曲线。


%% PTDF
C_PTDF = [];
for t = 1:T
    C_PTDF = [C_PTDF; -data_l(:,6) <=( (PTDF(:,data_g(:,1)) * P_g(:,t))...
        - (PTDF * dem(:,t))...
        + (PTDF(:,data_r(:,1)) * P_r(:,t))...
        + (PTDF(:,data_s(:,1)) * (P_dc(:,t)-P_c(:,t)))  )...
        <= data_l(:,6)];
end



%% 发电机
C_gen = [];
%% 启停状态（1时段不要约束了）
for t = 1:1
    C_gen = [C_gen; b_0 - b_state(:,t) + b_up(:,t)>=0];  % 1时段启动
    C_gen = [C_gen; b_state(:,t) - b_0 + b_down(:,t)>=0];  % 1时段停机
end
for t = 2:T
    C_gen = [C_gen; b_state(:,t-1) - b_state(:,t) + b_up(:,t)>=0]; % t时段启动
    C_gen = [C_gen; b_state(:,t) - b_state(:,t-1) + b_down(:,t)>=0]; % t时段停机
end
%% 最小开/停机时间

for i = 1:N_g
    for t = 1:1 %（1时段不要约束了）
        for tau = (t+1):(min(t+minT(i)-1,T)) %这里用t和t+1都行,t实际上是自然满足的
            C_gen = [C_gen; b_state(:,t) - b_0 <= b_state(:,tau)];  % 1时段最小开机
            C_gen = [C_gen; b_0 - b_state(:,t) <= (1-b_state(:,tau))];  % 1时段最小停机
        end
    end
    
    for t = 2:(T-1)  %最后一个时段可以不管以后的状态
        for tau = (t+1):(min(t+minT(i)-1,T)) %这里用t和t+1都行,t实际上是自然满足的
            C_gen = [C_gen; b_state(i,t) - b_state(i,t-1) <= b_state(i,tau)];  % t时段后最小开机
            C_gen = [C_gen; b_state(i,t-1) - b_state(i,t) <= (1-b_state(i,tau))];  % t时段后最小停机
        end
    end
end

%% 爬坡约束
P_g_range = data_g(:,9)-data_g(:,10);
w = 10;
for t = 1:1 %（1时段不要约束了）
    C_gen = [C_gen; P_g(:,t)-P_g0 <= b_0.*P_g_range.*ramp + (1-b_0).*data_g(:,w)];
    C_gen = [C_gen; P_g0-P_g(:,t) <= b_state(:,t).*P_g_range.*ramp + (1-b_state(:,t)).*data_g(:,w)];
end
for t = 2:T
    C_gen = [C_gen; P_g(:,t)-P_g(:,t-1) <= b_state(:,t-1).*P_g_range.*ramp + (1-b_state(:,t-1)).*data_g(:,w)];
    C_gen = [C_gen; P_g(:,t-1)-P_g(:,t) <= b_state(:,t).*P_g_range.*ramp + (1-b_state(:,t)).*data_g(:,w)];
end

%% 功率上下限
for t = 1:T
    C_gen = [C_gen; data_g(:,10).*b_state(:,t) <= P_g(:,t) <= data_g(:,9).*b_state(:,t)];
end






%% 储能
C_sto = [];
%    C_sto = [C_sto; 0 <= P_s <= 0; 0<= E_s <= 0]; %建模时必须去掉这个约束；测试的时候才用。
bigM = 10^5;
for i = 1:N_s
    C_sto = [C_sto; E_s(i)*E_lower(i)/baseMVA<=E_init(i)<=E_s(i)/baseMVA;];
    for t = 1:T
        C_sto = [C_sto; E_s(i)*E_lower(i)/baseMVA<=(E_init(i)+sum(P_c(i,1:t))*eff_c(i)-sum(P_dc(i,1:t))/eff_dc(i))<=E_s(i)/baseMVA];
        C_sto = [C_sto; 0<=P_c(i,t)<=P_s(i)/baseMVA; 0<=P_dc(i,t)<=P_s(i)/baseMVA];
    end
end




%% 可再生能源发电, 后面直接用场景。实际上，上限没必要但是下限必要。
C_r = [];
for t = 1:T
    C_r = [C_r; 0<=P_r(:,t)<=data_r(:,2)];
    %    C_r = [C_r; 0<=P_r(:,t)<=0];%测试的时候才用。
end




%% 系统功率平衡
C_bal = [];
for t = 1:T
    C_bal = [C_bal; sum(dem(:,t))>=sum(P_g(:,t))+sum(P_dc(:,t)-P_c(:,t))+sum(P_r(:,t)) ];
    C_bal = [C_bal; sum(dem(:,t))<=sum(P_g(:,t))+sum(P_dc(:,t)-P_c(:,t))+sum(P_r(:,t)) ];
end



%% 线性化运行成本
C_cost = [];
for i = 1:N_g
    for t = 1:T
        for seg = 1:length(k_gen(1,:))
            C_cost = [C_cost; Gcost(i,t)>=P_g(i,t)*k_gen(i,seg)+b_gen(i,seg)];
        end
    end
end







%% 约束集
C = [C_PTDF;C_gen;C_sto;C_bal;C_cost;C_r];
%% 目标函数： obj_1: 运行成本，obj_2: startup成本, obj_3: shutdown成本, obj_:储能运行成本
obj_1 = sum(sum(Gcost));
obj_2 = 0;
for i = 1:N_g
    obj_2 = obj_2+sum(b_up(i,:))*data_gcost(i,2);
end
%关机本身没有成本，为了不使第一阶段b_dowm被忽略，增加这个目标函数吧。
obj_3 = 0;
for i = 1:N_g
    obj_3 = obj_3+sum(b_down(i,:))*data_gcost(i,2)/10;
end

obj_4 = 3*sum(sum(P_dc+P_c))*baseMVA;
obj = obj_1+obj_2+obj_3+obj_4;




%% 求解
ops = sdpsettings('verbose',0,'solver','cplex');
%% 可以检验模型是否正确,但是要注意没有加可再生能源的约束。可到main中检验。
%  result = optimize(C,obj,ops)
%  value(b_state)

[model] = export(C, obj, ops);

%% 模型基本形式：  x为连续变量,y为离散变量，w为分离出来的可再生能源出力变量
%     min   cx+dy
%     s.t.  Ax+By+Ew+GD <= b+Fθ

%% 各类变量索引：要注意！
idx_bi = model.binary_variables;
idx_conti = 1:(2*N_g*T+2*N_s*T+N_s);
idx_dem = (idx_conti(end)+1):(idx_conti(end)+N_b*T);
idx_para = (idx_dem(end)+1):(idx_dem(end)+2*N_s);
idx_renew = (idx_para(end)+1):(idx_para(end)+N_r*T);
%% 各类变量数目
N_conti = length(idx_conti);
N_bi = length(idx_bi);
N_renew = length(idx_renew);
N_para = length(idx_para);
N_dem = length(idx_dem);
%% 约束矩阵/向量
A = model.Aineq(:,idx_conti);
B = model.Aineq(:,idx_bi);
E = model.Aineq(:,idx_renew);
G = model.Aineq(:,idx_dem);

b = model.bineq;
F = -model.Aineq(:,idx_para);
%% 目标函数向量
c = model.f(idx_conti)';
d = model.f(idx_bi)';

end


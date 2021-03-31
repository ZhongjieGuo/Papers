function [ A, B, E, G, b, F, c, d, N_conti,N_bi, N_renew, N_para, N_dem ] = model(data_b, data_g, data_l, data_r, data_s, data_gcost,PTDF)

baseMVA = 100;
T = 24;
eff_c  = data_s(:,2);
eff_dc = data_s(:,3);

E_lower = data_s(:,5);
ramp =  data_g(:,11); % ����  0.5*(P_max-P_min)
N_g = length(data_g(:,1));
N_l = length(data_l(:,1));
N_b = length(data_b(:,1));
N_r = length(data_r(:,1));
N_s = length(data_s(:,1));
%% �������Ի��ɱ���ϵ����k��b����ֱ���ñ���ֵ�Ϳ��ԣ��õ�������Ԫ��
[ k_gen,b_gen ] = gencost(baseMVA,N_g,data_g,data_gcost);

%% ������� ע��˳��
% ��������
P_g = sdpvar(N_g,T,'full'); % ���繦��
Gcost = sdpvar(N_g,T,'full'); % �ɱ����Ի��ĸ�������
P_c = sdpvar(N_s,T,'full'); % ���ܳ��
P_dc = sdpvar(N_s,T,'full'); % ���ܷŵ�
E_init = sdpvar(N_s,1,'full'); %���ܳ�ʼ״̬
dem = sdpvar(N_b,T,'full'); % ����

P_s = sdpvar(N_s,1,'full'); % ��������
E_s = sdpvar(N_s,1,'full');% ��������
P_r = sdpvar(N_r,T,'full');% ����Դ����
% ������
b_state = binvar(N_g,T,'full');
b_up = binvar(N_g,T,'full');
b_down = binvar(N_g,T,'full');
%% ��ʾ0ʱ�η����״̬�������Լ����һ��Ҫ�ã����Էſ���һʱ�̵�״̬���������������

b_0 = data_g(:,12);
P_g0 = b_0.*data_g(:,9)*0.5;
%% ��С��ͣ��ʱ�䡣
minT = data_g(:,13);
%% �������ߡ�


%% PTDF
C_PTDF = [];
for t = 1:T
    C_PTDF = [C_PTDF; -data_l(:,6) <=( (PTDF(:,data_g(:,1)) * P_g(:,t))...
        - (PTDF * dem(:,t))...
        + (PTDF(:,data_r(:,1)) * P_r(:,t))...
        + (PTDF(:,data_s(:,1)) * (P_dc(:,t)-P_c(:,t)))  )...
        <= data_l(:,6)];
end



%% �����
C_gen = [];
%% ��ͣ״̬��1ʱ�β�ҪԼ���ˣ�
for t = 1:1
    C_gen = [C_gen; b_0 - b_state(:,t) + b_up(:,t)>=0];  % 1ʱ������
    C_gen = [C_gen; b_state(:,t) - b_0 + b_down(:,t)>=0];  % 1ʱ��ͣ��
end
for t = 2:T
    C_gen = [C_gen; b_state(:,t-1) - b_state(:,t) + b_up(:,t)>=0]; % tʱ������
    C_gen = [C_gen; b_state(:,t) - b_state(:,t-1) + b_down(:,t)>=0]; % tʱ��ͣ��
end
%% ��С��/ͣ��ʱ��

for i = 1:N_g
    for t = 1:1 %��1ʱ�β�ҪԼ���ˣ�
        for tau = (t+1):(min(t+minT(i)-1,T)) %������t��t+1����,tʵ��������Ȼ�����
            C_gen = [C_gen; b_state(:,t) - b_0 <= b_state(:,tau)];  % 1ʱ����С����
            C_gen = [C_gen; b_0 - b_state(:,t) <= (1-b_state(:,tau))];  % 1ʱ����Сͣ��
        end
    end
    
    for t = 2:(T-1)  %���һ��ʱ�ο��Բ����Ժ��״̬
        for tau = (t+1):(min(t+minT(i)-1,T)) %������t��t+1����,tʵ��������Ȼ�����
            C_gen = [C_gen; b_state(i,t) - b_state(i,t-1) <= b_state(i,tau)];  % tʱ�κ���С����
            C_gen = [C_gen; b_state(i,t-1) - b_state(i,t) <= (1-b_state(i,tau))];  % tʱ�κ���Сͣ��
        end
    end
end

%% ����Լ��
P_g_range = data_g(:,9)-data_g(:,10);
w = 10;
for t = 1:1 %��1ʱ�β�ҪԼ���ˣ�
    C_gen = [C_gen; P_g(:,t)-P_g0 <= b_0.*P_g_range.*ramp + (1-b_0).*data_g(:,w)];
    C_gen = [C_gen; P_g0-P_g(:,t) <= b_state(:,t).*P_g_range.*ramp + (1-b_state(:,t)).*data_g(:,w)];
end
for t = 2:T
    C_gen = [C_gen; P_g(:,t)-P_g(:,t-1) <= b_state(:,t-1).*P_g_range.*ramp + (1-b_state(:,t-1)).*data_g(:,w)];
    C_gen = [C_gen; P_g(:,t-1)-P_g(:,t) <= b_state(:,t).*P_g_range.*ramp + (1-b_state(:,t)).*data_g(:,w)];
end

%% ����������
for t = 1:T
    C_gen = [C_gen; data_g(:,10).*b_state(:,t) <= P_g(:,t) <= data_g(:,9).*b_state(:,t)];
end






%% ����
C_sto = [];
%    C_sto = [C_sto; 0 <= P_s <= 0; 0<= E_s <= 0]; %��ģʱ����ȥ�����Լ�������Ե�ʱ����á�
bigM = 10^5;
for i = 1:N_s
    C_sto = [C_sto; E_s(i)*E_lower(i)/baseMVA<=E_init(i)<=E_s(i)/baseMVA;];
    for t = 1:T
        C_sto = [C_sto; E_s(i)*E_lower(i)/baseMVA<=(E_init(i)+sum(P_c(i,1:t))*eff_c(i)-sum(P_dc(i,1:t))/eff_dc(i))<=E_s(i)/baseMVA];
        C_sto = [C_sto; 0<=P_c(i,t)<=P_s(i)/baseMVA; 0<=P_dc(i,t)<=P_s(i)/baseMVA];
    end
end




%% ��������Դ����, ����ֱ���ó�����ʵ���ϣ�����û��Ҫ�������ޱ�Ҫ��
C_r = [];
for t = 1:T
    C_r = [C_r; 0<=P_r(:,t)<=data_r(:,2)];
    %    C_r = [C_r; 0<=P_r(:,t)<=0];%���Ե�ʱ����á�
end




%% ϵͳ����ƽ��
C_bal = [];
for t = 1:T
    C_bal = [C_bal; sum(dem(:,t))>=sum(P_g(:,t))+sum(P_dc(:,t)-P_c(:,t))+sum(P_r(:,t)) ];
    C_bal = [C_bal; sum(dem(:,t))<=sum(P_g(:,t))+sum(P_dc(:,t)-P_c(:,t))+sum(P_r(:,t)) ];
end



%% ���Ի����гɱ�
C_cost = [];
for i = 1:N_g
    for t = 1:T
        for seg = 1:length(k_gen(1,:))
            C_cost = [C_cost; Gcost(i,t)>=P_g(i,t)*k_gen(i,seg)+b_gen(i,seg)];
        end
    end
end







%% Լ����
C = [C_PTDF;C_gen;C_sto;C_bal;C_cost;C_r];
%% Ŀ�꺯���� obj_1: ���гɱ���obj_2: startup�ɱ�, obj_3: shutdown�ɱ�, obj_:�������гɱ�
obj_1 = sum(sum(Gcost));
obj_2 = 0;
for i = 1:N_g
    obj_2 = obj_2+sum(b_up(i,:))*data_gcost(i,2);
end
%�ػ�����û�гɱ���Ϊ�˲�ʹ��һ�׶�b_dowm�����ԣ��������Ŀ�꺯���ɡ�
obj_3 = 0;
for i = 1:N_g
    obj_3 = obj_3+sum(b_down(i,:))*data_gcost(i,2)/10;
end

obj_4 = 3*sum(sum(P_dc+P_c))*baseMVA;
obj = obj_1+obj_2+obj_3+obj_4;




%% ���
ops = sdpsettings('verbose',0,'solver','cplex');
%% ���Լ���ģ���Ƿ���ȷ,����Ҫע��û�мӿ�������Դ��Լ�����ɵ�main�м��顣
%  result = optimize(C,obj,ops)
%  value(b_state)

[model] = export(C, obj, ops);

%% ģ�ͻ�����ʽ��  xΪ��������,yΪ��ɢ������wΪ��������Ŀ�������Դ��������
%     min   cx+dy
%     s.t.  Ax+By+Ew+GD <= b+F��

%% �������������Ҫע�⣡
idx_bi = model.binary_variables;
idx_conti = 1:(2*N_g*T+2*N_s*T+N_s);
idx_dem = (idx_conti(end)+1):(idx_conti(end)+N_b*T);
idx_para = (idx_dem(end)+1):(idx_dem(end)+2*N_s);
idx_renew = (idx_para(end)+1):(idx_para(end)+N_r*T);
%% ���������Ŀ
N_conti = length(idx_conti);
N_bi = length(idx_bi);
N_renew = length(idx_renew);
N_para = length(idx_para);
N_dem = length(idx_dem);
%% Լ������/����
A = model.Aineq(:,idx_conti);
B = model.Aineq(:,idx_bi);
E = model.Aineq(:,idx_renew);
G = model.Aineq(:,idx_dem);

b = model.bineq;
F = -model.Aineq(:,idx_para);
%% Ŀ�꺯������
c = model.f(idx_conti)';
d = model.f(idx_bi)';

end


function [data_b, data_l, data_g, data_r] = IEEE33

baseMVA = 1;   %MVA
baseKV  = 12.66;  %KV
baseKA  = baseMVA/(sqrt(3)*baseKV);   %KA
baseZ   = baseKV^2/baseMVA;  %Ohms
%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
data_b = [  %% (Pd and Qd are specified in kW & kVAr here)
	1	3	0	0	0	0	1	1	0	12.66	1	1	1;
	2	1	100	60	0	0	1	1	0	12.66	1	1.1	0.9;
	3	1	90	40	0	0	1	1	0	12.66	1	1.1	0.9;
	4	1	120	80	0	0	1	1	0	12.66	1	1.1	0.9;
	5	1	60	30	0	0	1	1	0	12.66	1	1.1	0.9;
	6	1	60	20	0	0	1	1	0	12.66	1	1.1	0.9;
	7	1	200	100	0	0	1	1	0	12.66	1	1.1	0.9;
	8	1	200	100	0	0	1	1	0	12.66	1	1.1	0.9;
	9	1	60	20	0	0	1	1	0	12.66	1	1.1	0.9;
	10	1	60	20	0	0	1	1	0	12.66	1	1.1	0.9;
	11	1	45	30	0	0	1	1	0	12.66	1	1.1	0.9;
	12	1	60	35	0	0	1	1	0	12.66	1	1.1	0.9;
	13	1	60	35	0	0	1	1	0	12.66	1	1.1	0.9;
	14	1	120	80	0	0	1	1	0	12.66	1	1.1	0.9;
	15	1	60	10	0	0	1	1	0	12.66	1	1.1	0.9;
	16	1	60	20	0	0	1	1	0	12.66	1	1.1	0.9;
	17	1	60	20	0	0	1	1	0	12.66	1	1.1	0.9;
	18	1	90	40	0	0	1	1	0	12.66	1	1.1	0.9;
	19	1	90	40	0	0	1	1	0	12.66	1	1.1	0.9;
	20	1	90	40	0	0	1	1	0	12.66	1	1.1	0.9;
	21	1	90	40	0	0	1	1	0	12.66	1	1.1	0.9;
	22	1	90	40	0	0	1	1	0	12.66	1	1.1	0.9;
	23	1	90	50	0	0	1	1	0	12.66	1	1.1	0.9;
	24	1	420	200	0	0	1	1	0	12.66	1	1.1	0.9;
	25	1	420	200	0	0	1	1	0	12.66	1	1.1	0.9;
	26	1	60	25	0	0	1	1	0	12.66	1	1.1	0.9;
	27	1	60	25	0	0	1	1	0	12.66	1	1.1	0.9;
	28	1	60	20	0	0	1	1	0	12.66	1	1.1	0.9;
	29	1	120	70	0	0	1	1	0	12.66	1	1.1	0.9;
	30	1	200	600	0	0	1	1	0	12.66	1	1.1	0.9;
	31	1	150	70	0	0	1	1	0	12.66	1	1.1	0.9;
	32	1	210	100	0	0	1	1	0	12.66	1	1.1	0.9;
	33	1	60	40	0	0	1	1	0	12.66	1	1.1	0.9;
];
data_b (:,3:4) = data_b (:,3:4)/1000; %化成MW
data_b (:,3:4) = data_b (:,3:4)/baseMVA;


%% generator data   第一个slack节点看作一个发电机就行;cost跟电价相同
%	bus	P_max(MW) P_min	Qmax Qmin Cost($/MWh) % node / ^2/^1/^0
data_g = [
      1	   6	  0	  4  -4   0     80    0
      18   1.5    0   1  -1   5	    35    0
      33   1.5    0   1  -1   8	    30    0
];
data_g(:,2:5) = data_g(:,2:5)/baseMVA;
%% slack bus



%% REG   
%  bus	P_max(MW) P_min
data_r = [
     3   1.5   0
     30  1.5   0
];
data_r (:,2:3) = data_r (:,2:3)/baseMVA;




%% branch data  此问题用I的上限 
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
I_1 = 0.595/baseKA;
I_2 = 0.160/baseKA;
I_3 = 0.325/baseKA;
I_4 = 0.460/baseKA;
data_l = [  %% (r and x specified in ohms here, converted to p.u. below)
	1	2	0.0922	0.0470	0	I_1	0	0	0	0	1	-360	360;
	2	3	0.4930	0.2511	0	I_1	0	0	0	0	1	-360	360;
	3	4	0.3660	0.1864	0	I_1	0	0	0	0	1	-360	360;
	4	5	0.3811	0.1941	0	I_1	0	0	0	0	1	-360	360;
	5	6	0.8190	0.7070	0	I_1	0	0	0	0	1	-360	360;
	6	7	0.1872	0.6188	0	I_1	0	0	0	0	1	-360	360;
	7	8	0.7114	0.2351	0	I_1	0	0	0	0	1	-360	360;
	8	9	1.0300	0.7400	0	I_1	0	0	0	0	1	-360	360;
	9	10	1.0440	0.7400	0	I_1	0	0	0	0	1	-360	360;
	10	11	0.1966	0.0650	0	I_1	0	0	0	0	1	-360	360;
	11	12	0.3744	0.1238	0	I_1	0	0	0	0	1	-360	360;
	12	13	1.4680	1.1550	0	I_1	0	0	0	0	1	-360	360;
	13	14	0.5416	0.7129	0	I_1	0	0	0	0	1	-360	360;
	14	15	0.5910	0.5260	0	I_1	0	0	0	0	1	-360	360;
	15	16	0.7463	0.5450	0	I_1	0	0	0	0	1	-360	360;
	16	17	1.2890	1.7210	0	I_1	0	0	0	0	1	-360	360;
	17	18	0.7320	0.5740	0	I_1	0	0	0	0	1	-360	360;
	2	19	0.1640	0.1565	0	I_2	0	0	0	0	1	-360	360;
	19	20	1.5042	1.3554	0	I_2	0	0	0	0	1	-360	360;
	20	21	0.4095	0.4784	0	I_2	0	0	0	0	1	-360	360;
	21	22	0.7089	0.9373	0	I_2	0	0	0	0	1	-360	360;
	3	23	0.4512	0.3083	0	I_4	0	0	0	0	1	-360	360;
	23	24	0.8980	0.7091	0	I_4	0	0	0	0	1	-360	360;
	24	25	0.8960	0.7011	0	I_4	0	0	0	0	1	-360	360;
	6	26	0.2030	0.1034	0	I_3	0	0	0	0	1	-360	360;
	26	27	0.2842	0.1447	0	I_3	0	0	0	0	1	-360	360;
	27	28	1.0590	0.9337	0	I_3	0	0	0	0	1	-360	360;
	28	29	0.8042	0.7006	0	I_3	0	0	0	0	1	-360	360;
	29	30	0.5075	0.2585	0	I_3	0	0	0	0	1	-360	360;
	30	31	0.9744	0.9630	0	I_3	0	0	0	0	1	-360	360;
	31	32	0.3105	0.3619	0	I_3	0	0	0	0	1	-360	360;
	32	33	0.3410	0.5302	0	I_3	0	0	0	0	1	-360	360;
];
data_l(:,3:4) = data_l(:,3:4)/baseZ;



end


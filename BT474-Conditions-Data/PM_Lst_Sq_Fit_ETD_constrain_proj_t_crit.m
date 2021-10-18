clear all
close all
clc

%This program generates a curve using the Phenom_Model_Forward function,
%adds noise with the Add_Percent_Noise function, then uses least squares
%fitting to extract the input parameter values.
num = xlsread('GH1915_2treat_working_mod_4_27_v2.xlsx');
g_0_vector = xlsread('GH1915_2treat_g0.xlsx');
t_r_vector = xlsread('GH1915_2treat_t_r.xlsx');

t_vector = num(:,1);
[length,wid] = size(num);

trunc_points = ones(12,1);
trunc_points = trunc_points.*length;
for i = 1:12
    for j = 1:length
       if isnan(num(j,i+1)) == 1
           trunc_points(i,1) = j-1;
           break
       end
    end
end

%z_guess = [f_r,g_r,k,t_d,N_max];
z_guess = [0.5;0.025;0.005;0.02;60000];

RSS = zeros(12,1);

%To vary the t_d half-life constraint, use this as the t_d min:
%For 231:
%t_d_constraint_vector = [ones(6,1)*0.0069;ones(6,1)*0.0069;ones(6,1)*0.0075;ones(6,1)*0.0085;ones(6,1)*0.0095;ones(6,1)*0.0115;ones(6,1)*0.0125;ones(6,1)*0.015;ones(6,1)*0.023;ones(6,1)*0.035];

%For BT474:
%t_d_constraint_vector = [ones(6,1)*0.0078;ones(6,1)*0.0060;ones(6,1)*0.0060;ones(6,1)*0.0064;ones(6,1)*0.0081;ones(6,1)*0.0102;ones(6,1)*0.0124;ones(6,1)*0.0146;ones(6,1)*0.0171;ones(6,1)*0.0185];

%Constraints used for MCF7
%    z_output(m,:) = lsqnonlin(@(z) PM_Lst_Sq_Function_v15_known_t_r(z,count_vector,time_vector,num(1,m+1),g_0,t_r),z_guess,[0;0.003;0.001;0.0069;4*N_0_vector(1,m)],[1;0.1;0.01;2;120000]);

N_0_vector = num(1,2:wid);


for m = 1:12
    g_0 = g_0_vector(m);
    t_r = t_r_vector(m);
    time_vector = num(1:trunc_points(m,1),1);
    count_vector = num(1:trunc_points(m,1),m+1);
    z_output(m,:) = lsqnonlin(@(z) PM_Lst_Sq_Function_v15_known_t_r(z,count_vector,time_vector,num(1,m+1),g_0,t_r),z_guess,[0;0.003;0.001;0.0069;5000],[1;0.1;0.01;2;120000]);
%    z_output(m,:) = lsqnonlin(@(z) PM_Lst_Sq_Function_v15_known_t_r(z,count_vector,time_vector,num(1,m+1),g_0,t_r),z_guess,[0;0.003;0.001;0.0069;4*N_0_vector(1,m)],[1;0.1;0.01;2;120000]);
%    z_output(m,:) = lsqnonlin(@(z) PM_Lst_Sq_Function_v15_known_t_r(z,count_vector,time_vector,num(1,m+1),g_0,t_r),z_guess,[0;0.003;0.001;t_d_constraint_vector;5000],[1;0.1;0.01;2;120000]);
    extended_t_vector = transpose([0:2016]);
    intermediate = PM_Forward_v15(num(1,m+1),extended_t_vector,g_0,z_output(m,5),z_output(m,1),z_output(m,2),z_output(m,3),z_output(m,4),t_r);
    final_forward(:,1) = extended_t_vector;
    final_forward(:,m+1) = intermediate(:,2);
    RSS_intermediate = PM_Lst_Sq_Function_v15_known_t_r([z_output(m,1);z_output(m,2);z_output(m,3);z_output(m,4);z_output(m,5)],count_vector,time_vector,num(1,m+1),g_0,t_r);
    RSS_int_2 = RSS_intermediate.^2;
    for n = 1:trunc_points(m,1)
        RSS(m) = RSS(m)+RSS_int_2(n);
    end
end

AIC_vector = zeros(12,1);
for o = 1:12
    AIC_vector(o) = 2*5+length*log(RSS(o)/length);
end

z_output_adj = [z_output,AIC_vector];

%Find the critical time, where N = 2*N_0, for each column
numcols = wid-1;
t_crit_vector = zeros(numcols,1);
for p = 2:wid
    for q = 1:2017
        if final_forward(q,p)>=2*N_0_vector(1,p-1)
            t_crit_vector(p-1,1) = final_forward(q,1);
            break
        end
    end
end

f1 = figure;
f2 = figure;
f3 = figure;
f4 = figure;
f5 = figure;
f6 = figure;
%plot(t_vector,num(:,2),t_vector,final_forward(:,2),t_vector,num(:,3),t_vector,final_forward(:,3),t_vector,num(:,4),t_vector,final_forward(:,4),t_vector,num(:,5),t_vector,final_forward(:,5),t_vector,num(:,6),t_vector,final_forward(:,6),t_vector,num(:,7),t_vector,final_forward(:,7));
figure(f1);
plot(t_vector,num(:,8),extended_t_vector,final_forward(:,8));
figure(f2);
plot(t_vector,num(:,9),extended_t_vector,final_forward(:,9));
figure(f3);
plot(t_vector,num(:,10),extended_t_vector,final_forward(:,10));
figure(f4);
plot(t_vector,num(:,11),extended_t_vector,final_forward(:,11));
figure(f5);
plot(t_vector,num(:,12),extended_t_vector,final_forward(:,12));
figure(f6);
plot(t_vector,num(:,13),extended_t_vector,final_forward(:,13));
%plot(t_vector,num(:,3),t_vector,final_forward(:,3));
clear all
close all
clc

%This program generates a curve using the Phenom_Model_Forward function,
%adds noise with the Add_Percent_Noise function, then uses least squares
%fitting to extract the input parameter values.
num = xlsread('GH1825_working_mod_10_7_v2.xlsx');
g_0_vector = xlsread('GH1825_g_0_values.xlsx');
t_r_vector = xlsread('GH1825_t_r_values_v2.xlsx');

t_vector = num(:,1);
length = length(num);
trunc_points = ones(60,1);
trunc_points = trunc_points.*length;
for i = 1:60
    for j = 1:length
       if isnan(num(j,i+1)) == 1
           trunc_points(i,1) = j-1;
           break
       end
    end
end

%z_guess = [f_r,g_r,k,t_d,N_max];
z_guess = [0.5;0.025;0.005;0.02;60000];

RSS = zeros(60,1);

%To vary the t_d half-life constraint, use this as the t_d min:

%t_d_constraint_vector = [ones(6,1)*0.0069;ones(6,1)*0.0069;ones(6,1)*0.0075;ones(6,1)*0.0085;ones(6,1)*0.0095;ones(6,1)*0.0115;ones(6,1)*0.0125;ones(6,1)*0.015;ones(6,1)*0.023;ones(6,1)*0.035];

%For BT474:
t_d_constraint_vector = [ones(6,1)*0.0078;ones(6,1)*0.0060;ones(6,1)*0.0060;ones(6,1)*0.0064;ones(6,1)*0.0081;ones(6,1)*0.0102;ones(6,1)*0.0124;ones(6,1)*0.0146;ones(6,1)*0.0171;ones(6,1)*0.0185];



for m = 1:60
    g_0 = g_0_vector(m);
    t_r = t_r_vector(m);
    time_vector = num(1:trunc_points(m,1),1);
    count_vector = num(1:trunc_points(m,1),m+1);
%    z_output(m,:) = lsqnonlin(@(z) PM_Lst_Sq_Function_v15_known_t_r(z,count_vector,time_vector,num(1,m+1),g_0,t_r),z_guess,[0;0.003;0.001;0.0069;5000],[1;0.1;0.01;2;120000]);
    z_output(m,:) = lsqnonlin(@(z) PM_Lst_Sq_Function_v15_known_t_r(z,count_vector,time_vector,num(1,m+1),g_0,t_r),z_guess,[0;0.003;0.001;t_d_constraint_vector;5000],[1;0.1;0.01;2;120000]);
    intermediate = PM_Forward_v15(num(1,m+1),time_vector,g_0,z_output(m,5),z_output(m,1),z_output(m,2),z_output(m,3),z_output(m,4),t_r);
    final_forward(:,1) = t_vector;
    final_forward(1:trunc_points(m,1),m+1) = intermediate(:,2);
    RSS_intermediate = PM_Lst_Sq_Function_v15_known_t_r([z_output(m,1);z_output(m,2);z_output(m,3);z_output(m,4);z_output(m,5)],count_vector,time_vector,num(1,m+1),g_0,t_r);
    RSS_int_2 = RSS_intermediate.^2;
    for n = 1:trunc_points(m,1)
        RSS(m) = RSS(m)+RSS_int_2(n);
    end
end

AIC_vector = zeros(60,1);
for o = 1:60
    AIC_vector(o) = 2*5+length*log(RSS(o)/length);
end

z_output_adj = [z_output,AIC_vector];

f1 = figure;
f2 = figure;
f3 = figure;
f4 = figure;
f5 = figure;
f6 = figure;
%plot(t_vector,num(:,2),t_vector,final_forward(:,2),t_vector,num(:,3),t_vector,final_forward(:,3),t_vector,num(:,4),t_vector,final_forward(:,4),t_vector,num(:,5),t_vector,final_forward(:,5),t_vector,num(:,6),t_vector,final_forward(:,6),t_vector,num(:,7),t_vector,final_forward(:,7));
figure(f1);
plot(t_vector,num(:,2),t_vector,final_forward(:,2));
figure(f2);
plot(t_vector,num(:,3),t_vector,final_forward(:,3));
figure(f3);
plot(t_vector,num(:,4),t_vector,final_forward(:,4));
figure(f4);
plot(t_vector,num(:,5),t_vector,final_forward(:,5));
figure(f5);
plot(t_vector,num(:,6),t_vector,final_forward(:,6));
figure(f6);
plot(t_vector,num(:,7),t_vector,final_forward(:,7));
%plot(t_vector,num(:,3),t_vector,final_forward(:,3));
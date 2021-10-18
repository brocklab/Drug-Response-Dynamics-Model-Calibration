clear all
close all
clc

%This program generates a curve using the Phenom_Model_Forward function,
%adds noise with the Add_Percent_Noise function, then uses least squares
%fitting to extract the input parameter values.
num = xlsread('GH1910_16dInt_working_modified_9_11.xlsx');
%g_0_vector = xlsread('g_0_1910_16dInt.xlsx');
t_r_vector = xlsread('GH1910_16dInt_t_r_values.xlsx');

t_vector = num(:,1);
length = length(num);
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

z_guess = [0.5;0.025;0.005;60000];

RSS = zeros(12,1);

for m = 1:12
    t_r = t_r_vector(m);
    time_vector = num(1:trunc_points(m,1),1);
    count_vector = num(1:trunc_points(m,1),m+1);
    z_output(m,:) = lsqnonlin(@(z) PM_Lst_Sq_Function_v17_known_t_r(z,count_vector,time_vector,num(1,m+1),t_r),z_guess,[0;0;0;0],[1;0.1;0.5;200000]);
    intermediate = PM_Forward_v17(num(1,m+1),time_vector,z_output(m,4),z_output(m,1),z_output(m,2),z_output(m,3),t_r);
    final_forward(:,1) = t_vector;
    final_forward(1:trunc_points(m,1),m+1) = intermediate(:,2);
    RSS_intermediate = PM_Lst_Sq_Function_v17_known_t_r([z_output(m,1);z_output(m,2);z_output(m,3);z_output(m,4)],num(:,m+1),num(:,1),num(1,m+1),t_r);
    RSS_int_2 = RSS_intermediate.^2;
    for n = 1:trunc_points(m,1)
        RSS(m) = RSS(m)+RSS_int_2(n);
    end
end

AIC_vector = zeros(12,1);
for o = 1:12
    AIC_vector(o) = 2*4+length*log(RSS(o)/length);
end

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
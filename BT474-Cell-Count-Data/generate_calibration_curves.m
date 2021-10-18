clear all
close all
clc

%Read parameter values in

num = xlsread('GH1818_75nM_param_values.xlsx');

f_r_vec=num(:,1);
g_r_vec=num(:,2);
k_d_vec=num(:,3);
t_d_vec=num(:,4);
t_r_vec=num(:,6);
g_0_vec=num(:,7);
N_max_vec=num(:,5);
N_0_vec=num(:,8);

param_ex_1=num(2:6,:);
param_ex_2=[num(1,:);num(3:6,:)];
param_ex_3=[num(1:2,:);num(4:6,:)];
param_ex_4=[num(1:3,:);num(5:6,:)];
param_ex_5=[num(1:4,:);num(6,:)];
param_ex_6=num(1:5,:);

for i=1:6
    f_r_ave_vec(i,1)=(f_r_vec(1)+f_r_vec(2)+f_r_vec(3)+f_r_vec(4)+f_r_vec(5)+f_r_vec(6)-f_r_vec(i))/5;
    g_r_ave_vec(i,1)=(g_r_vec(1)+g_r_vec(2)+g_r_vec(3)+g_r_vec(4)+g_r_vec(5)+g_r_vec(6)-g_r_vec(i))/5;
    g_0_ave_vec(i,1)=(g_0_vec(1)+g_0_vec(2)+g_0_vec(3)+g_0_vec(4)+g_0_vec(5)+g_0_vec(6)-g_0_vec(i))/5;
    k_d_ave_vec(i,1)=(k_d_vec(1)+k_d_vec(2)+k_d_vec(3)+k_d_vec(4)+k_d_vec(5)+k_d_vec(6)-k_d_vec(i))/5;
    N_0_ave_vec(i,1)=(N_0_vec(1)+N_0_vec(2)+N_0_vec(3)+N_0_vec(4)+N_0_vec(5)+N_0_vec(6)-N_0_vec(i))/5;
    N_max_ave_vec(i,1)=(N_max_vec(1)+N_max_vec(2)+N_max_vec(3)+N_max_vec(4)+N_max_vec(5)+N_max_vec(6)-N_max_vec(i))/5;
    t_d_ave_vec(i,1)=(t_d_vec(1)+t_d_vec(2)+t_d_vec(3)+t_d_vec(4)+t_d_vec(5)+t_d_vec(6)-t_d_vec(i))/5;
end

for j=1:8
    param_std_devs(1,j)=std(param_ex_1(:,j));
    param_std_devs(2,j)=std(param_ex_2(:,j));
    param_std_devs(3,j)=std(param_ex_3(:,j));
    param_std_devs(4,j)=std(param_ex_4(:,j));
    param_std_devs(5,j)=std(param_ex_5(:,j));
    param_std_devs(6,j)=std(param_ex_6(:,j));
end

f_r_std_vec=param_std_devs(:,1);
g_r_std_vec=param_std_devs(:,2);
k_d_std_vec=param_std_devs(:,3);
t_d_std_vec=param_std_devs(:,4);
t_r_std_vec=param_std_devs(:,6);
g_0_std_vec=param_std_devs(:,7);
N_max_std_vec=param_std_devs(:,5);
N_0_std_vec=param_std_devs(:,8);

%Now use these leave-one-out average and standard deviation values to
%generate 1000 data sets per replicate - 6000 total
m=0;
for k=1:6
    for l=1:1000
        m=m+1;
        parameters(m,1) = random('Normal',f_r_ave_vec(k,1),f_r_std_vec(k,1));
           if parameters(m,1)<0
       parameters(m,1)=0;
           end
       parameters(m,2) = random('Normal',g_r_ave_vec(k,1),g_r_std_vec(k,1));
           if parameters(m,2)<0
       parameters(m,2)=0;
           end
       parameters(m,3) = random('Normal',k_d_ave_vec(k,1),k_d_std_vec(k,1));
           if parameters(m,3)<0
       parameters(m,3)=0;
           end
       parameters(m,4) = random('Normal',t_d_ave_vec(k,1),t_d_std_vec(k,1));
           if parameters(m,4)<0
       parameters(m,4)=0;
           end
       parameters(m,5) = random('Normal',g_0_ave_vec(k,1),g_0_std_vec(k,1));
           if parameters(m,5)<0
       parameters(m,5)=0;
           end
       parameters(m,6) = random('Normal',N_max_ave_vec(k,1),N_max_std_vec(k,1));
           if parameters(m,6)<0
       parameters(m,6)=0;
           end
    end
end

%load the time vector

t_vec= xlsread('GH1818_75nM_t_vec.xlsx');

n=0;
for o=1:6
    for p=1:1000
        n=n+1;
   curve=PM_Forward_v15(N_0_vec(1,1),t_vec,parameters(n,5),parameters(n,6),parameters(n,1),parameters(n,2),parameters(n,3),parameters(n,4),t_r_vec(o,1)); 
    output_curves(:,n)=curve(:,2);
    end
end

%strip infinite values

for q=1:6000
   for r=1:206
       if isnan(output_curves(r,q))==0
           final_output(r,q)=output_curves(r,q);
       end
   end
end

final_output_t=transpose(final_output);

clear all
close all
clc

%This program is intended to randomly generate 100 curves based on the
%phenomenological model with randomly generated parameter values.  This
%version will have no t_r, proliferation for the resistant cells will start
%immediately.

t_vector = [0:1:1000];
N_init = 5000;
N_max = 50000;
f_r_mean = 0.50;
f_r_stdev = 0.15;
g_r_mean = 0.025;
g_r_stdev = 0.005;
k_mean = 0.005;
k_stdev = 0.0015;
t_d_mean = 0.02;
t_d_stdev = 0.005;
t_r_mean = 200;
t_r_stdev = 30;

%generate 100 sets of parameters
for i = 1:100
   parameters(i,1) = random('Normal',f_r_mean,f_r_stdev); 
   parameters(i,2) = random('Normal',g_r_mean,g_r_stdev); 
   parameters(i,3) = random('Normal',k_mean,k_stdev);
   parameters(i,4) = random('Normal',t_d_mean,t_d_stdev);
   parameters(i,5) = random('Normal',t_r_mean,t_r_stdev);
   
   if parameters(i,1)<0
       parameters(i,1)=0;
   end
   if parameters(i,2)<0
       parameters(i,2)=0;
   end
   if parameters(i,3)<0
       parameters(i,3)=0;
   end
   if parameters(i,4)<0
       parameters(i,4)=0;
   end
   if parameters(i,5)<0
       parameters(i,5)=0;
   end
end


%Generate exact cell number curves from parameters.
cell_number = zeros(1001,101);
cell_number(:,1) = t_vector;
for j = 1:100
    f_r = parameters(j,1);
    g_r = parameters(j,2);
    k_D = parameters(j,3);
    t_d = parameters(j,4);
    t_r = parameters(j,5);
    intermediate = PM_Forward_v13(N_init,t_vector,N_max,f_r,g_r,k_D,t_d,t_r);
    cell_number(:,j+1) = intermediate(:,2);
end

%Add noise to the curves.
noisy_data(:,1) = t_vector;
for l = 1:100
    intermediate = Add_Percent_Noise([cell_number(:,1),cell_number(:,l+1)]);
    noisy_data(:,l+1) = intermediate(:,2);
end


plot(cell_number(:,1),cell_number(:,2),cell_number(:,1),cell_number(:,3),cell_number(:,1),cell_number(:,4),cell_number(:,1),noisy_data(:,2),cell_number(:,1),noisy_data(:,3),cell_number(:,1),noisy_data(:,4))
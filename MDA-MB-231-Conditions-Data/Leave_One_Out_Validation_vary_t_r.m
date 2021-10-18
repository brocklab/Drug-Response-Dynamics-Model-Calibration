clear all
close all
clc

%Read in the data files, including the cell number curves, the conditions
%file, and the calibration output which has parameter values.

data = xlsread('GH1818_75nM_working_mod_6_15.xlsx');
conditions = xlsread('GH1818_75nM_conditions.xlsx');
calibration_filename = 'GH1818_75nM_model_calibration_output.xlsx';
% model_1_parameters = xlsread(calibration_filename,2);
% model_2_parameters = xlsread(calibration_filename,4);
% model_3_parameters = xlsread(calibration_filename,6);
% model_1_death_only_parameters = xlsread(calibration_filename,8);
% model_2_death_only_parameters = xlsread(calibration_filename,10);
% model_3_death_only_parameters = xlsread(calibration_filename,12);
% posttreat_time_vector = xlsread(calibration_filename,16);

% Use readmatrix for more compliant import
model_1_parameters = readmatrix(calibration_filename, 'UseExcel',1,'Sheet',2);
model_2_parameters = readmatrix(calibration_filename, 'UseExcel',1,'Sheet',4);
model_3_parameters = readmatrix(calibration_filename, 'UseExcel',1,'Sheet',6);
model_1_death_only_parameters = readmatrix(calibration_filename, 'UseExcel',1,'Sheet',4);
model_2_death_only_parameters = readmatrix(calibration_filename, 'UseExcel',1,'Sheet',10);
model_3_death_only_parameters = readmatrix(calibration_filename, 'UseExcel',1,'Sheet',12);
posttreat_time_vector = readmatrix(calibration_filename, 'UseExcel',1,'Sheet',16);


%Specify the name you want to use for the output file

output_filename = 'GH1818_75nM_loo_validation_v1_output.xlsx';

%Extract key values from the input files
num_rep = conditions(1,1);
num_time_points = length(data);
recovery_boolean = conditions(:,2);
t_r_vec = conditions(:,3);
time_of_treat = conditions(10,1);

%Use time of treatment to divide data into pre and post treatment sections;
%the pretreatment section will include up to 48 hours of data immediately
%preceding the start of the final drug exposure, and is used to calculate
%g_0, the pre-treatment net growth rate for the replicate.

pretreat_index = 0;
for i = 1:num_time_points
    if data(i,1)>=(time_of_treat-48)
        pretreat_index = i;
        break
    end
end
treatment_index = 0;
for i = 1:num_time_points
    if data(i,1)>=time_of_treat
        treatment_index = i;
        break
    end
end

pretreat_data = data(pretreat_index:(treatment_index-1),:);
posttreat_data = data(treatment_index:num_time_points,:);
N_0_vector = posttreat_data(1,2:(num_rep+1));

pretreat_t_vector = pretreat_data(:,1);
pretreat_length = length(pretreat_data);
pretreat_z_guess = [0.025;60000];

%Use the pretreatment data to find the net growth rate before the
%treatment.

for m = 1:num_rep
    pretreat_logistic_parameters(m,:) = lsqnonlin(@(z) logistic_growth_lst_sq_fn(z,pretreat_data(:,m+1),pretreat_data(:,1)),pretreat_z_guess,[-0.1;5000],[0.1;200000]);
end

g_0_vector = pretreat_logistic_parameters(:,1);

%Replicates may vary in the time point at which they are truncated; we need
%to find the length of each data vector.

trunc_points = ones(num_rep,1);
posttreat_length=length(posttreat_data);
trunc_points = trunc_points.*posttreat_length;
for i = 1:num_rep
    for j = 1:posttreat_length
       if isnan(posttreat_data(j,i+1)) == 1
           trunc_points(i,1) = j-1;
           break
       end
    end
end

%First we will use the recovery boolean to create separate parameter rables
%for the replicates which do recover and the replicates which do not.

recovery_reps = 0;
dying_reps = 0;

for n=1:num_rep
    if recovery_boolean(n,1) == 1
        recovery_reps = recovery_reps+1;
    else
        dying_reps = dying_reps+1;
    end
end

model_1_recovering_parameters = zeros(recovery_reps,10);
model_1_dying_parameters = zeros(dying_reps,6);
model_2_recovering_parameters = zeros(recovery_reps,10);
model_2_dying_parameters = zeros(dying_reps,6);
model_3_recovering_parameters = zeros(recovery_reps,9);
model_3_dying_parameters = zeros(dying_reps,5);

recovering_index=1;
dying_index=1;

for n=1:num_rep
    if recovery_boolean(n,1) == 1
        model_1_recovering_parameters(recovering_index,:)=[model_1_parameters(n,:),n,t_r_vec(n,1),g_0_vector(n,1),trunc_points(n,1),N_0_vector(1,n)];
        model_2_recovering_parameters(recovering_index,:)=[model_2_parameters(n,:),n,t_r_vec(n,1),g_0_vector(n,1),trunc_points(n,1),N_0_vector(1,n)];
        model_3_recovering_parameters(recovering_index,:)=[model_3_parameters(n,:),n,t_r_vec(n,1),g_0_vector(n,1),trunc_points(n,1),N_0_vector(1,n)];
        recovering_index=recovering_index+1;
    else
        model_1_dying_parameters(dying_index,:)=[model_1_death_only_parameters(n,:),n,g_0_vector(n,1),trunc_points(n,1),N_0_vector(1,n)];
        model_2_dying_parameters(dying_index,:)=[model_2_death_only_parameters(n,:),n,g_0_vector(n,1),trunc_points(n,1),N_0_vector(1,n)];
        model_3_dying_parameters(dying_index,:)=[model_3_death_only_parameters(n,:),n,g_0_vector(n,1),trunc_points(n,1),N_0_vector(1,n)];
        dying_index=dying_index+1;
    end
end

%Now we want to find averages and standard deviations for each parameter,
%in a leave-one-out pattern; that is, we want an average value for each
%parameter, for each replicate, where that average is calculated from the
%values of the other replicates, excluding the current replicate, and the
%same for standard deviations.  So if there are 6 replicates, the average
%and standard deviation used to validate replicate 1 will be based on the
%average and standard deviation of each parameter for replicates 2 through
%5.

%To do this, we construct a matrix for each parameter, such that the ith
%column of the matrix consists of the values for that parameter of all
%other replicates, excluding the ith replicate.  We will need one matrix
%per parameter, per model:

model_1_recovering_loo_f_r_values = zeros((recovery_reps-1),recovery_reps);
model_1_recovering_loo_g_r_values = zeros((recovery_reps-1),recovery_reps);
model_1_recovering_loo_k_d_values = zeros((recovery_reps-1),recovery_reps);
model_1_recovering_loo_N_max_values = zeros((recovery_reps-1),recovery_reps);
model_1_recovering_loo_t_d_values = zeros((recovery_reps-1),recovery_reps);
model_1_recovering_loo_t_r_values = zeros((recovery_reps-1),recovery_reps);

model_2_recovering_loo_f_r_values = zeros((recovery_reps-1),recovery_reps);
model_2_recovering_loo_g_r_values = zeros((recovery_reps-1),recovery_reps);
model_2_recovering_loo_k_d_values = zeros((recovery_reps-1),recovery_reps);
model_2_recovering_loo_N_max_values = zeros((recovery_reps-1),recovery_reps);
model_2_recovering_loo_t_d_values = zeros((recovery_reps-1),recovery_reps);
model_2_recovering_loo_t_r_values = zeros((recovery_reps-1),recovery_reps);

model_3_recovering_loo_f_r_values = zeros((recovery_reps-1),recovery_reps);
model_3_recovering_loo_g_r_values = zeros((recovery_reps-1),recovery_reps);
model_3_recovering_loo_k_d_values = zeros((recovery_reps-1),recovery_reps);
model_3_recovering_loo_N_max_values = zeros((recovery_reps-1),recovery_reps);
model_3_recovering_loo_t_r_values = zeros((recovery_reps-1),recovery_reps);

model_1_dying_loo_k_d_values = zeros((dying_reps-1),dying_reps);
model_1_dying_loo_t_d_values = zeros((dying_reps-1),dying_reps);

model_2_dying_loo_k_d_values = zeros((dying_reps-1),dying_reps);
model_2_dying_loo_t_d_values = zeros((dying_reps-1),dying_reps);

model_3_dying_loo_k_d_values = zeros((dying_reps-1),dying_reps);

%Now we will use a loop to fill in the values of the non-excluded replicate
%parameters in each column.


for n=1:recovery_reps
    row_index = 1;
    for m=1:recovery_reps
        if n ~= m
            model_1_recovering_loo_f_r_values(row_index,n) = model_1_recovering_parameters(m,1);
            model_1_recovering_loo_g_r_values(row_index,n) = model_1_recovering_parameters(m,2);
            model_1_recovering_loo_k_d_values(row_index,n) = model_1_recovering_parameters(m,3);
            model_1_recovering_loo_N_max_values(row_index,n) = model_1_recovering_parameters(m,4);
            model_1_recovering_loo_t_d_values(row_index,n) = model_1_recovering_parameters(m,5);
            model_1_recovering_loo_t_r_values(row_index,n) = model_1_recovering_parameters(m,7);

            model_2_recovering_loo_f_r_values(row_index,n) = model_2_recovering_parameters(m,1);
            model_2_recovering_loo_g_r_values(row_index,n) = model_2_recovering_parameters(m,2);
            model_2_recovering_loo_k_d_values(row_index,n) = model_2_recovering_parameters(m,3);
            model_2_recovering_loo_N_max_values(row_index,n) = model_2_recovering_parameters(m,4);
            model_2_recovering_loo_t_d_values(row_index,n) = model_2_recovering_parameters(m,5);
            model_2_recovering_loo_t_r_values(row_index,n) = model_2_recovering_parameters(m,7);

            model_3_recovering_loo_f_r_values(row_index,n) = model_3_recovering_parameters(m,1);
            model_3_recovering_loo_g_r_values(row_index,n) = model_3_recovering_parameters(m,2);
            model_3_recovering_loo_k_d_values(row_index,n) = model_3_recovering_parameters(m,3);
            model_3_recovering_loo_N_max_values(row_index,n) = model_3_recovering_parameters(m,4);
            model_3_recovering_loo_t_r_values(row_index,n) = model_3_recovering_parameters(m,6);
            row_index=row_index+1;
        end
    end
end

for n=1:dying_reps
    row_index = 1;
    for m=1:dying_reps
        if n~=m
            model_1_dying_loo_k_d_values(row_index,n) = model_1_dying_parameters(m,1);
            model_1_dying_loo_t_d_values(row_index,n) = model_1_dying_parameters(m,2);

            model_2_dying_loo_k_d_values(row_index,n) = model_2_dying_parameters(m,1);
            model_2_dying_loo_t_d_values(row_index,n) = model_2_dying_parameters(m,2);

            model_3_dying_loo_k_d_values(row_index,n) = model_3_dying_parameters(m,1);
            row_index=row_index+1;
        end
    end
end

%Now we want to use each vector to find an average and standard deviation
%for that parameter population.  Matrices to contain these parameter
%statistics:

model_1_parameter_averages_recovering=zeros(recovery_reps,7);
model_1_parameter_averages_dying=zeros(dying_reps,3);
model_2_parameter_averages_recovering=zeros(recovery_reps,7);
model_2_parameter_averages_dying=zeros(dying_reps,3);
model_3_parameter_averages_recovering=zeros(recovery_reps,7);
model_3_parameter_averages_dying=zeros(dying_reps,2);

model_1_parameter_stddev_recovering=zeros(recovery_reps,7);
model_1_parameter_stddev_dying=zeros(dying_reps,3);
model_2_parameter_stddev_recovering=zeros(recovery_reps,7);
model_2_parameter_stddev_dying=zeros(dying_reps,3);
model_3_parameter_stddev_recovering=zeros(recovery_reps,7);
model_3_parameter_stddev_dying=zeros(dying_reps,2);

%And now use a loop to calculate those values:

for n=1:recovery_reps
    model_1_parameter_averages_recovering(n,1)=mean(model_1_recovering_loo_f_r_values(:,n));
    model_1_parameter_averages_recovering(n,2)=mean(model_1_recovering_loo_g_r_values(:,n));
    model_1_parameter_averages_recovering(n,3)=mean(model_1_recovering_loo_k_d_values(:,n));
    model_1_parameter_averages_recovering(n,4)=mean(model_1_recovering_loo_N_max_values(:,n));
    model_1_parameter_averages_recovering(n,5)=mean(model_1_recovering_loo_t_d_values(:,n));
    model_1_parameter_averages_recovering(n,6)=model_1_recovering_parameters(n,6);
    model_1_parameter_averages_recovering(n,7)=mean(model_1_recovering_loo_t_r_values(:,n));
    
    model_2_parameter_averages_recovering(n,1)=mean(model_2_recovering_loo_f_r_values(:,n));
    model_2_parameter_averages_recovering(n,2)=mean(model_2_recovering_loo_g_r_values(:,n));
    model_2_parameter_averages_recovering(n,3)=mean(model_2_recovering_loo_k_d_values(:,n));
    model_2_parameter_averages_recovering(n,4)=mean(model_2_recovering_loo_N_max_values(:,n));
    model_2_parameter_averages_recovering(n,5)=mean(model_2_recovering_loo_t_d_values(:,n));
    model_2_parameter_averages_recovering(n,6)=model_2_recovering_parameters(n,6);
    model_2_parameter_averages_recovering(n,7)=mean(model_2_recovering_loo_t_r_values(:,n));
    
    model_3_parameter_averages_recovering(n,1)=mean(model_3_recovering_loo_f_r_values(:,n));
    model_3_parameter_averages_recovering(n,2)=mean(model_3_recovering_loo_g_r_values(:,n));
    model_3_parameter_averages_recovering(n,3)=mean(model_3_recovering_loo_k_d_values(:,n));
    model_3_parameter_averages_recovering(n,4)=mean(model_3_recovering_loo_N_max_values(:,n));
    model_3_parameter_averages_recovering(n,5)=model_3_recovering_parameters(n,5);
    model_2_parameter_averages_recovering(n,6)=mean(model_3_recovering_loo_t_r_values(:,n));
    
    model_1_parameter_stddev_recovering(n,1)=std(model_1_recovering_loo_f_r_values(:,n));
    model_1_parameter_stddev_recovering(n,2)=std(model_1_recovering_loo_g_r_values(:,n));
    model_1_parameter_stddev_recovering(n,3)=std(model_1_recovering_loo_k_d_values(:,n));
    model_1_parameter_stddev_recovering(n,4)=std(model_1_recovering_loo_N_max_values(:,n));
    model_1_parameter_stddev_recovering(n,5)=std(model_1_recovering_loo_t_d_values(:,n));
    model_1_parameter_stddev_recovering(n,6)=model_1_recovering_parameters(n,6);
    model_1_parameter_stddev_recovering(n,7)=std(model_1_recovering_loo_t_r_values(:,n));
    
    model_2_parameter_stddev_recovering(n,1)=std(model_2_recovering_loo_f_r_values(:,n));
    model_2_parameter_stddev_recovering(n,2)=std(model_2_recovering_loo_g_r_values(:,n));
    model_2_parameter_stddev_recovering(n,3)=std(model_2_recovering_loo_k_d_values(:,n));
    model_2_parameter_stddev_recovering(n,4)=std(model_2_recovering_loo_N_max_values(:,n));
    model_2_parameter_stddev_recovering(n,5)=std(model_2_recovering_loo_t_d_values(:,n));
    model_2_parameter_stddev_recovering(n,6)=model_2_recovering_parameters(n,6);
    model_2_parameter_stddev_recovering(n,7)=std(model_2_recovering_loo_t_r_values(:,n));
    
    model_3_parameter_stddev_recovering(n,1)=std(model_3_recovering_loo_f_r_values(:,n));
    model_3_parameter_stddev_recovering(n,2)=std(model_3_recovering_loo_g_r_values(:,n));
    model_3_parameter_stddev_recovering(n,3)=std(model_3_recovering_loo_k_d_values(:,n));
    model_3_parameter_stddev_recovering(n,4)=std(model_3_recovering_loo_N_max_values(:,n));
    model_3_parameter_stddev_recovering(n,5)=model_1_recovering_parameters(n,5);
    model_3_parameter_stddev_recovering(n,6)=std(model_3_recovering_loo_t_r_values(:,n));
end

for n=1:dying_reps
    model_1_parameter_averages_dying(n,1)=mean(model_1_dying_loo_k_d_values(:,n));
    model_1_parameter_averages_dying(n,2)=mean(model_1_dying_loo_t_d_values(:,n));
    model_1_parameter_averages_dying(n,3)=model_1_dying_parameters(n,3);
    
    model_2_parameter_averages_dying(n,1)=mean(model_2_dying_loo_k_d_values(:,n));
    model_2_parameter_averages_dying(n,2)=mean(model_2_dying_loo_t_d_values(:,n));
    model_2_parameter_averages_dying(n,3)=model_2_dying_parameters(n,3);
    
    model_3_parameter_averages_dying(n,1)=mean(model_3_dying_loo_k_d_values(:,n));
    model_3_parameter_averages_dying(n,2)=model_3_dying_parameters(n,2);
    
    model_1_parameter_stddev_dying(n,1)=std(model_1_dying_loo_k_d_values(:,n));
    model_1_parameter_stddev_dying(n,2)=std(model_1_dying_loo_t_d_values(:,n));
    model_1_parameter_stddev_dying(n,3)=model_1_dying_parameters(n,3);
    
    model_2_parameter_stddev_dying(n,1)=std(model_2_dying_loo_k_d_values(:,n));
    model_2_parameter_stddev_dying(n,2)=std(model_2_dying_loo_t_d_values(:,n));
    model_2_parameter_stddev_dying(n,3)=model_2_dying_parameters(n,3);
    
    model_3_parameter_stddev_dying(n,1)=std(model_3_dying_loo_k_d_values(:,n));
    model_3_parameter_stddev_dying(n,2)=model_3_dying_parameters(n,2);
end

%Now use these leave-one-out average and standard deviation values to
%generate 1000 data sets per replicate, per model 

m=0;
for k=1:recovery_reps
    for l=1:1000
        m=m+1;
        model_1_recovering_generated_parameters(m,1) = random('Normal',model_1_parameter_averages_recovering(k,1),model_1_parameter_stddev_recovering(k,1));
           if model_1_recovering_generated_parameters(m,1)<0
       model_1_recovering_generated_parameters(m,1)=0;
           end
       model_1_recovering_generated_parameters(m,2) = random('Normal',model_1_parameter_averages_recovering(k,2),model_1_parameter_stddev_recovering(k,2));
           if model_1_recovering_generated_parameters(m,2)<0
       model_1_recovering_generated_parameters(m,2)=0;
           end
       model_1_recovering_generated_parameters(m,3) = random('Normal',model_1_parameter_averages_recovering(k,3),model_1_parameter_stddev_recovering(k,3));
           if model_1_recovering_generated_parameters(m,3)<0
       model_1_recovering_generated_parameters(m,3)=0;
           end
       model_1_recovering_generated_parameters(m,4) = random('Normal',model_1_parameter_averages_recovering(k,4),model_1_parameter_stddev_recovering(k,4));
           if model_1_recovering_generated_parameters(m,4)<5000
       model_1_recovering_generated_parameters(m,4)=5000;
           end
       model_1_recovering_generated_parameters(m,5) = random('Normal',model_1_parameter_averages_recovering(k,5),model_1_parameter_stddev_recovering(k,5));
           if model_1_recovering_generated_parameters(m,5)<0
       model_1_recovering_generated_parameters(m,5)=0;
           end
       model_1_recovering_generated_parameters(m,6) = random('Normal',model_1_parameter_averages_recovering(k,7),model_1_parameter_stddev_recovering(k,7));
           if model_1_recovering_generated_parameters(m,6)<0
       model_1_recovering_generated_parameters(m,6)=0;
           end    
           
           
        model_2_recovering_generated_parameters(m,1) = random('Normal',model_2_parameter_averages_recovering(k,1),model_2_parameter_stddev_recovering(k,1));
           if model_2_recovering_generated_parameters(m,1)<0
       model_2_recovering_generated_parameters(m,1)=0;
           end
       model_2_recovering_generated_parameters(m,2) = random('Normal',model_2_parameter_averages_recovering(k,2),model_2_parameter_stddev_recovering(k,2));
           if model_2_recovering_generated_parameters(m,2)<0
       model_2_recovering_generated_parameters(m,2)=0;
           end
       model_2_recovering_generated_parameters(m,3) = random('Normal',model_2_parameter_averages_recovering(k,3),model_2_parameter_stddev_recovering(k,3));
           if model_2_recovering_generated_parameters(m,3)<0
       model_2_recovering_generated_parameters(m,3)=0;
           end
       model_2_recovering_generated_parameters(m,4) = random('Normal',model_2_parameter_averages_recovering(k,4),model_2_parameter_stddev_recovering(k,4));
           if model_2_recovering_generated_parameters(m,4)<5000
       model_2_recovering_generated_parameters(m,4)=5000;
           end
       model_2_recovering_generated_parameters(m,5) = random('Normal',model_2_parameter_averages_recovering(k,5),model_2_parameter_stddev_recovering(k,5));
           if model_2_recovering_generated_parameters(m,5)<0
       model_2_recovering_generated_parameters(m,5)=0;
           end
       model_2_recovering_generated_parameters(m,6) = random('Normal',model_2_parameter_averages_recovering(k,7),model_2_parameter_stddev_recovering(k,7));
           if model_2_recovering_generated_parameters(m,6)<0
       model_2_recovering_generated_parameters(m,6)=0;
           end
           
       model_3_recovering_generated_parameters(m,1) = random('Normal',model_3_parameter_averages_recovering(k,1),model_3_parameter_stddev_recovering(k,1));
           if model_3_recovering_generated_parameters(m,1)<0
       model_3_recovering_generated_parameters(m,1)=0;
           end
       model_3_recovering_generated_parameters(m,2) = random('Normal',model_3_parameter_averages_recovering(k,2),model_3_parameter_stddev_recovering(k,2));
           if model_3_recovering_generated_parameters(m,2)<0
       model_3_recovering_generated_parameters(m,2)=0;
           end
       model_3_recovering_generated_parameters(m,3) = random('Normal',model_3_parameter_averages_recovering(k,3),model_3_parameter_stddev_recovering(k,3));
           if model_3_recovering_generated_parameters(m,3)<0
       model_3_recovering_generated_parameters(m,3)=0;
           end
       model_3_recovering_generated_parameters(m,4) = random('Normal',model_3_parameter_averages_recovering(k,4),model_3_parameter_stddev_recovering(k,4));
           if model_3_recovering_generated_parameters(m,4)<5000
       model_3_recovering_generated_parameters(m,4)=5000;
           end
       model_3_recovering_generated_parameters(m,5) = random('Normal',model_3_parameter_averages_recovering(k,6),model_3_parameter_stddev_recovering(k,6));
           if model_3_recovering_generated_parameters(m,5)<5000
       model_3_recovering_generated_parameters(m,5)=5000;
           end
    end
end

m=0;
for k=1:dying_reps
    for l=1:1000
        m=m+1;
        model_1_dying_generated_parameters(m,1) = random('Normal',model_1_parameter_averages_dying(k,1),model_1_parameter_stddev_dying(k,1));
           if model_1_dying_generated_parameters(m,1)<0
       model_1_dying_generated_parameters(m,1)=0;
           end
       model_1_dying_generated_parameters(m,2) = random('Normal',model_1_parameter_averages_dying(k,2),model_1_parameter_stddev_dying(k,2));
           if model_1_dying_generated_parameters(m,2)<0
       model_1_dying_generated_parameters(m,2)=0;
           end
           
        model_2_dying_generated_parameters(m,1) = random('Normal',model_2_parameter_averages_dying(k,1),model_2_parameter_stddev_dying(k,1));
           if model_2_dying_generated_parameters(m,1)<0
       model_2_dying_generated_parameters(m,1)=0;
           end
       model_2_dying_generated_parameters(m,2) = random('Normal',model_2_parameter_averages_dying(k,2),model_2_parameter_stddev_dying(k,2));
           if model_2_dying_generated_parameters(m,2)<0
       model_2_dying_generated_parameters(m,2)=0;
           end
           
        model_3_dying_generated_parameters(m,1) = random('Normal',model_3_parameter_averages_dying(k,1),model_3_parameter_stddev_dying(k,1));
           if model_3_dying_generated_parameters(m,1)<0
       model_3_dying_generated_parameters(m,1)=0;
           end
    end
end

%Now we have all the random parameter sets generated, so we need to create
%a projected trajectory for each parameter set:

n=0;
for o=1:recovery_reps
    for p=1:1000
        n=n+1;
    curve=Model_1_Forward(model_1_recovering_parameters(o,10),posttreat_time_vector(1:model_1_recovering_parameters(o,9),1),model_1_recovering_parameters(o,8),model_1_recovering_generated_parameters(n,4),model_1_recovering_generated_parameters(n,1),model_1_recovering_generated_parameters(n,2),model_1_recovering_generated_parameters(n,3),model_1_recovering_generated_parameters(n,5),model_1_recovering_generated_parameters(n,6)); 
    model_1_simulated_recovering_curves_output(1:model_1_recovering_parameters(o,9),n)=curve(:,2);
    curve=Model_2_Forward(model_2_recovering_parameters(o,10),posttreat_time_vector(1:model_2_recovering_parameters(o,9),1),model_2_recovering_parameters(o,8),model_2_recovering_generated_parameters(n,4),model_2_recovering_generated_parameters(n,1),model_2_recovering_generated_parameters(n,2),model_2_recovering_generated_parameters(n,3),model_2_recovering_generated_parameters(n,5),model_2_recovering_generated_parameters(n,6)); 
    model_2_simulated_recovering_curves_output(1:model_2_recovering_parameters(o,9),n)=curve(:,2);
    curve=Model_3_Forward(model_3_recovering_parameters(o,9),posttreat_time_vector(1:model_3_recovering_parameters(o,8),1),model_3_recovering_generated_parameters(n,4),model_3_recovering_generated_parameters(n,1),model_3_recovering_generated_parameters(n,2),model_3_recovering_generated_parameters(n,3),model_3_recovering_generated_parameters(n,5)); 
    model_3_simulated_recovering_curves_output(1:model_3_recovering_parameters(o,8),n)=curve(:,2);
    end
end

n=0;
for o=1:dying_reps
    for p=1:1000
        n=n+1;
    curve=Model_1_Death_Only_Forward(model_1_dying_parameters(o,6),posttreat_time_vector(1:model_1_dying_parameters(o,5),1),model_1_dying_parameters(o,4),model_1_dying_generated_parameters(n,1),model_1_dying_generated_parameters(n,2)); 
    model_1_simulated_dying_curves_output(1:model_1_dying_parameters(o,5),n)=curve(:,2);
    curve=Model_2_Death_Only_Forward(model_2_dying_parameters(o,6),posttreat_time_vector(1:model_2_dying_parameters(o,5),1),model_2_dying_parameters(o,4),model_2_dying_generated_parameters(n,1),model_2_dying_generated_parameters(n,2)); 
    model_2_simulated_dying_curves_output(1:model_2_dying_parameters(o,5),n)=curve(:,2);
    curve=Model_3_Death_Only_Forward(model_3_dying_parameters(o,5),posttreat_time_vector(1:model_3_dying_parameters(o,4),1),model_3_dying_generated_parameters(n,1)); 
    model_3_simulated_dying_curves_output(1:model_3_dying_parameters(o,4),n)=curve(:,2);
    end
end



%strip infinite values

for q=1:(1000*recovery_reps)
   for r=1:length(posttreat_time_vector)
       if isnan(model_1_simulated_recovering_curves_output(r,q))==0
           model_1_cleaned_simulated_recovering_curves_output(r,q)=model_1_simulated_recovering_curves_output(r,q);
       end
       if isnan(model_2_simulated_recovering_curves_output(r,q))==0
           model_2_cleaned_simulated_recovering_curves_output(r,q)=model_2_simulated_recovering_curves_output(r,q);
       end
       if isnan(model_3_simulated_recovering_curves_output(r,q))==0
           model_3_cleaned_simulated_recovering_curves_output(r,q)=model_3_simulated_recovering_curves_output(r,q);
       end
   end
end

for q=1:(1000*dying_reps)
   for r=1:length(posttreat_time_vector)
       if isnan(model_1_simulated_dying_curves_output(r,q))==0
           model_1_cleaned_simulated_dying_curves_output(r,q)=model_1_simulated_dying_curves_output(r,q);
       end
       if isnan(model_2_simulated_dying_curves_output(r,q))==0
           model_2_cleaned_simulated_dying_curves_output(r,q)=model_2_simulated_dying_curves_output(r,q);
       end
       if isnan(model_3_simulated_dying_curves_output(r,q))==0
           model_3_cleaned_simulated_dying_curves_output(r,q)=model_3_simulated_dying_curves_output(r,q);
       end
   end
end

%Find the mean and standard deviation for each ensemble of 1000 curves
%which correspond to a particular replicate, at each time point:
if recovery_reps>1
    
for n=1:recovery_reps
    for m=1:model_1_recovering_parameters(n,9)
        model_1_simulation_mean_recovering(m,n)=mean(model_1_cleaned_simulated_recovering_curves_output(m,(1000*(n-1)+1):(1000*n)));
        model_1_simulation_std_recovering(m,n)=std(model_1_cleaned_simulated_recovering_curves_output(m,(1000*(n-1)+1):(1000*n)));
        model_2_simulation_mean_recovering(m,n)=mean(model_2_cleaned_simulated_recovering_curves_output(m,(1000*(n-1)+1):(1000*n)));
        model_2_simulation_std_recovering(m,n)=std(model_2_cleaned_simulated_recovering_curves_output(m,(1000*(n-1)+1):(1000*n)));
        model_3_simulation_mean_recovering(m,n)=mean(model_3_cleaned_simulated_recovering_curves_output(m,(1000*(n-1)+1):(1000*n)));
        model_3_simulation_std_recovering(m,n)=std(model_3_cleaned_simulated_recovering_curves_output(m,(1000*(n-1)+1):(1000*n)));
    end
end

end


if dying_reps>1
    
for n=1:dying_reps
    for m=1:model_1_dying_parameters(n,5)
        model_1_simulation_mean_dying(m,n)=mean(model_1_cleaned_simulated_dying_curves_output(m,(1000*(n-1)+1):(1000*n)));
        model_1_simulation_std_dying(m,n)=std(model_1_cleaned_simulated_dying_curves_output(m,(1000*(n-1)+1):(1000*n)));
        model_2_simulation_mean_dying(m,n)=mean(model_2_cleaned_simulated_dying_curves_output(m,(1000*(n-1)+1):(1000*n)));
        model_2_simulation_std_dying(m,n)=std(model_2_cleaned_simulated_dying_curves_output(m,(1000*(n-1)+1):(1000*n)));
        model_3_simulation_mean_dying(m,n)=mean(model_3_cleaned_simulated_dying_curves_output(m,(1000*(n-1)+1):(1000*n)));
        model_3_simulation_std_dying(m,n)=std(model_3_cleaned_simulated_dying_curves_output(m,(1000*(n-1)+1):(1000*n)));
    end
end

end

%Use the mean and standard deviation values obtained above to generate 50,
%80, and 95% confidence intervals for each validation curve
if recovery_reps>1
    
for n=1:recovery_reps
    for m=1:model_1_recovering_parameters(n,9)
        model_1_recovering_50_CI_min(m,n)=model_1_simulation_mean_recovering(m,n)-0.67449*model_1_simulation_std_recovering(m,n);
        model_2_recovering_50_CI_min(m,n)=model_2_simulation_mean_recovering(m,n)-0.67449*model_2_simulation_std_recovering(m,n);
        model_3_recovering_50_CI_min(m,n)=model_3_simulation_mean_recovering(m,n)-0.67449*model_3_simulation_std_recovering(m,n);
        model_1_recovering_50_CI_max(m,n)=model_1_simulation_mean_recovering(m,n)+0.67449*model_1_simulation_std_recovering(m,n);
        model_2_recovering_50_CI_max(m,n)=model_2_simulation_mean_recovering(m,n)+0.67449*model_2_simulation_std_recovering(m,n);
        model_3_recovering_50_CI_max(m,n)=model_3_simulation_mean_recovering(m,n)+0.67449*model_3_simulation_std_recovering(m,n);
        
        model_1_recovering_80_CI_min(m,n)=model_1_simulation_mean_recovering(m,n)-1.281552*model_1_simulation_std_recovering(m,n);
        model_2_recovering_80_CI_min(m,n)=model_2_simulation_mean_recovering(m,n)-1.281552*model_2_simulation_std_recovering(m,n);
        model_3_recovering_80_CI_min(m,n)=model_3_simulation_mean_recovering(m,n)-1.281552*model_3_simulation_std_recovering(m,n);
        model_1_recovering_80_CI_max(m,n)=model_1_simulation_mean_recovering(m,n)+1.281552*model_1_simulation_std_recovering(m,n);
        model_2_recovering_80_CI_max(m,n)=model_2_simulation_mean_recovering(m,n)+1.281552*model_2_simulation_std_recovering(m,n);
        model_3_recovering_80_CI_max(m,n)=model_3_simulation_mean_recovering(m,n)+1.2815529*model_3_simulation_std_recovering(m,n);        
        
        model_1_recovering_95_CI_min(m,n)=model_1_simulation_mean_recovering(m,n)-1.95964*model_1_simulation_std_recovering(m,n);
        model_2_recovering_95_CI_min(m,n)=model_2_simulation_mean_recovering(m,n)-1.95964*model_2_simulation_std_recovering(m,n);
        model_3_recovering_95_CI_min(m,n)=model_3_simulation_mean_recovering(m,n)-1.95964*model_3_simulation_std_recovering(m,n);
        model_1_recovering_95_CI_max(m,n)=model_1_simulation_mean_recovering(m,n)+1.95964*model_1_simulation_std_recovering(m,n);
        model_2_recovering_95_CI_max(m,n)=model_2_simulation_mean_recovering(m,n)+1.95964*model_2_simulation_std_recovering(m,n);
        model_3_recovering_95_CI_max(m,n)=model_3_simulation_mean_recovering(m,n)+1.95964*model_3_simulation_std_recovering(m,n);
    end
end
end
if dying_reps>1
for n=1:dying_reps
    for m=1:model_1_dying_parameters(n,5)
        model_1_dying_50_CI_min(m,n)=model_1_simulation_mean_dying(m,n)-0.67449*model_1_simulation_std_dying(m,n);
        model_2_dying_50_CI_min(m,n)=model_2_simulation_mean_dying(m,n)-0.67449*model_2_simulation_std_dying(m,n);
        model_3_dying_50_CI_min(m,n)=model_3_simulation_mean_dying(m,n)-0.67449*model_3_simulation_std_dying(m,n);
        model_1_dying_50_CI_max(m,n)=model_1_simulation_mean_dying(m,n)+0.67449*model_1_simulation_std_dying(m,n);
        model_2_dying_50_CI_max(m,n)=model_2_simulation_mean_dying(m,n)+0.67449*model_2_simulation_std_dying(m,n);
        model_3_dying_50_CI_max(m,n)=model_3_simulation_mean_dying(m,n)+0.67449*model_3_simulation_std_dying(m,n);
        
        model_1_dying_80_CI_min(m,n)=model_1_simulation_mean_dying(m,n)-1.281552*model_1_simulation_std_dying(m,n);
        model_2_dying_80_CI_min(m,n)=model_2_simulation_mean_dying(m,n)-1.281552*model_2_simulation_std_dying(m,n);
        model_3_dying_80_CI_min(m,n)=model_3_simulation_mean_dying(m,n)-1.281552*model_3_simulation_std_dying(m,n);
        model_1_dying_80_CI_max(m,n)=model_1_simulation_mean_dying(m,n)+1.281552*model_1_simulation_std_dying(m,n);
        model_2_dying_80_CI_max(m,n)=model_2_simulation_mean_dying(m,n)+1.281552*model_2_simulation_std_dying(m,n);
        model_3_dying_80_CI_max(m,n)=model_3_simulation_mean_dying(m,n)+1.2815529*model_3_simulation_std_dying(m,n);        
        
        model_1_dying_95_CI_min(m,n)=model_1_simulation_mean_dying(m,n)-1.95964*model_1_simulation_std_dying(m,n);
        model_2_dying_95_CI_min(m,n)=model_2_simulation_mean_dying(m,n)-1.95964*model_2_simulation_std_dying(m,n);
        model_3_dying_95_CI_min(m,n)=model_3_simulation_mean_dying(m,n)-1.95964*model_3_simulation_std_dying(m,n);
        model_1_dying_95_CI_max(m,n)=model_1_simulation_mean_dying(m,n)+1.95964*model_1_simulation_std_dying(m,n);
        model_2_dying_95_CI_max(m,n)=model_2_simulation_mean_dying(m,n)+1.95964*model_2_simulation_std_dying(m,n);
        model_3_dying_95_CI_max(m,n)=model_3_simulation_mean_dying(m,n)+1.95964*model_3_simulation_std_dying(m,n);
    end
end
end

%Now compare the calculated LOO max and min for each confidence interval to
%the data value at each time point, for each replicate.
if recovery_reps>1
for n=1:recovery_reps
    for m=1:model_1_recovering_parameters(n,9)
        if posttreat_data(m,(n+1))>model_1_recovering_50_CI_max(m,n)
            model_1_recovering_50_CI_check(m,n)=0;
        elseif posttreat_data(m,(n+1))<model_1_recovering_50_CI_min(m,n)
            model_1_recovering_50_CI_check(m,n)=0;
        else
            model_1_recovering_50_CI_check(m,n)=1;
        end
        
        if posttreat_data(m,(n+1))>model_2_recovering_50_CI_max(m,n)
            model_2_recovering_50_CI_check(m,n)=0;
        elseif posttreat_data(m,(n+1))<model_2_recovering_50_CI_min(m,n)
            model_2_recovering_50_CI_check(m,n)=0;
        else
            model_2_recovering_50_CI_check(m,n)=1;
        end
        
        if posttreat_data(m,(n+1))>model_3_recovering_50_CI_max(m,n)
            model_3_recovering_50_CI_check(m,n)=0;
        elseif posttreat_data(m,(n+1))<model_3_recovering_50_CI_min(m,n)
            model_3_recovering_50_CI_check(m,n)=0;
        else
            model_3_recovering_50_CI_check(m,n)=1;
        end
        
        if posttreat_data(m,(n+1))>model_1_recovering_80_CI_max(m,n)
            model_1_recovering_80_CI_check(m,n)=0;
        elseif posttreat_data(m,(n+1))<model_1_recovering_80_CI_min(m,n)
            model_1_recovering_80_CI_check(m,n)=0;
        else
            model_1_recovering_80_CI_check(m,n)=1;
        end
        
        if posttreat_data(m,(n+1))>model_2_recovering_80_CI_max(m,n)
            model_2_recovering_80_CI_check(m,n)=0;
        elseif posttreat_data(m,(n+1))<model_2_recovering_80_CI_min(m,n)
            model_2_recovering_80_CI_check(m,n)=0;
        else
            model_2_recovering_80_CI_check(m,n)=1;
        end
        
        if posttreat_data(m,(n+1))>model_3_recovering_80_CI_max(m,n)
            model_3_recovering_80_CI_check(m,n)=0;
        elseif posttreat_data(m,(n+1))<model_3_recovering_80_CI_min(m,n)
            model_3_recovering_80_CI_check(m,n)=0;
        else
            model_3_recovering_80_CI_check(m,n)=1;
        end
        
        if posttreat_data(m,(n+1))>model_1_recovering_95_CI_max(m,n)
            model_1_recovering_95_CI_check(m,n)=0;
        elseif posttreat_data(m,(n+1))<model_1_recovering_95_CI_min(m,n)
            model_1_recovering_95_CI_check(m,n)=0;
        else
            model_1_recovering_95_CI_check(m,n)=1;
        end
        
        if posttreat_data(m,(n+1))>model_2_recovering_95_CI_max(m,n)
            model_2_recovering_95_CI_check(m,n)=0;
        elseif posttreat_data(m,(n+1))<model_2_recovering_95_CI_min(m,n)
            model_2_recovering_95_CI_check(m,n)=0;
        else
            model_2_recovering_95_CI_check(m,n)=1;
        end
        
        if posttreat_data(m,(n+1))>model_3_recovering_95_CI_max(m,n)
            model_3_recovering_95_CI_check(m,n)=0;
        elseif posttreat_data(m,(n+1))<model_3_recovering_95_CI_min(m,n)
            model_3_recovering_95_CI_check(m,n)=0;
        else
            model_3_recovering_95_CI_check(m,n)=1;
        end
    end
end
end

if dying_reps>1
for n=1:dying_reps
    for m=1:model_1_dying_parameters(n,5)
        if posttreat_data(m,(n+1))>model_1_dying_50_CI_max(m,n)
            model_1_dying_50_CI_check(m,n)=0;
        elseif posttreat_data(m,(n+1))<model_1_dying_50_CI_min(m,n)
            model_1_dying_50_CI_check(m,n)=0;
        else
            model_1_dying_50_CI_check(m,n)=1;
        end
        
        if posttreat_data(m,(n+1))>model_2_dying_50_CI_max(m,n)
            model_2_dying_50_CI_check(m,n)=0;
        elseif posttreat_data(m,(n+1))<model_2_dying_50_CI_min(m,n)
            model_2_dying_50_CI_check(m,n)=0;
        else
            model_2_dying_50_CI_check(m,n)=1;
        end
        
        if posttreat_data(m,(n+1))>model_3_dying_50_CI_max(m,n)
            model_3_dying_50_CI_check(m,n)=0;
        elseif posttreat_data(m,(n+1))<model_3_dying_50_CI_min(m,n)
            model_3_dying_50_CI_check(m,n)=0;
        else
            model_3_dying_50_CI_check(m,n)=1;
        end
        
        if posttreat_data(m,(n+1))>model_1_dying_80_CI_max(m,n)
            model_1_dying_80_CI_check(m,n)=0;
        elseif posttreat_data(m,(n+1))<model_1_dying_80_CI_min(m,n)
            model_1_dying_80_CI_check(m,n)=0;
        else
            model_1_dying_80_CI_check(m,n)=1;
        end
        
        if posttreat_data(m,(n+1))>model_2_dying_80_CI_max(m,n)
            model_2_dying_80_CI_check(m,n)=0;
        elseif posttreat_data(m,(n+1))<model_2_dying_80_CI_min(m,n)
            model_2_dying_80_CI_check(m,n)=0;
        else
            model_2_dying_80_CI_check(m,n)=1;
        end
        
        if posttreat_data(m,(n+1))>model_3_dying_80_CI_max(m,n)
            model_3_dying_80_CI_check(m,n)=0;
        elseif posttreat_data(m,(n+1))<model_3_dying_80_CI_min(m,n)
            model_3_dying_80_CI_check(m,n)=0;
        else
            model_3_dying_80_CI_check(m,n)=1;
        end
        
        if posttreat_data(m,(n+1))>model_1_dying_95_CI_max(m,n)
            model_1_dying_95_CI_check(m,n)=0;
        elseif posttreat_data(m,(n+1))<model_1_dying_95_CI_min(m,n)
            model_1_dying_95_CI_check(m,n)=0;
        else
            model_1_dying_95_CI_check(m,n)=1;
        end
        
        if posttreat_data(m,(n+1))>model_2_dying_95_CI_max(m,n)
            model_2_dying_95_CI_check(m,n)=0;
        elseif posttreat_data(m,(n+1))<model_2_dying_95_CI_min(m,n)
            model_2_dying_95_CI_check(m,n)=0;
        else
            model_2_dying_95_CI_check(m,n)=1;
        end
        
        if posttreat_data(m,(n+1))>model_3_dying_95_CI_max(m,n)
            model_3_dying_95_CI_check(m,n)=0;
        elseif posttreat_data(m,(n+1))<model_3_dying_95_CI_min(m,n)
            model_3_dying_95_CI_check(m,n)=0;
        else
            model_3_dying_95_CI_check(m,n)=1;
        end
    end
end
end
%Output all of the relevant data and calculations for downstream use:
if recovery_reps>=1
writematrix(model_1_recovering_generated_parameters,output_filename,'Sheet','Mod 1 Recovery Sim Params','Range','B2');
writematrix(model_2_recovering_generated_parameters,output_filename,'Sheet','Mod 2 Recovery Sim Params','Range','B2');
writematrix(model_3_recovering_generated_parameters,output_filename,'Sheet','Mod 3 Recovery Sim Params','Range','B2');
writematrix(model_1_cleaned_simulated_recovering_curves_output,output_filename,'Sheet','Mod 1 Recovery Sim Curves','Range','B2');
writematrix(model_2_cleaned_simulated_recovering_curves_output,output_filename,'Sheet','Mod 2 Recovery Sim Curves','Range','B2');
writematrix(model_3_cleaned_simulated_recovering_curves_output,output_filename,'Sheet','Mod 2 Recovery Sim Curves','Range','B2');
end

if recovery_reps>1
writematrix([model_1_recovering_50_CI_check,model_1_recovering_80_CI_check,model_1_recovering_95_CI_check],output_filename,'Sheet','Mod 1 Rec CI Check','Range','B2');
writematrix([model_2_recovering_50_CI_check,model_2_recovering_80_CI_check,model_2_recovering_95_CI_check],output_filename,'Sheet','Mod 2 Rec CI Check','Range','B2');
writematrix([model_3_recovering_50_CI_check,model_3_recovering_80_CI_check,model_3_recovering_95_CI_check],output_filename,'Sheet','Mod 3 Rec CI Check','Range','B2');
end

if dying_reps>=1
writematrix(model_1_dying_generated_parameters,output_filename,'Sheet','Mod 1 Dying Sim Params','Range','B2');
writematrix(model_2_dying_generated_parameters,output_filename,'Sheet','Mod 2 Dying Sim Params','Range','B2');
writematrix(model_3_dying_generated_parameters,output_filename,'Sheet','Mod 3 Dying Sim Params','Range','B2');
writematrix(model_1_cleaned_simulated_dying_curves_output,output_filename,'Sheet','Mod 1 Dying Sim Curves','Range','B2');
writematrix(model_2_cleaned_simulated_dying_curves_output,output_filename,'Sheet','Mod 2 Dying Sim Curves','Range','B2');
writematrix(model_3_cleaned_simulated_dying_curves_output,output_filename,'Sheet','Mod 3 Dying Sim Curves','Range','B2');
end
if dying_reps>1
writematrix([model_1_dying_50_CI_check,model_1_dying_80_CI_check,model_1_dying_95_CI_check],output_filename,'Sheet','Mod 1 Die CI Check','Range','B2');
writematrix([model_2_dying_50_CI_check,model_2_dying_80_CI_check,model_2_dying_95_CI_check],output_filename,'Sheet','Mod 2 Die CI Check','Range','B2');
writematrix([model_3_dying_50_CI_check,model_3_dying_80_CI_check,model_3_dying_95_CI_check],output_filename,'Sheet','Mod 3 Die CI Check','Range','B2');
end
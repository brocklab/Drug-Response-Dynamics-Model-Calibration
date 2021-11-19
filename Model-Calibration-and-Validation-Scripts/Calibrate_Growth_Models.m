clear all
close all
clc

%Read in the data file with the elapsed time vs. cell number data.  This
%should be an excel sheet with the time vector in the first column and the
%n cell number vectors in columns 2-(n+1).  There should be no headers

data = xlsread('GH1909_200nM_working_mod_11_09.xlsx');
num_time_points = length(data);

%Read in all of the auxiliary information:  t_r values, etc.

conditions = xlsread('GH1909_200nM_conditions.xlsx');

%Specify the name you want to use for the output file

output_filename = 'GH1909_200nM_model_calibration_output_v2.xlsx';

%Specify the name you want to use for a separate output file containing the
%cell number curves for the resistant and sensitive cell compartments.

RS_output_filename = 'GH1909_200nM_model_calibration_RS_output_v2.xlsx';

%Extract key data from conditions sheet

num_rep = conditions(1,1);
time_of_treat = conditions(10,1);
recovery_boolean = conditions(:,2);
t_r_vec = conditions(:,3);
model_1_z_init_guess = conditions(1:5,6);
model_2_z_init_guess = conditions(1:5,9);
model_3_z_init_guess = conditions(1:4,12);
model_1_lower_bounds = conditions(1:5,5);
model_2_lower_bounds = conditions(1:5,8);
model_3_lower_bounds = conditions(1:4,11);
model_1_upper_bounds = conditions(1:5,7);
model_2_upper_bounds = conditions(1:5,10);
model_3_upper_bounds = conditions(1:4,13);
model_1_Death_Only_z_init_guess = conditions(1:2,16);
model_2_Death_Only_z_init_guess = conditions(1:2,19);
model_3_Death_Only_z_init_guess = conditions(1,22);
model_1_Death_Only_lower_bounds = conditions(1:2,15);
model_2_Death_Only_lower_bounds = conditions(1:2,18);
model_3_Death_Only_lower_bounds = conditions(1,21);
model_1_Death_Only_upper_bounds = conditions(1:2,17);
model_2_Death_Only_upper_bounds = conditions(1:2,20);
model_3_Death_Only_upper_bounds = conditions(1,23);

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

%Create a matrix for the residual sum of squares for each replicate, for each
%model for use in calculating the AIC

RSS=zeros(num_rep,6);

%Calibrate each cell number curve to the three models, saving the calibrated
%parameters and RSS to output at the end.

for m = 1:num_rep
    g_0 = g_0_vector(m);
    t_r = t_r_vec(m);
    
    %Set the lower bound of N_max such that it is at least 2.1 times the
    %initial cell number; this facilitates finding Critical Time later.
    N_max_lower_bound_for_t_crit = 2.1*N_0_vector(1,m);
    if model_1_lower_bounds(4,1)<N_max_lower_bound_for_t_crit
        model_1_lower_bounds(4,1)=N_max_lower_bound_for_t_crit;
    end
    if model_2_lower_bounds(4,1)<(2.1*N_0_vector(1,m))
        model_2_lower_bounds(4,1)=(2.1*N_0_vector(1,m));
    end
    if model_3_lower_bounds(4,1)<(2.1*N_0_vector(1,m))
        model_3_lower_bounds(4,1)=(2.1*N_0_vector(1,m));
    end
    
    time_vector = posttreat_data(1:trunc_points(m,1),1);
    count_vector = posttreat_data(1:trunc_points(m,1),m+1);
    model_1_z_output(m,:) = lsqnonlin(@(z) Model_1_Lst_Sq_Function(z,count_vector,time_vector,N_0_vector(1,m),g_0,t_r),model_1_z_init_guess,model_1_lower_bounds,model_1_upper_bounds);
    model_2_z_output(m,:) = lsqnonlin(@(z) Model_2_Lst_Sq_Function(z,count_vector,time_vector,N_0_vector(1,m),g_0,t_r),model_2_z_init_guess,model_2_lower_bounds,model_2_upper_bounds);
    model_3_z_output(m,:) = lsqnonlin(@(z) Model_3_Lst_Sq_Function(z,count_vector,time_vector,N_0_vector(1,m),t_r),model_3_z_init_guess,model_3_lower_bounds,model_3_upper_bounds);
    model_1_Death_Only_z_output(m,:) = lsqnonlin(@(z) Model_1_Death_Only_Lst_Sq_Function(z,count_vector,time_vector,N_0_vector(1,m),g_0),model_1_Death_Only_z_init_guess,model_1_Death_Only_lower_bounds,model_1_Death_Only_upper_bounds);
    model_2_Death_Only_z_output(m,:) = lsqnonlin(@(z) Model_2_Death_Only_Lst_Sq_Function(z,count_vector,time_vector,N_0_vector(1,m),g_0),model_2_Death_Only_z_init_guess,model_2_Death_Only_lower_bounds,model_2_Death_Only_upper_bounds);
    model_3_Death_Only_z_output(m,:) = lsqnonlin(@(z) Model_3_Death_Only_Lst_Sq_Function(z,count_vector,time_vector,N_0_vector(1,m)),model_3_Death_Only_z_init_guess,model_3_Death_Only_lower_bounds,model_3_Death_Only_upper_bounds);
    extended_t_vector = transpose([0:8760]);
    extended_t_vector = extended_t_vector+posttreat_data(1,1);
    model_1_intermediate = Model_1_Forward(N_0_vector(m),extended_t_vector,g_0,model_1_z_output(m,4),model_1_z_output(m,1),model_1_z_output(m,2),model_1_z_output(m,3),model_1_z_output(m,5),t_r);
    model_2_intermediate = Model_2_Forward(N_0_vector(m),extended_t_vector,g_0,model_2_z_output(m,4),model_2_z_output(m,1),model_2_z_output(m,2),model_2_z_output(m,3),model_2_z_output(m,5),t_r);
    model_3_intermediate = Model_3_Forward(N_0_vector(m),extended_t_vector,model_3_z_output(m,4),model_3_z_output(m,1),model_3_z_output(m,2),model_3_z_output(m,3),t_r);
    model_1_Death_Only_intermediate = Model_1_Death_Only_Forward(N_0_vector(m),extended_t_vector,g_0,model_1_Death_Only_z_output(m,1),model_1_Death_Only_z_output(m,2));
    model_2_Death_Only_intermediate = Model_2_Death_Only_Forward(N_0_vector(m),extended_t_vector,g_0,model_2_Death_Only_z_output(m,1),model_2_Death_Only_z_output(m,2));
    model_3_Death_Only_intermediate = Model_3_Death_Only_Forward(N_0_vector(m),extended_t_vector,model_3_Death_Only_z_output(m,1));
    model_1_final_forward(:,1) = extended_t_vector;
    model_1_final_forward(:,m+1) = model_1_intermediate(:,2);
    model_2_final_forward(:,1) = extended_t_vector;
    model_2_final_forward(:,m+1) = model_2_intermediate(:,2);
    model_3_final_forward(:,1) = extended_t_vector;
    model_3_final_forward(:,m+1) = model_3_intermediate(:,2);

    model_1_RSintermediate = Model_1_RSForward(N_0_vector(m),extended_t_vector,g_0,model_1_z_output(m,4),model_1_z_output(m,1),model_1_z_output(m,2),model_1_z_output(m,3),model_1_z_output(m,5),t_r);
    model_2_RSintermediate = Model_2_RSForward(N_0_vector(m),extended_t_vector,g_0,model_2_z_output(m,4),model_2_z_output(m,1),model_2_z_output(m,2),model_2_z_output(m,3),model_2_z_output(m,5),t_r);
    model_3_RSintermediate = Model_3_RSForward(N_0_vector(m),extended_t_vector,model_3_z_output(m,4),model_3_z_output(m,1),model_3_z_output(m,2),model_3_z_output(m,3),t_r);


    model_1_R_forward(:,1) = extended_t_vector;
    model_1_R_forward(:,m+1) = model_1_RSintermediate(:,3);
    model_2_R_forward(:,1) = extended_t_vector;
    model_2_R_forward(:,m+1) = model_2_RSintermediate(:,3);
    model_3_R_forward(:,1) = extended_t_vector;
    model_3_R_forward(:,m+1) = model_3_RSintermediate(:,3);

    model_1_S_forward(:,1) = extended_t_vector;
    model_1_S_forward(:,m+1) = model_1_RSintermediate(:,4);
    model_2_S_forward(:,1) = extended_t_vector;
    model_2_S_forward(:,m+1) = model_2_RSintermediate(:,4);
    model_3_S_forward(:,1) = extended_t_vector;
    model_3_S_forward(:,m+1) = model_3_RSintermediate(:,4);


    model_1_Death_Only_final_forward(:,1) = extended_t_vector;
    model_1_Death_Only_final_forward(:,m+1) = model_1_Death_Only_intermediate(:,2);
    model_2_Death_Only_final_forward(:,1) = extended_t_vector;
    model_2_Death_Only_final_forward(:,m+1) = model_2_Death_Only_intermediate(:,2);
    model_3_Death_Only_final_forward(:,1) = extended_t_vector;
    model_3_Death_Only_final_forward(:,m+1) = model_3_Death_Only_intermediate(:,2);
    model_1_RSS_intermediate = Model_1_Lst_Sq_Function(model_1_z_output(m,1:5),count_vector,time_vector,N_0_vector(m),g_0,t_r);
    model_1_RSS_int_2 = model_1_RSS_intermediate.^2;
    model_2_RSS_intermediate = Model_2_Lst_Sq_Function(model_2_z_output(m,1:5),count_vector,time_vector,N_0_vector(m),g_0,t_r);
    model_2_RSS_int_2 = model_2_RSS_intermediate.^2;
    model_3_RSS_intermediate = Model_3_Lst_Sq_Function(model_3_z_output(m,1:4),count_vector,time_vector,N_0_vector(m),t_r);
    model_3_RSS_int_2 = model_3_RSS_intermediate.^2;
    model_1_Death_Only_RSS_intermediate = Model_1_Death_Only_Lst_Sq_Function(model_1_Death_Only_z_output(m,1:2),count_vector,time_vector,N_0_vector(m),g_0);
    model_1_Death_Only_RSS_int_2 = model_1_Death_Only_RSS_intermediate.^2;
    model_2_Death_Only_RSS_intermediate = Model_2_Death_Only_Lst_Sq_Function(model_2_Death_Only_z_output(m,1:2),count_vector,time_vector,N_0_vector(m),g_0);
    model_2_Death_Only_RSS_int_2 = model_2_Death_Only_RSS_intermediate.^2;
    model_3_Death_Only_RSS_intermediate = Model_3_Death_Only_Lst_Sq_Function(model_3_Death_Only_z_output(m,1),count_vector,time_vector,N_0_vector(m));
    model_3_Death_Only_RSS_int_2 = model_3_Death_Only_RSS_intermediate.^2;
    
    for n = 1:trunc_points(m,1)
        RSS(m,1) = RSS(m,1)+model_1_RSS_int_2(n);
        RSS(m,2) = RSS(m,2)+model_2_RSS_int_2(n);
        RSS(m,3) = RSS(m,3)+model_3_RSS_int_2(n);
        RSS(m,4) = RSS(m,4)+model_1_Death_Only_RSS_int_2(n);
        RSS(m,5) = RSS(m,5)+model_2_Death_Only_RSS_int_2(n);
        RSS(m,6) = RSS(m,6)+model_3_Death_Only_RSS_int_2(n);
    end
    
    %Reset the lower bound of N_max for the next iteration
    model_1_lower_bounds(4,1)=conditions(4,5);
    model_2_lower_bounds(4,1)=conditions(4,8);
    model_3_lower_bounds(4,1)=conditions(4,11);
    
end

%Use the RSS values to calculate the AIC for model 1, for each
%replicate:

AIC = zeros(num_rep,6);
for o = 1:num_rep
    AIC(o,1) = 2*5+trunc_points(o)*log(RSS(o,1)/trunc_points(o));
    AIC(o,2) = 2*5+trunc_points(o)*log(RSS(o,2)/trunc_points(o));
    AIC(o,3) = 2*4+trunc_points(o)*log(RSS(o,3)/trunc_points(o));
    AIC(o,4) = 2*2+trunc_points(o)*log(RSS(o,4)/trunc_points(o));
    AIC(o,5) = 2*2+trunc_points(o)*log(RSS(o,5)/trunc_points(o));
    AIC(o,6) = 2*1+trunc_points(o)*log(RSS(o,6)/trunc_points(o));
end

%Use the AIC values to check which model performs best for each replicate;
%we will use the recovery boolean input to check the recovery models in
%replicates where the recovery occurred and the death only model in
%replicates where recovery did not occur.  Note that it is possible for
%multiple models to have the same AIC value (within rounding errors) in
%some cases, so the total number of "best" cases identified may exceed the
%number of replicates.

best_model = zeros(6,1);

minimum_AIC_recovery = zeros(num_rep,1);

for n=1:num_rep
    minimum_AIC_recovery(n,1)=AIC(n,1);
    if AIC(n,2)<minimum_AIC_recovery(n,1)
        minimum_AIC_recovery(n,1)=AIC(n,2);
    end
    if AIC(n,3)<minimum_AIC_recovery(n,1)
        minimum_AIC_recovery(n,1)=AIC(n,3);
    end
end

minimum_AIC_dying = zeros(num_rep,1);

for n=1:num_rep
    minimum_AIC_dying(n,1)=AIC(n,4);
    if AIC(n,5)< minimum_AIC_dying(n,1)
         minimum_AIC_dying(n,1)=AIC(n,5);
    end
    if AIC(n,6)< minimum_AIC_dying(n,1)
         minimum_AIC_dying(n,1)=AIC(n,6);
    end
end

recovery_reps = 0;
dying_reps = 0;

for n=1:num_rep
    if recovery_boolean(n,1) == 1
        recovery_reps = recovery_reps+1;
        if AIC(n,1)<=minimum_AIC_recovery(n,1)
            best_model(1,1)=best_model(1,1)+1;
        end
        if AIC(n,2)<=minimum_AIC_recovery(n,1)
            best_model(2,1)=best_model(2,1)+1;
        end
        if AIC(n,3)<=minimum_AIC_recovery(n,1)
            best_model(3,1)=best_model(3,1)+1;
        end
    else
        dying_reps = dying_reps+1;
        if AIC(n,4)<=minimum_AIC_dying(n,1)
            best_model(4,1)=best_model(4,1)+1;
        end
        if AIC(n,5)<=minimum_AIC_dying(n,1)
            best_model(5,1)=best_model(5,1)+1;
        end
        if AIC(n,6)<=minimum_AIC_dying(n,1)
            best_model(6,1)=best_model(6,1)+1;
        end
    end
end


%Find the critical time, where N = 2*N_0, for each replicate, for models 1,
%2, and 3 (but not the death only versions because they never recover)

t_crit = zeros(num_rep,3);
for p = 1:num_rep
    for q = 1:8760
        if model_1_final_forward(q,p+1)>=2*N_0_vector(1,p)
            t_crit(p,1) = model_1_final_forward(q,1);
            break
        end
    end
end

for p = 1:num_rep
    for q = 1:8760
        if model_2_final_forward(q,p+1)>=2*N_0_vector(1,p)
            t_crit(p,2) = model_2_final_forward(q,1);
            break
        end
    end
end

for p = 1:num_rep
    for q = 1:8760
        if model_3_final_forward(q,p+1)>=2*N_0_vector(1,p)
            t_crit(p,3) = model_3_final_forward(q,1);
            break
        end
    end
end

%Write all the calibrated values into an output file:

writematrix( model_1_R_forward,RS_output_filename,'Sheet','Model 1 R Curves','Range','B2');
writematrix( model_2_R_forward,RS_output_filename,'Sheet','Model 2 R Curves','Range','B2');
writematrix( model_3_R_forward,RS_output_filename,'Sheet','Model 3 R Curves','Range','B2');
writematrix( model_1_S_forward,RS_output_filename,'Sheet','Model 1 S Curves','Range','B2');
writematrix( model_2_S_forward,RS_output_filename,'Sheet','Model 2 S Curves','Range','B2');
writematrix( model_3_S_forward,RS_output_filename,'Sheet','Model 3 S Curves','Range','B2');


writematrix(model_1_final_forward,output_filename,'Sheet','Model 1 Curves','Range','B2');
writematrix(model_1_z_output,output_filename,'Sheet','Model 1 Parameters','Range','B2');
writematrix(model_2_final_forward,output_filename,'Sheet','Model 2 Curves','Range','B2');
writematrix(model_2_z_output,output_filename,'Sheet','Model 2 Parameters','Range','B2');
writematrix(model_3_final_forward,output_filename,'Sheet','Model 3 Curves','Range','B2');
writematrix(model_3_z_output,output_filename,'Sheet','Model 3 Parameters','Range','B2');

writematrix(model_1_Death_Only_final_forward,output_filename,'Sheet','Model 1 Death Only Curves','Range','B2');
writematrix(model_1_Death_Only_z_output,output_filename,'Sheet','Model 1 Death Only Parameters','Range','B2');
writematrix(model_2_Death_Only_final_forward,output_filename,'Sheet','Model 2 Death Only Curves','Range','B2');
writematrix(model_2_Death_Only_z_output,output_filename,'Sheet','Model 2 Death Only Parameters','Range','B2');
writematrix(model_3_Death_Only_final_forward,output_filename,'Sheet','Model 3 Death Only Curves','Range','B2');
writematrix(model_3_Death_Only_z_output,output_filename,'Sheet','Model 3 Death Only Parameters','Range','B2');

writematrix(AIC,output_filename,'Sheet','AIC Values','Range','B2');
writematrix(best_model,output_filename,'Sheet','AIC Values','Range','J2');
writematrix([recovery_reps;dying_reps],output_filename,'Sheet','AIC Values','Range','J9');

writematrix(t_crit,output_filename,'Sheet','Critical Times','Range','B2');
writematrix(g_0_vector,output_filename,'Sheet','g_0 Values','Range','B2');

writematrix(posttreat_data(:,1),output_filename,'Sheet','Post Treat Time Vec','Range','A1');

%Add labels to the output file:

writecell({'f_r','g_r','k_d','N_max','t_d'},output_filename,'Sheet','Model 1 Parameters','Range','B1')
writecell({'f_r', 'g_r', 'k_d', 'N_max', 't_d'},output_filename,'Sheet','Model 2 Parameters','Range','B1:F1')
writecell({'f_r', 'g_r', 'k_d', 'N_max'},output_filename,'Sheet','Model 3 Parameters','Range','B1:E1')

writecell({'k_d','t_d'},output_filename,'Sheet','Model 1 Death Only Parameters','Range','B1')
writecell({'k_d','t_d'},output_filename,'Sheet','Model 2 Death Only Parameters','Range','B1')
writecell({'k_d'},output_filename,'Sheet','Model 3 Death Only Parameters','Range','B1')

writecell({'Model 1 AIC', 'Model 2 AIC', 'Model 3 AIC', 'Model 1 Death Only AIC', 'Model 2 Death Only AIC', 'Model 3 Death Only AIC'},output_filename,'Sheet','AIC Values','Range','B1:G1')
writecell({'Model 1 Best, Recovering'; 'Model 2 Best, Recovering'; 'Model 3 Best, Recovering'; 'Model 1 Best, Dying'; 'Model 2 Best, Dying'; 'Model 3 Best, Dying'; '-'; 'Replicates which recover'; 'Replicates which die'},output_filename,'Sheet','AIC Values','Range','I2:I10')


writecell({'Model 1 Critical Time', 'Model 2 Critical Time', 'Model 3 Critical Time'},output_filename,'Sheet','Critical Times','Range','B1:D1')


%Output a figure as confirmation that the script is finished running
f1 = figure;
figure(f1);
plot(posttreat_data(:,1),posttreat_data(:,2),extended_t_vector,model_1_final_forward(:,2));
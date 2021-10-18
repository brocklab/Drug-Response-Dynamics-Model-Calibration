%This function accepts the values for parameters and experimental
%conditions and generates a predicted cell-number-over-time curve for that
%set of values.  This version uses a linear delay on sensitive cell
%death (model 2), and assumes that 100% of the population is sensitive.
%This is used when no cells in the replicate recovered during the course of
%the experiment, as determined by a manual check of the final microscope
%images.

function [forward_model] = Model_2_Death_Only_Forward(N_init,t_vector,g_0,k_D,t_d)

%Make a matrix to hold cell number values
lgth = length(t_vector);
cell_number = zeros(lgth,1);

%identify initial time, so that the linear decay in the death term
%starts at the time of treatment.
t_init = t_vector(1,1);

%Fill in the first row with the initial conditions, based on N_init and
%f_r.  First column is total, second is R cells and third is S cells
cell_number(1,1) = N_init;

%Use a for loop to fill in the rest of the time points based on the model
for i = 2:lgth
    if (t_vector(i)-t_init)<t_d
        k = k_D-(k_D+g_0)*(1-(t_vector(i)-t_init)/t_d);
    else
        k = k_D;
    end
    cell_number(i,1) = cell_number(i-1,1)-cell_number(i-1,1)*k*(t_vector(i)-t_vector(i-1));
end

%Output the time and total cell number vectors
output = zeros(lgth,2);
output(:,1) = t_vector;
output(:,2) = cell_number(:,1);
%In this output, column 1 is the time vector, column 2 is the total cell
%number

[forward_model] = output;
end
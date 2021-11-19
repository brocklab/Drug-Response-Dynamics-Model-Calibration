%This function accepts the values for parameters and experimental
%conditions and generates a predicted cell-number-over-time curve for that
%set of values.  This version includes a proliferation delay term t_r and
%an exponential time delay on the cell death term (model 1).

function [forward_model] = Model_1_RSForward(N_init,t_vector,g_0,N_max,f_r,g_r,k_D,t_d,t_r)

%Make a matrix to hold cell number values
lgth = length(t_vector);
cell_number = zeros(lgth,3);

%identify initial time, so that the exponential decay in the death term
%starts at the time of treatment.
t_init = t_vector(1,1);

%Fill in the first row with the initial conditions, based on N_init and
%f_r.  First column is total, second is R cells and third is S cells
cell_number(1,1) = N_init;
cell_number(1,2) = f_r*N_init;
cell_number(1,3) = (1-f_r)*N_init;

%Use a for loop to fill in the rest of the time points based on the model
for i = 2:lgth
    k = k_D-(k_D+g_0)*exp(-t_d*(t_vector(i)-t_init));
        if t_vector(i)<(t_r+t_init)
            cell_number(i,2) = cell_number(i-1,2);
            cell_number(i,3) = cell_number(i-1,3)-cell_number(i-1,3)*k*(t_vector(i)-t_vector(i-1));
            cell_number(i,1) = cell_number(i,2)+cell_number(i,3);
        else
            cell_number(i,2) = cell_number(i-1,2)+cell_number(i-1,2)*g_r*(1-cell_number(i-1,1)/N_max)*(t_vector(i)-t_vector(i-1));
            cell_number(i,3) = cell_number(i-1,3)-cell_number(i-1,3)*k*(t_vector(i)-t_vector(i-1));
            cell_number(i,1) = cell_number(i,2)+cell_number(i,3);
        end
end

%Output the time and total cell number vectors
output = zeros(lgth,4);
output(:,1) = t_vector;
output(:,2) = cell_number(:,1);
output(:,3) = cell_number(:,2);
output(:,4) = cell_number(:,3);
%In this output, column 1 is the time vector, column 2 is the total cell
%number

[forward_model] = output;
end
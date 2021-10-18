function [forward_model] = logistic_growth(N_init,t_vector,g_0,N_max)

%Make a matrix to hold cell number values
lgth = length(t_vector);
cell_number = zeros(lgth,1);

%Fill in the first row with the initial conditions, based on N_init and
%f_r.  First column is total, second is R cells and third is S cells
cell_number(1,1) = N_init;

%Use a for loop to fill in the rest of the time points based on the model
for i = 2:lgth
    cell_number(i,1) = cell_number(i-1,1)+cell_number(i-1,1)*g_0*(1-cell_number(i-1,1)/N_max)*(t_vector(i)-t_vector(i-1));
end

%Output the time and total cell number vectors
output = zeros(lgth,2);
output(:,1) = t_vector;
output(:,2) = cell_number(:,1);
%In this output, column 1 is the time vector, column 2 is the total cell
%number

[forward_model] = output;
end
%This function accepts parameter guesses for the phenomenological model, as
%well as count and time vectors, uses the parameter guesses to generate a
%forward model, and subtracts the data from the model to give an error
%vector that can be used with lsqnonlin to optimize the parameter guess.

function [error_vector] = logistic_growth_lst_sq_fn(z,count,t_vector)
g_0 = z(1);
N_max = z(2);
N_init = count(1,1);

prediction = logistic_growth(N_init,t_vector,g_0,N_max);
prediction_vector = prediction(:,2);

[error_vector] = prediction_vector-count;
end
%This function accepts parameter guesses for the phenomenological model, as
%well as count and time vectors, uses the parameter guesses to generate a
%forward model, and subtracts the data from the model to give an error
%vector that can be used with lsqnonlin to optimize the parameter guess.
%This version is for model 3, which incorporates a proliferation delay
%(t_r) and a constant sensitive cell death rate.

function [error_vector] = Model_3_Lst_Sq_Function(z,count,t_vector,N_init,t_r)
f_r = z(1);
g_r = z(2);
k_D = z(3);
N_max = z(4);

prediction = Model_3_Forward(N_init,t_vector,N_max,f_r,g_r,k_D,t_r);
prediction_vector = prediction(:,2);

[error_vector] = prediction_vector-count;
end
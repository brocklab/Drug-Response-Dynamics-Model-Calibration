%This function accepts parameter guesses for the phenomenological model, as
%well as count and time vectors, uses the parameter guesses to generate a
%forward model, and subtracts the data from the model to give an error
%vector that can be used with lsqnonlin to optimize the parameter guess.
%This version is for model 2, which incorporates a proliferation delay
%(t_r) and a linear delay on sensitive cell death.

function [error_vector] = Model_2_Lst_Sq_Function(z,count,t_vector,N_init,g_0,t_r)
f_r = z(1);
g_r = z(2);
k_D = z(3);
N_max = z(4);
t_d = z(5);


prediction = Model_2_Forward(N_init,t_vector,g_0,N_max,f_r,g_r,k_D,t_d,t_r);
prediction_vector = prediction(:,2);

[error_vector] = prediction_vector-count;
end
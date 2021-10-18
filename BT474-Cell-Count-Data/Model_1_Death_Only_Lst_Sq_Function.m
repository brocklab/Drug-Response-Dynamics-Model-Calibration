%This function accepts parameter guesses for the phenomenological model, as
%well as count and time vectors, uses the parameter guesses to generate a
%forward model, and subtracts the data from the model to give an error
%vector that can be used with lsqnonlin to optimize the parameter guess.
%This version is for model 1, with 100% of the cells assumed to be
%sensitive, and is used only when manual observation confirms that no cells
%recovered during the experiment.

function [error_vector] = Model_1_Death_Only_Lst_Sq_Function(z,count,t_vector,N_init,g_0)
k_D = z(1);
t_d = z(2);

prediction = Model_1_Death_Only_Forward(N_init,t_vector,g_0,k_D,t_d);
prediction_vector = prediction(:,2);

[error_vector] = prediction_vector-count;
end
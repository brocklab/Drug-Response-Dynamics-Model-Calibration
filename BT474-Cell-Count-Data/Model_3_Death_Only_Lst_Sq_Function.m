%This function accepts parameter guesses for the phenomenological model, as
%well as count and time vectors, uses the parameter guesses to generate a
%forward model, and subtracts the data from the model to give an error
%vector that can be used with lsqnonlin to optimize the parameter guess.
%This version is for model 3, with 100% of the cells assumed to be
%sensitive, and is used only when manual observation confirms that no cells
%recovered during the experiment.

function [error_vector] = Model_3_Death_Only_Lst_Sq_Function(z,count,t_vector,N_init)
k_D = z(1);

prediction = Model_3_Death_Only_Forward(N_init,t_vector,k_D);
prediction_vector = prediction(:,2);

[error_vector] = prediction_vector-count;
end
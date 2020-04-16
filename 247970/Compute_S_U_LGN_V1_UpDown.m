function [ S1, U1, S_L, U_L, S1_history] = Compute_S_U_LGN_V1_UpDown( S1, U1, S_L, U_L, ...
                    X, A_Up, A_Down, lambda1, s_b, uEta, nU, threshType, s1Max, sL_Max, histFlag)
% [ S1, U1, S_L, U_L] = Compute_S_U_LGN_V1_UpDown( S1, U1, S_L, U_L, X, A_Up, A_Down, lambda1, s_b, uEta, nU, threshType, s1Max, sL_Max)
% This function computes the membrane potentials and firing rates of simple
% cells and LGN cells in the model
%
% S1, U1: firing rates and membrane potentials of V1 simple cells
% S_L, U_L: firing rates and membrane potentials of LGN cells
% X: ON and OFF input for LGN cells from early visual system
% A_Up: feedforward connections between LGN and V1 simple cells
% A_Down: feedback connections between LGN and V1 simple cells
% lambda1: threshold of simple cells that controls the sparseness of simple cells
% s_b: background firing rate
% uEta: updating rate of membrane potentials U
% nU: number of iterations of calculating membrane potentials U
% threshType: type of thresholding function that computes firing rates of simple cells from membrane potentials


for n = 1 : nU
    
    V_Leak = - A_Up' * repmat(s_b, size(S_L,1), size(S_L,2)); % membrane potential caused by leakage current
    
    % Dynamics of V1 simple cells: compute membrane potentials
    U1 = (1-uEta) * U1 + uEta * ( V_Leak + A_Up' * S_L + 1*S1 );
    
    % record the trajectory of the first simple cell
    if exist('histFlag','var')
        if histFlag==1
            %         S_history(:,i,:)=S;
            %         U_history(:,i,:)=U;
            S1_history(n,:) = S1(:,1);
        end
    end
    
    % Compute firing rates of simple cells S1 by thresholding U1
    if isequal( threshType, 'soft' ) % soft thresholding at lambda
        S1 = wthresh( U1, 's', lambda1 );
    elseif isequal( threshType, 'hard' ) % hard thresholding at lambda
        S1 = wthresh( U1, 'h', lambda1);
    elseif isequal( threshType, 'sigmoid' ) % sigmoid thresholding at lambda
        alpha = 1/20; beta = 20;
        S1 = alpha * log( 1 + exp( beta*(U1-lambda1) ) );
    else
        S1 = max(U1-lambda1,0); % non-negative soft thresholding at lambda
    end
    S1 = min(S1, s1Max); % keep the response below maximal firing rate
    
    % Dynamics of LGN cells: compute membrane potentials
    U_L = (1-uEta) * U_L + uEta * (X + A_Down * S1 + s_b);
    
    % Compute firing rates of LGN cells S_L by rectifying U_L
    S_L = max(U_L,0);
    
    S_L = min(S_L, sL_Max); % keep the response below maximal firing rate
    
end
    
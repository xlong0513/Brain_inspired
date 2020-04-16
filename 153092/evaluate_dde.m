function [t, r_r, t_zero] = evaluate_dde(lgn_struct, in_struct)

% EVALUATE_DDE - Function to evaluate a linear/nonlinear DDE with a
% given drive function. This function is only needed for feedback models. 
%
% Input parameters are: 
%   in_struct - struct with fields:
%          form: 'vec', 'fnc' or 'ilfnc'
%          t_in: column vector
%          r_in: column vector
%           t_g: column vector
%           r_g: column vector
%
%   lgn_struct - struct with fields and default values 
%       eta_ffi: 0.5000
%        tau_rg: 10
%       tau_rig: 50
%     Delta_rig: 0
%        w_fbON: -1.2000
%      w_fbOFFx: 1.2000
%        tau_rc: 10
%      Delta_fb: 0
%         l_cON: 0
%        l_cOFF: 0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Extracting parameters
  w_fbON = lgn_struct.w_fbON;
  w_fbOFFx = lgn_struct.w_fbOFFx;
  tau_rc = lgn_struct.tau_rc;
  Delta_fb = lgn_struct.Delta_fb;
  l_cON = lgn_struct.l_cON;
  l_cOFF = lgn_struct.l_cOFF;
  % Extracting input/misc parameters
  if isfield(in_struct,'mode')
      mode = in_struct.mode;
  else
      mode='N/A';
  end
  form = in_struct.form;
  if strcmp(form,'vec')
      t_g = in_struct.t_g;
      r_g = in_struct.r_g;
  elseif strcmp(form,'fnc')
      h = in_struct.h;
%      B = in_struct.B;
  end  
  
  if isfield(in_struct,'tstop')
      tstop=in_struct.tstop;
  else
      tstop = t_g(end);
  end
    
  tspan = [0,tstop];

  % Setting some options for the DDE-solver. The 'MaxStep' option is very
  % demanding. The options 'RelTol' and 'AbsTol' are less demanding.
  
  opt = ddeset('RelTol',1e-4,'MaxStep', 1,'Events',@events);
    
  % opt = ddeset('InitialStep',0.05,);

  % Calling the DDE-solver with the nested function 'dde_eval' and history 
  % 'dde_history'. 

  sol = dde23(@dde_eval, Delta_fb, @dde_history, tspan, opt);

  % Extracting values from the 'sol' struct.

  t = (sol.x)'; % Time
  z = (sol.y)'; % Feedback contributions to the LGN firing rate

  % The feedforward contribution to the LGN firing rate is an analytical
  % expression contained in 'ff_contrib'.
  
  if strcmp(form,'vec')
      r_bar = interp1(t_g,r_g,t,'linear');
  elseif strcmp(form,'fnc')
      r_bar = feval(h, lgn_struct, in_struct, t);
  end
          
  %r_bar = feval(@ff_contrib, ff_params, B, t_out);

  % r_r is the LGN relay cell response, which is the sum of feedforward
  % (r_bar) and feedback (z) contributions. 

  r_r = z + r_bar;
  t_zero=crossings;
  
  if strcmp(mode,'script')
      % If the response shows oscillatory tendencies -> simulate for a
      % longer time
      if length(t_zero) >=3
          i_bound = length(t);      % This is the last index in the previous solution
          tspan = [tstop,2*tstop];
          sol = dde23(@dde_eval, Delta_fb, sol, tspan, opt);

          t = (sol.x)';          
          z = (sol.y)';
          t(i_bound)=[];   % delete double entry
          z(i_bound)=[];   % delete double entry
          r_bar = feval(h, lgn_struct, in_struct, t);
          r_r = z+r_bar;
          t_zero = crossings;
      end
  end

% ----------- EVENTS -----------

    function [value,isterminal,direction] = events(t,y,Z) %#ok<INUSD>
        % Nested function to determine events.
        
        if strcmp(form,'vec')
            r_bar = interp1(t_g,r_g,t,'linear');
        elseif strcmp(form,'fnc')
            r_bar = feval(h, lgn_struct, in_struct, t);
        end
        
        value = y + r_bar;
        isterminal = 0; % do not break
        direction = 0; % all zeros are to be computed     
        
    end

% --------- CROSSINGS -------------

    function t_cross = crossings       
        t_cross = sol.xe;
        k = 0;        
        for i = 1:length(t_cross)    
            if t_cross(i) > 1
                k = k+1;
                t_cross(k) = t_cross(i);        
            end    
        end
        t_cross((k+1):i) = [];        
    end

% -------------- History-function ---------------
  function s = dde_history(t)

  % DDE_HISTORY - function which defines the history, in general this is the background firing rate

    s = 0;

  end


% ------------- Evaulation of DDE --------------

  function dzdt = dde_eval(t,z,z_lag)

  % DDE_EVAL - This function defines the model
    
    % Updating the waitbar
    if strcmp('mode','gui')
        waitbar(t/tstop)
    end
    
    % Extracting parameter values    
    
    %Delta_fb = lgn_struct.Delta_fb;
       
    % Setting lagged time
    u = t - Delta_fb;

    % -- Assign value to the retinal input ---------
    
    if strcmp(form,'vec')
        r_bar_lag = interp1(t_g,r_g,u,'linear','extrap');
    elseif strcmp(form,'fnc')
        r_bar_lag = feval(h, lgn_struct, in_struct, u);
    end        

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
  %
  % This part defines the model in use
  % 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          
    % This is the normal model. If w_fbOFFx = w_fbON and l_cOFF = -l_cON, the
    % feedback is linear. 
                
    if (r_bar_lag + z_lag - l_cON) >= 0
        
        dzdt = 1/tau_rc*(w_fbON*(r_bar_lag + z_lag - l_cON) - z);
        
    elseif ( -(r_bar_lag + z_lag - l_cOFF)) >= 0
        
        dzdt = 1/tau_rc*(w_fbOFFx*(l_cOFF - (r_bar_lag + z_lag)) - z);
        
    else
        dzdt = -z/tau_rc;
    end
  end
end

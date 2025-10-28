function [negLL, struct_out] = LogLike_Stage1_dd_UC(starting_values, data_input, other_inputs)
% MODEL IS NOW THE DOUBLE DRIFT MODEL OF CLARK ETC. (see my notes)
% ------------------------------------------------------------------------------------------------------
% 					y(t) = trend(t) + cycle(t)
% 			trend(t) = trend(t-1) + g(t-1) + w(t)
% 					g(t) = g(t-1) + e(t)
%  a(L)*cycle(t) = u(t)
% where y*(t) = trend(t) = y~(t) = cycle
% ------------------------------------------------------------------------------------------------------
% STATE SPACE MODEL:
% 		Observed:	y_t			= D_t + M*alpha_t			+ e_t;		Var(e_t) = H.
% 		State:		alpha_t = C_t + Phi*alpha_t-1	+ R*n_t;	Var(n_t) = Q.
% CALL AS: 
% 		[LogLik, att, Ptt] = kalmanfilter(y, D, M, H, C, Phi, Q, R, a1, P1) 
% Pmean.att = kalmanfilter(Y, Pmean.Dt, Pmean.M, Pmean.H, Ct, Phi, Pmean.Q, R, a00, P00); 
% ******************************************************************************************************
% DIFFUSE PRIOR parsed through other_parameters
% ******************************************************************************************************

% STARTING VALUES
a1 = starting_values(1);
a2 = starting_values(2);

Sigma_ystar		= starting_values(3); 
Sigma_gtrend	= starting_values(4); 
Sigma_ytilde	= starting_values(5); 

% a00_in  = starting_values(7); 
a00_in = other_inputs.a00; 
P00_in = other_inputs.P00;
% OTHER INPUTS
% Sigma_n = other_inputs(1);

% ------------------------------------------------------------------------------------------------------
% MAKE MEASURMENT EQUATION PARAMETERS [D M H]
D = 0;
M = [1 0 1 0];
H = 0;
% MAKE STATE EQUATION PARAMTERS [C Phi R Q]
C = 0;
Phi = [	1 1 0 0;
			  0 1 0 0;
				0 0 a1 a2;
				0 0 1 0];

% SELECTION VECTOR
S = [ 1 0 0;
		  0 1 0;
			0 0 1;
			0 0 0];
		
% VARIANCE OF STATES	
Q = zeros(3,3);
Q(1,1) = (Sigma_ystar)^2;
Q(2,2) = (Sigma_gtrend)^2;
Q(3,3) = (Sigma_ytilde)^2;

% ------------------------------------------------------------------------------------------------------
% INPUT DATA
% ------------------------------------------------------------------------------------------------------
Y = data_input;

% number of states
Ns = length(Phi);

% INITIALIZE THE STATE VECTOR
a00 = zeros(Ns,1);
P00	= zeros(Ns,Ns);  

% set a00(1) at initival value to be estimated.
a00(:) = a00_in;

% FIND THE VARIANCE OF THE STATIONARY AR(4) PART.
SQS= S*Q*S';
phi_2		= Phi(3:end,3:end);
ns = length(phi_2);
% P00_tmp = reshape(inv( eye((Ns-1)^2) - kron(phi_2,phi_2) ) * vec( RQR(2:end,2:end) ),4,4);
P00_tmp = reshape( ( eye((Ns-ns)^2) - kron(phi_2,phi_2) ) \ vec( SQS(1+ns:end,1+ns:end) ),2,2);
P00(ns+1:end,ns+1:end) = P00_tmp;

% MAKE DIFFUSE PRIOR FOR NON-STATIONARY STATE
P00(1:2,1:2) = P00_in*eye(2,2);

% LIKELIHOOD FUNCTION FROM KF RECURSIONS
LL = kalman_filter_loglike(Y, D, M, H, C, Phi, Q, S, a00, P00);

% RETURN NEGATIVE FOR MINIMIZATION
negLL = -LL;

% ------------------------------------------------------------------------------------------------------
% RETURN THE OPTIONAL STRUCTURE OF PARMATERS
% ------------------------------------------------------------------------------------------------------
struct_out.D = D;
struct_out.M = M;
struct_out.H = H;

struct_out.C = C;
struct_out.Phi = Phi;
struct_out.S = S;
struct_out.Q = Q;

struct_out.a00 = a00;
struct_out.P00 = P00;
% loglikelihood, not the negative, ie. -LL>
struct_out.LL	 = LL;	 

% end

















































%%%  EOF
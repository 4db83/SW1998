function [negLL, struct_out] = SW1998_LL_wrapper_G13_ARMA22(starting_values, data_input, other_inputs, KFS_on)
% initial values are = [a(L) = 1:4, Sigma_e a00]
% other inptus are [Sigma_n]
% Sigma_n is fixed at 0.13 , or whatever is put in with Lambd fixed.
% a00(1) is estimated with the rest of the parameters.
% P00 is at stationary part of the model
% ------------------------------------------------------------------------------------------------------
% STATE SPACE MODEL:
% 		Observed:	y_t			= D_t + M*alpha_t			+ e_t;		Var(e_t) = H.
% 		State:		alpha_t = C_t + Phi*alpha_t-1	+ R*n_t;	Var(n_t) = Q.
% CALL AS: 
% 		[LogLik, att, Ptt] = kalmanfilter(y, D, M, H, C, Phi, Q, R, a1, P1) 
% Pmean.att = kalmanfilter(Y, Pmean.Dt, Pmean.M, Pmean.H, Ct, Phi, Pmean.Q, R, a00, P00); 
% ******************************************************************************************************
% NOTE: P00(1,1) IS THE SAME AS IN STOCK AND WATSON LNAIRC.PRC LNAIR1.PRC FILES
% ******************************************************************************************************

SetDefaultValue(4,'KFS_on',0);

% STARTING VALUES
a1 = starting_values(1);
a2 = starting_values(2);
a3 = starting_values(3);

b1 = starting_values(4);
b2 = starting_values(5);

Sigma_e	= starting_values(6); 
a00_in  = starting_values(end); 
% Sigma_n = startingvalues(6); 

% OTHER INPUTS (sigma_zeta = sigma_\Delta_\beta = Lambda/T*LR-stdev(eps(t)).
Sigma_z = other_inputs(1);

% ------------------------------------------------------------------------------------------------------
% MAKE MEASURMENT EQUATION PARAMETERS [D M H]
D = 0;
M = [1 1 b1 b2];
H = 0;
% MAKE STATE EQUATION PARAMTERS [C Phi R Q]
C = 0;
Phi = [	1 0 0 0;
				0 a1 a2 a3;
			  0 1 0 0;
				0 0 1 0];

% SELECTION VECTOR
S = [ 1 0;
		  0 1;
			0 0;
			0 0];
% VARIANCE OF STATES	
Q = zeros(2,2);
Q(1,1) = (Sigma_z)^2;
Q(2,2) = (Sigma_e)^2;

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
a00(1) = a00_in;

% FIND THE VARIANCE OF THE STATIONARY AR(4) PART.
SQS		= S*Q*S';
phi_2	= Phi(2:end,2:end);

% P00_tmp = reshape(inv( eye((Ns-1)^2) - kron(phi_2,phi_2) ) * vec( RQR(2:end,2:end) ),4,4);
P00_tmp = reshape( ( eye((Ns-1)^2) - kron(phi_2,phi_2) ) \ vec( SQS(2:end,2:end) ),3,3);
P00(2:end,2:end) = P00_tmp;
% ---------------------------------------------------------------------------------------------
% NOTE HERE THAT AS IN STOCK AND WATSON, IN MUE SECOND STAGE KF, P00 MATRIX HAS P00(1,1) = 0.
% ---------------------------------------------------------------------------------------------

% LIKELIHOOD FUNCTION FROM KF RECURSIONS
LL = kalman_filter_loglike(Y, D, M, H, C, Phi, Q, S, a00, P00)

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

% Kalman Filter and Smooth the results now
if KFS_on
	struct_out.KFS = kalman_filter_smoother(Y, D, M, H, C, Phi, Q, S, a00, P00, 0);
end



%%  EOF







function [negLL, struct_out] = LogLike_MNZ0(starting_values, data_input, other_inputs, KFS_on)
% ------------------------------------------------------------------------------------------------------
% MODEL IS:
% ------------------------------------------------------------------------------------------------------
%						y(t) = trend(t) + cycle(t)
%				trend(t) = mu + trend(t-1) + w(t)
%	 a(L)*cycle(t) = u(t)
% where y*(t) = trend(t) = y~(t) = cycle
% ------------------------------------------------------------------------------------------------------
% STATE SPACE MODEL:
% 		Observed:	Y(t)			= D(t) + M*alpha(t)			+ e(t);		Var(e_t) = H.
% 		State:		alpha(t)	= C(t) + Phi*alpha(t-1)	+ S*n(t);	Var(n_t) = Q.
% CALL AS: 
% 		[LogLik, att, Ptt] = kalmanfilter(y, D, M, H, C, Phi, Q, R, a1, P1) 
% Pmean.att = kalmanfilter(Y, Pmean.Dt, Pmean.M, Pmean.H, Ct, Phi, Pmean.Q, R, a00, P00); 
% ******************************************************************************************************
% DIFFUSE PRIOR parsed through other_parameters
% ******************************************************************************************************

SetDefaultValue(4,'KFS_on',0);

% STARTING VALUES
a1 = starting_values(1);
a2 = starting_values(2);
mu = starting_values(3);

% STANDARD DEVIATIONS
sigma_ystar = starting_values(4); 
sigma_ytild	= starting_values(5); 

% OTHER INPUTS
% prior on states
a00_in	= other_inputs.a00; 
P00_in	= other_inputs.P00;
% correlation parameter fixed at some other value, could be zero, or in grid search
rho			= other_inputs.rho;

% ------------------------------------------------------------------------------------------------------
% MAKE MEASURMENT EQUATION PARAMETERS [D M H]
D = 0;
M = [1 1 0];
H = 0;
% MAKE STATE EQUATION PARAMTERS [C Phi R Q]
C = [mu; 0; 0];
Phi = [	1  0  0;
			  0 a1 a2;
				0  1  0];

% SELECTION VECTOR
S = [ 1 0;
		  0 1;
			0 0];
		
nQ = size(S,2);		
		
% VARIANCE OF STATES	
Q = zeros(nQ);
Q(1,1) = (sigma_ystar)^2;
Q(2,2) = (sigma_ytild)^2;
Q(1,2) = sigma_ystar*sigma_ytild*rho ;
Q(2,1) = sigma_ystar*sigma_ytild*rho ;

% ------------------------------------------------------------------------------------------------------
% INPUT DATA
% ------------------------------------------------------------------------------------------------------
Y = data_input;

% number of states
nPhi = length(Phi);

% INITIALIZE THE STATE VECTOR
a00 = zeros(nPhi,1);
P00	= zeros(nPhi);  

% set a00(1) at initival value to be estimated.
a00(:) = a00_in;

% Precompute SQS to avoid having to compute inside the loop
SQS		= S*Q*S';

% FIND THE VARIANCE OF THE STATIONARY AR PART.
nI1		= 1;													% number of I(1) states
phi_2	= Phi(nI1+1:end,nI1+1:end);		% I(0) State vector Phi
nI0		= size(phi_2,1);							% % number of I(0) states

% P00_tmp = reshape(inv( eye((Ns-1)^2) - kron(phi_2,phi_2) ) * vec( RQR(2:end,2:end) ),4,4);
P00_tmp = reshape( ( eye(nI0^2) - kron(phi_2,phi_2) ) \ vec( SQS(nI1+1:end,nI1+1:end) ), nI0, nI0);
P00(nI1+1:end,nI1+1:end) = P00_tmp;

% MAKE DIFFUSE PRIOR FOR NON-STATIONARY STATE
P00(1:nI1,1:nI1) = P00_in*eye(nI1);

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

% Kalman Filter and Smooth the results now
if KFS_on
	struct_out.KFS = kalman_filter_smoother(Y, D, M, H, C, Phi, Q, S, a00, P00);
end

% end


















































%%%  EOF
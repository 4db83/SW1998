function [negLL, struct_out] = LogLike_Stage1_DB_SSF(starting_values, data_input, other_inputs)
% ---------------------------------------------------------------------------------------------
% MODEL IS STAGE ONE MODEL OF HLW BUT IN MY STATE SPACE FORM:
% ---------------------------------------------------------------------------------------------
% STATE SPACE MODEL:
% 		Observed:	Y(t)			= D(t) + M*alpha(t)			+ e(t);		Var(e_t) = H.
% 		State:		alpha(t)	= C(t) + Phi*alpha(t-1)	+ S*n(t);	Var(n_t) = Q.
% CALL AS: 
% 		[LogLik, att, Ptt] = kalmanfilter(y, D, M, H, C, Phi, Q, R, a1, P1) 
% Pmean.att = kalmanfilter(Y, Pmean.Dt, Pmean.M, Pmean.H, Ct, Phi, Pmean.Q, R, a00, P00); 
% *********************************************************************************************
% DIFFUSE PRIOR parsed through other_parameters
% *********************************************************************************************

% ---------------------------------------------------------------------------------------------
% INPUT PARAMETERS
% ---------------------------------------------------------------------------------------------

% STARTING VALUES
AR2		= starting_values(1:2);
b_pi	= starting_values(3);
b_y		= starting_values(4);

mu		= starting_values(5);
pi_bar = 0;

% STANDARD DEVIATIONS
sigma_ystar = starting_values(end-2); 
sigma_pi		= starting_values(end-1); 
sigma_ytild	= starting_values(end); 

% CORRELATION PARAMETER TO BE ESTIMATED
% rho = starting_values(6);

% OTHER INPUTS
% prior on states
a00_in	= other_inputs.a00; 
P00_in	= other_inputs.P00;

% Sigma_n = other_inputs(1);
% rho			= other_inputs.rho;

% ---------------------------------------------------------------------------------------------
% MAKE KALMAN FILTER PARAMETER MATRICES
% ---------------------------------------------------------------------------------------------
nS	= 7;		% Number of rows in state vector alpha
nQ	= 3;		% Number of shocks in state vector alpha
nI1	= 5;		% Number of I(1) variables in state vector alpha

% MAKE MEASURMENT EQUATION PARAMETERS [D M H]
D = 0;

M = zeros(2,nS);
M(1,[1 6])	= 1;
M(2,2)			=	1;

H = 0;
% MAKE STATE EQUATION PARAMTERS [C Phi R Q]
C = zeros(nS,1);
C(1) = mu;
C(2) = pi_bar;

Phi = [	1  0 	0  0  0  0   0;
				0  0	0	 0	0	 b_y 0;
				0  1 	0  0  0  0   0;
				0  0 	1  0  0  0   0;
				0  0 	0  1  0  0   0;
				0  0 	0  0  0  0 	 0;
				0  0 	0  0  0  1   0];
			
Phi(2,[2])	 = b_pi;
Phi(2,[3:5]) = (1-b_pi)/3;

Phi(6,[6:7]) = AR2';

% SELECTION VECTOR
S = zeros(nS,nQ);
S(1,1) = 1;
S(2,2) = 1;
S(6,3) = 1;

		
% VARIANCE OF STATES	
Q = zeros(nQ);
Q(1,1) = (sigma_ystar)^2;
Q(2,2) = (sigma_pi)^2;
Q(3,3) = (sigma_ytild)^2;

% ---------------------------------------------------------------------------------------------
% INPUT DATA
% ---------------------------------------------------------------------------------------------
Y = data_input;

% INITIALIZE THE STATE VECTOR
a00 = zeros(nS,1);
P00	= zeros(nS);  

% set a00(1) at initival value to be estimated.
a00(:) = a00_in;

% Precompute SQS to avoid having to compute inside the loop
SQS		= S*Q*S';

% FIND THE VARIANCE OF THE STATIONARY AR PART.
phi_2	= Phi(nI1+1:end,nI1+1:end);		% I(0) State vector Phi
nI0		= size(phi_2,1);							% % number of I(0) states

%%% UNCOMMNET BELOW TO USE DIFFUSE IF EIG PHI_2 > .99
% % if max(eig(phi_2))>.99
% % 	P00_tmp = P00_in*eye(size(phi_2));
% % else
% % 	P00_tmp = reshape( ( eye(nI0^2) - kron(phi_2,phi_2) ) \ vec( SQS(nI1+1:end,nI1+1:end) ), nI0, nI0);
% % end

%%% UNCOMMNET BELOW TO USE STANDARD FORM
P00_tmp = reshape( ( eye(nI0^2) - kron(phi_2,phi_2) ) \ vec( SQS(nI1+1:end,nI1+1:end) ), nI0, nI0);
%%%%%	% P00_tmp = reshape(inv( eye((Ns-1)^2) - kron(phi_2,phi_2) ) * vec( RQR(2:end,2:end) ),4,4);

P00(nI1+1:end,nI1+1:end) = P00_tmp;

% MAKE DIFFUSE PRIOR FOR NON-STATIONARY STATE
P00(1:nI1,1:nI1) = P00_in*eye(nI1);

% LIKELIHOOD FUNCTION FROM KF RECURSIONS
LL = kalman_filter_loglike(Y, D, M, H, C, Phi, Q, S, a00, P00);

% RETURN NEGATIVE FOR MINIMIZATION
negLL = -LL;

% ---------------------------------------------------------------------------------------------
% RETURN THE OPTIONAL STRUCTURE OF PARMATERS
% ---------------------------------------------------------------------------------------------
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
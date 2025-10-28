function [negLL, struct_out] = LogLike_Stage1_HLW_SSF_g_in_trend(starting_values, data_input, other_inputs, KFS_on)
% Replicates Stage 1 model results of HLW(2017), using exactly their data as input.
% Stage 1 SSF of HWL, but now with g in the RW specification for the trend, so C = [g;0;0] 
% rather than 0 as before in LogLike_Stage1_HWL_SSF_detrend.
% ---------------------------------------------------------------------------------------------
% Their SSF is:
% ---------------------------------------------------------------------------------------------
% 	Y(t) = A*X(t) + H*Xi(t) + v(t), Var(v) = R
%  Xi(t) = F*Xi(t-1) + S*e(t),			Var(e) = Q
% ---------------------------------------------------------------------------------------------
% MY STATE SPACE FORM IS:
% 		Observed:	Y(t)			= D(t) + M*alpha(t)			+ e(t);		Var(e_t) = H.
% 		State:		alpha(t)	= C(t) + Phi*alpha(t-1)	+ S*n(t);	Var(n_t) = Q.
% ---------------------------------------------------------------------------------------------
% CALL AS: 
% ---------------------------------------------------------------------------------------------
% 		[LogLik, att, Ptt] = kalmanfilter(y, D, M, H, C, Phi, Q, R, a1, P1) 
% Pmean.att = kalmanfilter(Y, Pmean.Dt, Pmean.M, Pmean.H, Ct, Phi, Pmean.Q, R, a00, P00); 
% ---------------------------------------------------------------------------------------------

SetDefaultValue(4,'KFS_on',0);

% STARTING VALUES
% mean parameters
a = starting_values(1:2)';
b = starting_values(3:4)';
g = starting_values(5);

% variances
s2_ytld		= starting_values(end-2)^2;
s2_pi			= starting_values(end-1)^2;
s2_ystr		=	starting_values(end)^2;

% OTHER INPUTS
% initial state mean and variance
a00 = other_inputs.a00;
P00 = other_inputs.P00;
% g =	other_inputs.g;

% ---------------------------------------------------------------------------------------------
% INPUT DATA
% ---------------------------------------------------------------------------------------------
Y = data_input(1:2,:);
X = data_input(3:end,:);
TT = length(Y);

% ---------------------------------------------------------------------------------------------
% MAKE MEASURMENT EQUATION PARAMETERS [D M H]
% ---------------------------------------------------------------------------------------------
% MAKE THEIR A FOR EXOGENEOUS VARIABLES
A = [	a(1:2) 0 0	;
			b(2) 0 b(1) ( 1-b(1) )];

% MAKE Dt
D = A*X;

% THEIR H. NOTE THE RESCALING BY 4!, ie., -a_r/2 * 4 --> -a_r*2 for column 4:5 entries
M = [ 1 -a(1:2) ;
			0 -b(2) 0 ];

% VARIANCE OF MEASUREMENT
H	= zeros(2,2);
H(1,1) = s2_ytld;
H(2,2) = s2_pi;

% ---------------------------------------------------------------------------------------------
% MAKE STATE EQUATION PARAMTERS [C Phi R Q]
% ---------------------------------------------------------------------------------------------
C		= [g;0;0];

Phi = [	1 0 0;
				1 0 0;
			  0 1 0];

% NUMBER OF STATES
Ns = length(Phi);

% VARIANCE OF STATES	

% SELECTION MATRIX S: set to identity when using the hwl specifiacition with lambdas 
S = eye(Ns);
% VARIANCE OF STATES
Q = zeros(Ns);

Q(1,1) = s2_ystr;

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
% loglikelihood LL, not the negative, ie. -LL 
struct_out.LL	 = LL;	 

% Kalman Filter and Smooth the results now
if KFS_on
	struct_out.KFS = kalman_filter_smoother(Y, D, M, H, C, Phi, Q, S, a00, P00);
end










































%%%  EOF
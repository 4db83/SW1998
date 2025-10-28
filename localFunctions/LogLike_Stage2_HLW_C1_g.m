function [negLL, struct_out] = LogLike_Stage2_HLW_C1_g(starting_values, data_input, other_inputs, KFS_on)
% STAGE 2 MODEL of HLW(2017), using exactly their data and priors, with TWO LAGS of g(t) in Xi(t).
% MLE on sigma_g, rather than fixing it at MUE Stage 1 estimate.
% ------------------------------------------------------------------------------------------------------
%   STAGE 2 MODEL M(1)MLE.g: with MLE on sigma_g and INTERCEPT a_0 also estimate
% ------------------------------------------------------------------------------------------------------
%   ay(L)y~(t) = a_0 + a_r(L)[r(t)-4*g(t)], with 
%			  Dy*(t) = g(t-1) + e^y*(t).		[rather than with g(t-2)!, correct S specification]
% ------------------------------------------------------------------------------------------------------
%		Params to estimate:  [a_y1, a_y2, a_r, a_0,   b_pi, b_y, sig_y~, sig_pi, sig_y*, sig_g]
% 											 [ 1  ,   2 ,  3 ,  4 ,    6  ,  7 ,   8   ,   9   ,   10  ,  11  ]
%		This is the minimum misspecified (restricted) Stage 2 model to operationalize MUE for S2.
%		State Vector Xi = [y*(t) y*(t-1) y*(t-2) g(t-1) g(t-2)].
% ------------------------------------------------------------------------------------------------------
% Their SSF is:
% ------------------------------------------------------------------------------------------------------
% 	Y(t) = A*X(t) + H*Xi(t) + v(t), Var(v) = R
%  Xi(t) = F*Xi(t-1) + S*e(t),			Var(e) = Q
% ------------------------------------------------------------------------------------------------------
% MY STATE SPACE FORM IS:
% 		Observed:	Y(t)			= D(t) + M*alpha(t)			+ e(t);		Var(e_t) = H.
% 		State:		alpha(t)	= C(t) + Phi*alpha(t-1)	+ S*n(t);	Var(n_t) = Q.
% ------------------------------------------------------------------------------------------------------
% CALL AS: 
% ------------------------------------------------------------------------------------------------------
% 		[LogLik, att, Ptt] = kalmanfilter(y, D, M, H, C, Phi, Q, R, a1, P1) 
% Pmean.att = kalmanfilter(Y, Pmean.Dt, Pmean.M, Pmean.H, Ct, Phi, Pmean.Q, R, a00, P00); 
% ******************************************************************************************************
% PRIORs are parsed through other_inputs
% FULL parameter vector: [a_y1, a_y2, a_r,	a_0, a_g, b_pi, b_y, sig_y~, sig_pi, sig_y*, sig_g]
% 										   [ 1  ,   2 ,  3 ,	 4 ,  5 ,  6  ,  7 ,   8   ,   9   ,   10  ,  11  ] 
% ******************************************************************************************************

SetDefaultValue(4,'KFS_on',0);

nS	= 5;		% Number of rows in state vector alpha
nQ	= 2;		% Number of shocks in state vector alpha
nH  = 2;		% Number of shocks in measurement vector
nI1	= 5;		% Number of I(1) variables in state vector alpha

% OTHER INPUTS
% initial state mean and variance
a00 = other_inputs.a00;
P00 = other_inputs.P00;
Lambda_g = other_inputs.Lambda_g;
% Lambda_z = other_inputs.Lambda_z;

% STARTING VALUES
% mean parameters
a_y1	= starting_values(1);
a_y2	= starting_values(2);
a_r		= starting_values(3);
a_0	  = starting_values(4);							% intercept term included now

b_pi	= starting_values(5);
b_y		= starting_values(6);

% VARIANCES
s2_ytld	= starting_values(end-3)^2;                
s2_pi		= starting_values(end-2)^2;                
s2_ystr	=	starting_values(end-1)^2;                
s2_g		= starting_values(end  )^2;			% sigma2_g 

% ------------------------------------------------------------------------------------------------------
% INPUT DATA
% ------------------------------------------------------------------------------------------------------
Y = data_input(1:2,:);
X = data_input(3:end,:);

% ------------------------------------------------------------------------------------------------------
% MAKE MEASURMENT EQUATION PARAMETERS [D M H]
% ------------------------------------------------------------------------------------------------------
% MAKE THEIR A FOR EXOGENEOUS VARIABLES		(add a_0 for the constant)
A = [	a_y1 a_y2 a_r/2 a_r/2  0		 0        a_0 ;
			b_y  0    0     0      b_pi (1-b_pi)  0		];

% MAKE Dt
D = A*X;

% THEIR H. NOTE THE RESCALING BY 4 annualized as in HLW code!, ie., -a_r/2 * 4 --> -a_r*2 for column 4:5 entries. Here we impose the same coefficient and use the core ystar equation as from the model 3 trick
M = [ 1 -a_y1 -a_y2  -4*a_r/2  -4*a_r/2  ;
			0 -b_y   0      0					0				];
			
% normally we have this for the Stage 3 model where the coefficients on g(t) and z(t) are restricted to be -a_r/2 (and annualized for g(t) trend growth)
% M = [ 1 -a_y1 -a_y2  -4*a_r/2  -4*a_r/2  -a_r/2  -a_r/2;
% 			0 -b_y   0      0					0					0				0			];

% VARIANCE OF MEASUREMENT
H	= zeros(nH);
H(1,1) = s2_ytld;
H(2,2) = s2_pi;

% ------------------------------------------------------------------------------------------------------
% MAKE STATE EQUATION PARAMTERS [C Phi R Q]
% ------------------------------------------------------------------------------------------------------
C		= 0;

Phi = [1 0 0 1 0;
			 1 0 0 0 0;
			 0 1 0 0 0;
			 0 0 0 1 0;
			 0 0 0 1 0];

% SELECTION MATRIX S: set to identity when using the hlw specifiacition with lambdas 
S = [1 1; 
		 0 0; 
		 0 0; 
		 0 1; 
		 0 0];

% VARIANCE OF STATES
Q = zeros(nQ);
Q(1,1) = s2_ystr;
Q(2,2) = s2_g;

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
% loglikelihood LL, not the negative, ie. -LL 
struct_out.LL	 = LL;	 

% Kalman Filter and Smooth the results now
if KFS_on
	struct_out.KFS = kalman_filter_smoother(Y, D, M, H, C, Phi, Q, S, a00, P00);
end


% End of function
end







%%%  EOF
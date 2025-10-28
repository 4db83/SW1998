function [negLL, struct_out] = LogLike_Stage1_HWL_SSF_detrend(starting_values, data_input, other_inputs, AddTrend)
% Replicates Stage 1 model results of HLW(2017), using exactly their data as input.
% Stage 1 SSF of HWL with 'artificial' detrending of y.
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


SetDefaultValue(4,'AddTrend',0);

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

% what they call 'Make the data stationary', by detrending on a time trend
Y(1,:) = Y(1,:) - g*(1:TT);
X(1,:) = X(1,:) - g*(0:TT-1);
X(2,:) = X(2,:) - g*(-1:TT-2);

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
C		= 0;
Phi = [	1 0 0;
				1 0 0;
			  0 1 0];

% NUMBER OF STATES
Ns = length(Phi);

% VARIANCE OF STATES	

% SELECTION MATRIX S: set to identity when using the hwl specifiacition with lambdas 
% S = eye(Ns);
S = [1;0;0];
% VARIANCE OF STATES
% Q = zeros(Ns);
% Q(1,1) = s2_ystr;
Q = s2_ystr;

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

% theta out.
struct_out.theta_out = [a';b';g;sqrt(s2_ytld);sqrt(s2_pi);sqrt(s2_ystr)];

% NOW ADD THE TREND BACK IN AS THEY DO IN KALMAN.STATES.WRAPPER.R
% if (stage == 1) {
%       states$filtered$xi.tt <- states$filtered$xi.tt + cbind(1:T,0:(T-1),-1:(T-2)) * parameters[5]
%       states$smoothed$xi.tT <- states$smoothed$xi.tT + cbind(1:T,0:(T-1),-1:(T-2)) * parameters[5]
%   }
if AddTrend
	KFout = kalman_filter_smoother(Y, D, M, H, C, Phi, Q, S, a00, P00);
	% ADD TREND BACK 
	KFout_plus_trend = KFout;
	KFout_plus_trend.att = KFout.att + g*[(1:TT)' (0:TT-1)' (-1:TT-2)'];
	KFout_plus_trend.atT = KFout.atT + g*[(1:TT)' (0:TT-1)' (-1:TT-2)'];
	% now return to structure
	struct_out.KFS = KFout;
	struct_out.KFS_plus_trend = KFout_plus_trend;
end

%% end










































%%%  EOF
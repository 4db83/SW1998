function [negLL, struct_out] = LogLike_Stage2_HWL_SSF(starting_values, data_input, other_inputs, KFS_on)
% Stage 2 model results of HLW(2017), using exactly their data as input.
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
% DIFFUSE PRIOR parsed through other_inputs
% ******************************************************************************************************

SetDefaultValue(4,'KFS_on',0);

nS	= 4;		% Number of rows in state vector alpha
nQ	= 2;		% Number of shocks in state vector alpha
nI1	= 1;		% Number of I(1) variables in state vector alpha

% OTHER INPUTS
% initial state mean and variance
a00 = other_inputs.a00;
P00 = other_inputs.P00;
Lambda_g = other_inputs.Lambda_g;

% STARTING VALUES
% mean parameters
a_y1	= starting_values(1);
a_y2	= starting_values(2);
a_r		= starting_values(3);
a_0		= starting_values(4);
a_g		= starting_values(5);

b_pi	= starting_values(6);
b_y		= starting_values(7);

% variances
s2_ytld		= starting_values(end-2)^2;
s2_pi			= starting_values(end-1)^2;
s2_ystr		=	starting_values(end)^2;
s2_g			= ( Lambda_g*starting_values(end) )^2;

% ------------------------------------------------------------------------------------------------------
% INPUT DATA
% ------------------------------------------------------------------------------------------------------
Y = data_input(1:2,:);
X = data_input(3:end,:);
TT = length(Y);

% ------------------------------------------------------------------------------------------------------
% MAKE MEASURMENT EQUATION PARAMETERS [D M H]
% ------------------------------------------------------------------------------------------------------
% MAKE THEIR A FOR EXOGENEOUS VARIABLES
A = [	a_y1 a_y2 a_r/2 a_r/2  0    0				 a_0;
			b_y  0    0     0      b_pi (1-b_pi) 0];

% MAKE Dt
D = A*X;

% THEIR H. NOTE THE RESCALING BY 4!, ie., -a_r/2 * 4 --> -a_r*2 for column 4:5 entries
M = [ 1 -a_y1 -a_y2  a_g ;
			0 -b_y   0     0];

% VARIANCE OF MEASUREMENT
H	= zeros(2,2);
H(1,1) = s2_ytld;
H(2,2) = s2_pi;

% ------------------------------------------------------------------------------------------------------
% MAKE STATE EQUATION PARAMTERS [C Phi R Q]
% ------------------------------------------------------------------------------------------------------
C		= 0;

Phi = [	1 0 0 1;
				1 0 0 0;
			  0 1 0 0;
				0 0 0 1];

% SELECTION MATRIX S: set to identity when using the hlw specifiacition with lambdas 
S = zeros(nS,nQ);
S(1,1)	= 1;
S(nS,2) = 1;

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
	struct_out.KFS = kalman_filter_smoother(Y, D, M, H, C, Phi, Q, S, a00, P00, 0);
end

% end

% theta out.
% struct_out.theta_out = [a';b';g;sqrt(s2_ytld);sqrt(s2_pi);sqrt(s2_ystr)];

% NOW ADD THE TREND BACK IN AS THEY DO IN KALMAN.STATES.WRAPPER.R
% if (stage == 1) {
%       states$filtered$xi.tt <- states$filtered$xi.tt + cbind(1:T,0:(T-1),-1:(T-2)) * parameters[5]
%       states$smoothed$xi.tT <- states$smoothed$xi.tT + cbind(1:T,0:(T-1),-1:(T-2)) * parameters[5]
%   }
% % if AddTrend
% % 	KFout = kalman_filter_smoother(Y, D, M, H, C, Phi, Q, S, a00, P00);
% % 	% ADD TREND BACK 
% % 	KFout_plus_trend = KFout;
% % 	KFout_plus_trend.att = KFout.att + g*[(1:TT)' (0:TT-1)' (-1:TT-2)'];
% % 	KFout_plus_trend.atT = KFout.atT + g*[(1:TT)' (0:TT-1)' (-1:TT-2)'];
% % 	% now return to structure
% % 	struct_out.KFS = KFout;
% % 	struct_out.KFS_plus_trend = KFout_plus_trend;
% % end

%% end










































%%%  EOF
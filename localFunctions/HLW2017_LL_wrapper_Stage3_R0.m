function [negLL, struct_out]= HLW2017_LL_wrapper_Stage3_R0(starting_values, data_input, other_inputs)
% initial values are = [a(L) = 1:4, Sigma_e a00]
% other inptus are [Sigma_n]
% Sigma_n is fixed at 0.13 , or whatever is put in with Lambd fixed.
% a00(1) is estimated with the rest of the parameters.
% P00 is at stationary part of the model
% ------------------------------------------------------------------------------------------------------

% STARTING VALUES
% mean parameters
a = starting_values(1:3)';
b = starting_values(4:5)';
% variances
s2_ytld		= starting_values(6)^2;
s2_pi			= starting_values(7)^2;
s2_ystr		=	starting_values(8)^2;
s2_g			= starting_values(9)^2;
s2_z			= starting_values(10)^2;

% OTHER INPUTS
% predetermined lambda values
% LAMBDA_G = other_inputs.LAMBDA_G;
% LAMBDA_Z = other_inputs.LAMBDA_Z;
% initial state mean and variance
a00 = other_inputs.a00;
P00 = other_inputs.P00;

% ------------------------------------------------------------------------------------------------------
% INPUT DATA
% ------------------------------------------------------------------------------------------------------
Y = data_input(1:2,:);
X = data_input(3:end,:);

% ------------------------------------------------------------------------------------------------------
% MAKE MEASURMENT EQUATION PARAMETERS [D M H]
% ------------------------------------------------------------------------------------------------------
% MAKE THEIR A FOR EXOGENEOUS VARIABLES
A = [	a(1:2) a(3)/2 a(3)/2 0 0	;
			b(2) 0 0 0 b(1) ( 1-b(1) )];

% MAKE Dt
D = A*X;

% THEIR H. NOTE THE RESCALING BY 4!, ie., -a_r/2 * 4 --> -a_r*2 for column 4:5 entries
M = [ 1 -a(1:2) -a(3)*2 -a(3)*2 -a(3)/2 -a(3)/2;
			0 -b(2) 0 0 0 0 0];

% VARIANCE OF MEASUREMENT
H	= zeros(2,2);
H(1,1) = s2_ytld;
H(2,2) = s2_pi;

% ------------------------------------------------------------------------------------------------------
% MAKE STATE EQUATION PARAMTERS [C Phi R Q]
% ------------------------------------------------------------------------------------------------------
C		= 0;
Phi = [	1 0 0 1 0 0 0;
				1 0 0 0 0 0 0;
			  0 1 0 0 0 0 0;
				0 0 0 1 0 0 0;
				0 0 0 1 0 0 0;
				0 0 0 0 0 1 0;
				0 0 0 0 0 1 0];

% NUMBER OF STATES
Ns = length(Phi);

% VARIANCE OF STATES	

% SELECTION MATRIX S: set to identity when using the hwl specifiacition with lambdas 
% S = eye(Ns);
S = [	1 1 0;
			0 0 0;
			0 0 0;
			0 1 0;
			0 0 0;
			0 0 1;
			0 0 0];
			
% VARIANCE OF STATES
Q = zeros(3);

Q(1,1) = s2_ystr;
Q(2,2) = s2_g;
Q(3,3) = s2_z;

%%% INITIALIZE THE STATE VECTOR
% these are taken from their P00 and a00 values;
% % a00 = zeros(Ns,1);
% % P00	= zeros(Ns,Ns);  
% % 
% % % set a00(1) at initival value to be estimated.
% % a00(1) = a00_in;
% % 
% % % FIND THE VARIANCE OF THE STATIONARY AR(4) PART.
% % RQR= R*Q*R';
% % phi_2		= Phi(2:end,2:end);
% % % P00_tmp = reshape(inv( eye((Ns-1)^2) - kron(phi_2,phi_2) ) * vec( RQR(2:end,2:end) ),4,4);
% % P00_tmp = reshape( ( eye((Ns-1)^2) - kron(phi_2,phi_2) ) \ vec( RQR(2:end,2:end) ),4,4);
% % P00(2:end,2:end) = P00_tmp;

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

% end

















































%%%  EOF
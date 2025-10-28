function [LAMBDA_Z_db, Lambda_chow, Aout] = call_MUE_HLW_stage2_M0g(theta, S2object, YY, XX)
% this is a modified version of the call_MUE_HLW_stage2.m file because we now need to construct
% GY(t) = y~(t)-a_1y~(t)-a_1y~(t) - ar(r(t-1)+r(t-2)/2 - 4*(g(t-1)+g(t-2));
% get the theta parameter input vector with right dimension
theta4_MUE2 = theta(1:3); % get the order right [ay1,ay2,ar]

% MAKE THE Y AND X INPUT VARIABLES of the correct model now
% YY	  = [lnGDP INFL]; 
% [GDP(t-1) GDP(t-2) r(t-1) r(t-2) pi(t-1) mean(pi(t-2:t-4))]
% XX	= [ GPD_2_lags RR_2_lags INFL_4_lags(:,1) mean(INFL_4_lags(:,2:4),2) ];

YX_db = make_yx_4_stage2_M0g(YY, XX, S2object);
yy = YX_db.tT(:,1);			% y~(t) (T x 1) (y~(t) = cycle variable)
xx = YX_db.tT(:,2:end); % [y~(t-1) y~(t-2) rr_2-gg_2] (T x 3)
												%	or [y~(t-1) y~(t-2) {r(t-1)-4*g(t-1)+r(t-2)-4*g(t-2)}/2 ] (T x 3) 

% CALL TO MUE PROCEDURE TO COMPUTE LAMBDA_Z
[LAMBDA_Z_db, Lambda_chow, Stats_ols, Stats_chow, Aout]	= MUE_HLW_Stage2_M0g(yy, xx, theta4_MUE2);

% % -------------------------------------------------------------------------------------------------------
% % CALL: make_yx_4_stage2_M0g_unrestricted to return XX = [y~(t-1) y~(t-1) r(t-1) r(t-2) g(t-1) g(t-1) 1]
% % -------------------------------------------------------------------------------------------------------
% YX_ur = make_yx_4_stage2_M0g_unrestricted(YY, XX, S2object);
% yy_ur = YX_ur.tT(:,1);	% y~(t)
% x	= YX_ur.tT(:,2:end);	% [y~(t-1) y~(t-1) r(t-1) r(t-2) g(t-1) g(t-2) 1]
%  
% xx_ur = [x(:,1:2) x(:,3) x(:,5) x(:,end)];
%  
% [L1,~,~,~,A1]	= MUE_HLW_Stage2_M0g(yy_ur, xx_ur);
%  
% plot(Aout.Fstat_ols);
% hold on;
% plot(A1.Fstat_ols,'--');
% hold off;
% 
% printstructs(L1, LAMBDA_Z_db)

% adding y and x as well as parameter input vector
Aout.yy = yy;
Aout.xx = xx;
Aout.YX = YX_db;
Aout.stats_ols		= Stats_ols;
Aout.stats_chow		= Stats_chow;
Aout.theta4_MUE2	= theta4_MUE2;


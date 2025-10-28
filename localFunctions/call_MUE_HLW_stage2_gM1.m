function [LAMBDA_Z_db, Lambda_chow, Aout] = call_MUE_HLW_stage2_gM1(theta, S2object, YY, XX)
% SAME AS call_MUE_HLW_stage2_gM0, but now includes an intercept term as well
% this is a modified version of the call_MUE_HLW_stage2.m file because we now need to construct
% GY(t) = y~(t)-a_1y~(t)-a_1y~(t) - ar(r(t-1)+r(t-2)/2 - 4*(g(t-1)+g(t-2)) - a0;
% get the theta parameter input vector with right dimension
theta4_MUE2 = theta(1:4); % get the order right [ay1,ay2,ar,a0]

% MAKE THE Y AND X INPUT VARIABLES AS IN HLW(2017)
% YY	  = [lnGDP INFL]; 
% [GDP(t-1) GDP(t-2) r(t-1) r(t-2) pi(t-1) mean(pi(t-2:t-4))]
% XX	= [ GPD_2_lags RR_2_lags INFL_4_lags(:,1) mean(INFL_4_lags(:,2:4),2) ];

YX_db = make_yx_4_stage2_gM0(YY, XX, S2object);
yy = YX_db.tT(:,1);			% y~(t) (T x 1)
% xx = YX_db.tT(:,2:end); % [y~(t-1) y~(t-2) r_2-gg_2] (T x 3)

% now add a vector of ones
xx = [YX_db.tT(:,2:end) ones(length(YX_db.tT(:,2:end)),1) ]; % [y~(t-1) y~(t-2) r_2-gg_2 1] (T x 4)

% GY = yy - xx*theta4_MUE2;
% MUE_HLW_Stage1(GY)

% CALL TO MUE PROCEDURE TO COMPUTE LAMBDA_Z
[LAMBDA_Z_db, Lambda_chow, Stats_ols, Stats_chow, Aout]	= MUE_HLW_Stage2_gM0(yy, xx, theta4_MUE2);

% adding y and x as well as parameter input vector
Aout.yy = yy;
Aout.xx = xx;
Aout.YX = YX_db;
Aout.stats_ols		= Stats_ols;
Aout.stats_chow		= Stats_chow;
Aout.theta4_MUE2	= theta4_MUE2;


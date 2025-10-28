function Hout = setoutsideTicks(tickSize, doubleTicks, LabeloffSet, currentHandle)
%{F: Sets Y axis ticks on the outside when Yaxis labels are printed on boths sides
% ---------------------------------------------------------------------------------------------
% CALL AS: 
%		setoutsideTicks										or as
%		setoutsideTicks(tickSize)
% 	setoutsideTicks(tickSize, currentHandle)
%										
%	INPUT:	
%		ticksize = ticksize adjustment shrinks the ticksize to 3/4 of original length.
%		currentHandle = handle to current figure (Defaults to gca)
% OUTPUT: 
%		Hout = handle to axis objects. if not called not returned.
% ---------------------------------------------------------------------------------------------
% Created :		07.01.2019.
% Modified:		07.01.2019.
% Copyleft:		Daniel Buncic.
% ---------------------------------------------------------------------------------------------}

% SetDefaultValue(1,'tickSize'			,3/4);
SetDefaultValue(1,'tickSize'			, [.7 .5]); % this is about symmetric size, NOT [.5 .5];
SetDefaultValue(2,'doubleTicks'	  , 1);
SetDefaultValue(3,'LabeloffSet'	  , 0);
SetDefaultValue(4,'currentHandle'	, gca);

ax(1) = currentHandle;
ax(2) = copyobj(ax(1), gcf);
% set(h1, 'TickDir', 'out', 'YTick', []);
set(ax(1), 'TickDir', 'out',  'XTick'     , []);
set(ax(2), 'TickDir', 'in'	, 'YTickLabel', []);

tickshrink(tickSize, ax(2));

if length(tickSize)>1
	tickshrink(tickSize(2), ax(1));
	tickshrink(tickSize(1), ax(2));
else
	tickshrink(tickSize, ax(2));
end

if doubleTicks
  linkaxes(ax)
  % ax(2).YAxis.Visible = 'off';
  % ax(2).XAxisLocation = 'bottom';
  % ax(2).XAxis.TickDirection = 'both';
  ax(1).XAxis.TickLabels = {};
  ax(2).XRuler.TickLabelGapOffset = 4+LabeloffSet;  % <-- adjust ticklabels away from axes
  ax(1).XAxis.TickLength=[.005 0];                 % <-- change the length of the tick marks
  ax(1).XTick = [ax(2).XTick];                      % <-- Xticks positions in the 2nd axes as in the 1st axes
  set(ax(1),'box','off');
else
  ax(2).XRuler.TickLabelGapOffset = LabeloffSet;  % <-- adjust ticklabels away from axes
end

if nargout
	Hout = [ax(1) ax(2)];
end













% %% old stuff
% ax(1) = gca;
% set(ax(1),'box','on');
% % define the second axes  and set some properties
% ax(2) = axes('Position',ax(1).Position,'Color','none');
% ax(2).XAxisLocation = 'bottom';
% ax(2).YAxis.Visible = 'off';
% ax(2).XAxis.TickDirection = 'out';
% ax(2).XTickLabels = [];
% ax(2).XTick = [0 ax(1).XTick];        % <-- Xticks positions in the 2nd axes as in the 1st axes
% % link the first axes with the second axes
% linkaxes(ax)
% % other properties
% ax(1).XRuler.TickLabelGapOffset = 5;  % <-- adjust ticklabels away from axes
% ax(2).XAxis.TickLength=[0.005 0];      % <-- change the length of the tick marks
% set(ax(2),'box','off');
% clf; close all;
% % my plot
% figure
% plot(1:9)
% % define the first axes (from the current axes) and set some properties
% ax(1) = gca;
% set(ax(1),'Fontsize',15,'FontWeight','Bold','LineWidth',2, 'box','on');
% % define the second axes  and set some properties
% ax(2) = axes('Position',ax(1).Position,'Color','none');
% ax(2).XAxisLocation = 'bottom';
% ax(2).YAxis.Visible = 'off';
% ax(2).XAxis.TickDirection = 'out';
% ax(2).XTickLabels = [];
% ax(2).XTick = [0 ax(1).XTick];        % <-- Xticks positions in the 2nd axes as in the 1st axes
% % link the first axes with the second axes
% linkaxes(ax)
% % other properties
% ax(1).XRuler.TickLabelGapOffset = 20; % <-- adjust ticklabels away from axes
% ax(2).XAxis.TickLength=[0.05 0];      % <-- change the length of the tick marks
% set(ax(2),'Fontsize',15,'FontWeight','Bold','LineWidth',2, 'box','off');
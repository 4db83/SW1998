function [] = add2yaxislabel(rightAlign)

if (nargin<1)||isempty(rightAlign); rightAlign = 0; end

% function create a second y-axis label.
ax1           = gca; % ax.YAxis.Exponent = 4;
fns1          = ax1.YAxis.FontSize;
ylim_					= get(gca,'YLim');
ylim_ticks 		= get(gca,'YTick');
ylim_lables_0 = get(gca,'YTickLabel');
ylim_lables 	= strrep(ylim_lables_0,' ','');
% if rightAlign second axis lables
if rightAlign == 1
  ylim_lables = num2str(ylim_ticks');
end
% make color a bit light than black as Matlab does
xclr = 0*ones(3,1);
set(ax1,'YColor',xclr);
Exp_ax1 = ax1.YRuler.Exponent;

% now set the set second axis here
ax2 = gca;
yyaxis(ax2,'right');
set(ax2,'YColor',xclr); 
set(ax2,'Ylim'  , ylim_);
set(ax2,'YTick' , ylim_ticks);
set(ax2,'FontSize', fns1);
if ~(Exp_ax1==0)
  set(ax2.YRuler, 'Exponent', Exp_ax1);
else
  set(ax2,'YTickLabel'  , ylim_lables);
end

% set the Exponent
% set(ax2,'YAxis.Exponent', Exp_ax1);
% OLD
% set(ax2,'YColor',xclr, ...
%         'YTickLabel', ylim_lables, ...
%         'YTick',ylim_ticks);






% 
% % set(ax2,'YTick'       , yticks_left);
% % set(ax2,'YTickLabel'  , ylim_lables);
% % ylim(ylim_left)
% % yticks(yticks_left)
% set(ax2,'YColor',xclr)
% % ax2.YTickLabel
% ax2.YRuler.Exponent = 3;
% % OLD
% % set(ax2,'YColor',xclr, ...
% %         'YTickLabel', ylim_lables, ...
% %         'YTick',ylim_ticks);


function varargout = setplot(figdim,y_digits,fntsize,axLineWidth,fighandle) 
%{F: specifies plots dimensions, fontsize and fontname
%===============================================================================
% Makes it easy to view strings in excel which are actually numbers without 
% converting them to numbers`
%-------------------------------------------------------------------------------
% 	USGAGE:	(1) setplot([height], Fontsize)
%           (2) setplot([top-positon height], Fontsize)
%						(3) setplot([top-positon width height], Fontsize) (rightpos at 0)
%						(4) setplot([leftpos rightpos height width], Fontsize) 
%           [leftpos rightpos] is optional and will only be called if one 
%           needs to center the plot. Default font size is 10, fontname is Palation.
%-------------------------------------------------------------------------------
% 	INPUT : 
%	  figdim		=  (1x2) matrix with [height width] dimension of plot
%   fntsize		=  scalar, fontsize. default fontsize is 10.
%		fighandle =	 sclar, figure handle. 
% 	OUTPUT:       
%	  zero arguments: adjusted plot.
% OLD:
% function varargout = setplot(figdim,fntsize,y_digits,axLineWidth,fighandle) 
%===============================================================================
% 	NOTES :   Always call it at the end of the plotting commands, after hold off;
%-------------------------------------------------------------------------------
% Created :		27.07.2013.
% Modified:		27.07.2013.
% Copyleft:		Daniel Buncic.
%------------------------------------------------------------------------------%}
widh_def = .9;
% SetDefaultValue(1, 'figdim'			, [widh_def .2]);
pos0 = get(gca,'Position');
SetDefaultValue(1, 'figdim'			, pos0);

FNS_0 = get(gca,'FontSize');
FNS0	= FNS_0 - 0;

xd = size(figdim);
if ~(max(xd) == 1)
	Big1 = figdim > 1;
	if any(Big1) > 0
		figdim = figdim./10;
	end
end

SetDefaultValue(2, 'y_digits'		, []);
SetDefaultValue(3, 'fntsize'		, FNS0);
SetDefaultValue(4, 'axLineWidth', 1); 
SetDefaultValue(5, 'fighandle'	, gca);

% added the font name conversion now to print in TNR as opposed to Palation!
FName = 'Times New Roman';
top_	= .45;
left_ = .05;

% if nargin < 2
  % if only one input argument is given
  if max(xd) == 2
    % sets width and height, takes the centered location [.1 .2]
%     set(fighandle,'Position',[[left_ top_] figdim],'FontSize',fntsize,'FontName',FName);
    set(fighandle,'Position',[left_ figdim(1) widh_def figdim(2)],'FontSize',fntsize,'FontName',FName);
  elseif max(xd) == 3
    % sets width and height as well as the location
    set(fighandle,'Position',[left_ figdim(1:3)],'FontSize',fntsize,'FontName',FName);
	elseif max(xd) == 4
    % sets width and height as well as the location
    set(fighandle,'Position',figdim,'FontSize',fntsize,'FontName',FName);
  elseif max(xd) == 1
    % only scalar input means that only the font size is to be set.
%     fntsize = figdim;
%     set(fighandle,'FontSize',fntsize,'FontName',FName);
		set(fighandle,'Position',[[left_ top_ widh_def] figdim],'FontSize',fntsize,'FontName',FName);
  else 
    disp('Error in plot diminsion');
	end

% deal with the tick format display  
yt0 = get(gca,'YTick')';
yts = str2double(get(gca,'YTickLabel'));
zz  = log10(yt0(end)./yts(end));
ax  = gca; 
abc = ax.YRuler.Exponent;

isI = min(floor(yts)==yts);

ytcks = get(gca,'YTick');
if ~isempty(y_digits)
  setyticklabels(ytcks)
  setytick(y_digits,fntsize,fighandle)
%   % set the ylim digits 
%   ytcks = get(gca,'YTick');
%   if ~isempty(y_digits)
% 	  setyticklabels(ytcks)
% 	  setytick(y_digits,fntsize,fighandle)
else
  if zz < -2
    ax.YRuler.TickLabelFormat = '%.0f';
    ax.YRuler.Exponent = abc - 1*~isI;
  else
    yt0 = get(gca,'YTick');
    yts = get(gca,'YTickLabel')';
    tmp_tcks = yts';
    jj = zeros(length(tmp_tcks),1);
    for i = 1:length(tmp_tcks)
      tmpstri = strfind(tmp_tcks(i),'.');
      if ~isempty(tmpstri{1})
        jj(i) = length(char(tmp_tcks(i)))- tmpstri{1};
      end
    end
    
    if max(jj) == 0 
      dgts = 0; 
    else
      dgts = max(jj) ;
    end

    ytck = str2double(yts)';
    % find 0 value and set to 0 without digits
    f0   = find((abs(ytck)<10e-14));       
    s22f = ['%2.' num2str(dgts) 'f'];  %
    % this is the new y-ticke lable with formatting s22f
    ffyy = cellstr(num2str(ytck,s22f));
    % add normal 0 lable for 0 axis
    if ~isempty(f0); ffyy{f0} = '0'; end
    set(gca,'YTickLabel',ffyy);
  end
end

% make linewidth of box
ax.LineWidth = axLineWidth;

% plotting standards, can be removed
% tickshrink(2/3);
% box on; grid on;
% set(gca,'GridLineStyle',':','GridAlpha',1/3)

plot_position = get(gca,'Position');

if nargout > 0
  varargout{1} = plot_position;
end


% OLD STUUFF
% end

%   y_digits = numel(num2str(abs(yt0(end))))
%   y_digits = max(numel(yts{end-1}),numel(yts{end}))-2;
%   if max(yt0)<0
%     y_digits = max(numel(yts{end-1}),numel(yts{end}))-3;
%   else
%     y_digits = max(numel(yts{end-1}),numel(yts{end}))-3;
%   end
%   if abs(yt2 - 1) < 10e-6
%     y_digits = 2;
%   else 
%     y_digits = log10(yt2)+1;
%   end
%     setytick(y_digits,fntsize,fighandle)
% setyticklabels(ytcks)
% end

% get axes handle and increase the thickness of lnes and choose grid value for color


% 	ax.GridAlpha = .15;    
	% else
  % if two input arguments are given then same as above except for last bit
% %   pst0 = get(fighandle,'Position');
% % 	pst0(1:2) = [.05 .55];
% %   if max(xd) == 2
% %     set(fighandle,'Position',[pst0(1:2) figdim],'FontSize',fntsize,'FontName',FName);
% %   elseif max(xd) == 3
% %     set(fighandle,'Position',[figdim(1) ff(2) figdim(2:3)],'FontSize',fntsize,'FontName',FName);
% % 	elseif max(xd) == 4
% %     set(fighandle,'Position',figdim,'FontSize',fntsize,'FontName',FName);    
% %   else 
% %     disp('Error in plot diminsion');
% %   end
% % end
% set axis laver over the top of the plots.
% set(gca, 'Layer','top')	
%set(gca,'Position',[[.1 .2] figdim],'FontSize',fntsize,'FontName',FName);


% % % if nargin < 2
% % %   % if only one input argument is given
% % %   if max(xd) == 2
% % %     figdim = figdim;
% % %     fntsize = 10;     % fontsize is 10 as default value
% % %     % sets width and height, takes the centered location [.1 .2]
% % %     set(fighandle,'Position',[[.055 .2] figdim],'FontSize',fntsize,'FontName',FName);
% % %   elseif max(xd) == 3
% % %     figdim = figdim;
% % %     fntsize = 10;     % fontsize is 10 as default value
% % %     % sets width and height as well as the location
% % %     set(fighandle,'Position',[figdim(1) 0 figdim(2:3)],'FontSize',fntsize,'FontName',FName);
% % % 	elseif max(xd) == 4
% % %     figdim = figdim;
% % %     fntsize = 10;     % fontsize is 10 as default value
% % %     % sets width and height as well as the location
% % %     set(fighandle,'Position',figdim,'FontSize',fntsize,'FontName',FName);
% % %   elseif max(xd) == 1
% % %     % only scalar input means that only the font size is to be set.
% % %     fntsize = figdim;
% % %     set(fighandle,'FontSize',fntsize,'FontName',FName);
% % %   else 
% % %     disp('Error in plot diminsion');
% % %   end
% % % else 
% % %   % if two input arguments are given then same as above except for last bit
% % %   pst0 = get(fighandle,'Position');
% % % 	pst0(1:2) = [.05 .55];
% % %   if max(xd) == 2
% % %     set(fighandle,'Position',[pst0(1:2) figdim],'FontSize',fntsize,'FontName',FName);
% % %   elseif max(xd) == 3
% % %     set(fighandle,'Position',[figdim(1) ff(2) figdim(2:3)],'FontSize',fntsize,'FontName',FName);
% % % 	elseif max(xd) == 4
% % %     set(fighandle,'Position',figdim,'FontSize',fntsize,'FontName',FName);    
% % %   else 
% % %     disp('Error in plot diminsion');
% % %   end
% % % end

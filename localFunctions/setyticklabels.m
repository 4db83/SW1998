function varargout= setyticklabels(ticks,digitsN,FNTsize,axLineWidth,fig_handle)
% set or adjust the yticklables and ticks at desired location vector ticks.
% =======================================================================================
% input
% ticks:			(kx1), location vector of ticks which will also be the lables.
% fig_handle:	scalar, figure handle, default is gca.
%
% call as:		setyticklabel([-.1:.1.1]) to set ticks and ticklables as the vector.
%							[-.1:.1.1] in the y-axis ticks.
% =======================================================================================
% db 14.06.2012
% modfified on 10.07.2013
% =======================================================================================

FNS_0 = get(gca,'FontSize');
set(gca,'FontName','Times New Roman');
FNS0	= FNS_0 - 0;

% SetDefaultValue(2, 'digitsN'	  , Ndec);
SetDefaultValue(3, 'FNTsize'   , FNS0);
axLineWidth0 = get(gca,'LineWidth');
SetDefaultValue(4, 'axLineWidth', axLineWidth0);
SetDefaultValue(5, 'fig_handle' , gca);

box_clr = 0.0*ones(1,3);
set(gca,'YColor',box_clr);
set(gca,'XColor',box_clr);
set(gca,'LineWidth', axLineWidth);

% get default ytick digits and Lables
ax = gca;
Ytcks0      = get(gca,'YTick');
ytckLabls0  = get(gca,'YTickLabel');

% DEFAULT TICKS IF NO INPUT
if nargin == 0 || isempty(ticks)
  ticks = Ytcks0; 
else 
  set(gca,'YTick'			, ticks)
  set(gca,'YTickLabel', cellstr(num2str(ticks')))
  Ytcks0     = get(gca, 'YTick');
  ytckLabls0 = get(gca, 'YTickLabel');
end

% make sure the two lenghts are the same
if length(ytckLabls0) == length(ytckLabls0)
  % find positive values 
  I = Ytcks0 > 0;
  Isum = sum(I);
  if Isum > 0
    IytckLabls = ytckLabls0(I);
    nc = size(char( IytckLabls ),2);
    for ii=1:length(IytckLabls)
      FP = find( char(IytckLabls(ii)) == '.', 1);
      if isempty(FP); FP = 0; end
      Iplus(ii) = FP;
    end
    if max(Iplus) == 0
      Ndec = 0;
    else
      Ndec = nc - max(Iplus);
    end
  else % Isum = 0
    IytckLabls = ytckLabls0(~I);
    nc = size(char( IytckLabls ),2);
    for ii=1:length(IytckLabls)
      FP = find( char(IytckLabls(ii)) == '.', 1);
      if isempty(FP); FP = 0; end
      Iminus(ii) = FP;
    end
    if max(Iminus) == 0
      Ndec = 0; 
    else
      Ndec = nc - max(Iminus);
    end
  end
else
  Ndec = 2;
end

% set default digits
SetDefaultValue(2, 'digitsN'	  , Ndec);

if ticks(1)<ticks(end)
	ylim(fig_handle,[ticks(1) ticks(end)]);
end

% digitsN
if (digitsN == -1)
  % set(fig_handle,'YTick'			,Ytcks0);
  % set(fig_handle,'YTickLabel'	,ytckLabls0);
else
  % FIND THE ZERO CROSSING
  % f0		= find(ticks==0);           % find 0 value.
  f0		= find(abs(ticks)<5*eps);
  s22f	= ['%2.' num2str(digitsN) 'f'];  %
  
  % this is the new y-ticke lable with formatting s22f
  ffyy = cellstr(num2str(ticks',s22f));
  
  % add normal 0 lable for 0 axis
  if ~isempty(f0)
  	ffyy{f0} = '0';
  end
  % ffyy
  % ticks
  % ytckLabls0
  set(fig_handle,'YTick'			, ticks);
  set(fig_handle,'YTickLabel'	, ffyy);
  set(fig_handle,'FontSize'		, FNTsize);
end

% set fontsize here
ax = gca;  % Get current axes
% Change Y-axis tick label font size
ax.YAxis.FontSize = FNTsize;

if nargout > 0
	varargout{1} = fig_handle;
end


end


% ytckLabls = get(gca,'Yticklabels')
% yend = char(ytckLabls(end))
% Ndot = find(yend == '.', 1)
% dot2end = yend(Ndot:end)
% Ndec = length(yend) - length(dot2end)
% % count_decimals0(x)
% 
% num = x;
% 
% % function decimalCount = count_decimals0(num)
%     % Convert the number to a string
%     numStr = num2str(num, '%.16f'); % Use a high precision to avoid truncation issues
%     % Find the position of the decimal point
%     decimalPointPos = find(numStr == '.', 1)
%     % If there is no decimal point, return 0
%     if isempty(decimalPointPos)
%         decimalCount = 0;
%     else
%         % Count the number of characters after the decimal point
%         decimalCount = length(numStr) - decimalPointPos;
%         % Remove trailing zeros
%         decimalCount = length(strrep(numStr(decimalPointPos+1:end), '0', ''));
%     end
% % end
% if Isum > 0
  % find all positive Ytcks
% ychar = char(ytckLabls);
% [nr,nc] = size(ychar)
% yend = char(ytckLabls(end));
% 
% for ii = 1:length(ytckLabls)
%   if ~isempty( find(char(ytckLabls(ii)) == '.', 1) )  
%     Ndot(ii) = find(char(ytckLabls(ii)) == '.', 1); 
%   end
%   % Ndot(ii) = if ~isempty(find(char(ytckLabls(ii)) == '.', 1)); end
% end
% % dot2end = ychar(Ndot:end);
% if isempty(Ndot)
%   Ndec = 0;
% else
%   Ndec = nc - Ndot;
% end

% Ndec
% SetDefaultValue(1, 'ticks'	    , Ytcks);
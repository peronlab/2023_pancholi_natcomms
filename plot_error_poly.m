%
% 2010 March S Peron
% 
% Function to make a line plot with error envelope
%
%   -time, vals, and err should be same size
%   -color is line color ; pcolor is color of surrounding envelope
%
% USAGE:
%  plot_error_poly (tvals, vals, err, color, pcolor)
%

function plot_error_poly (tvals, vals, err, color, pcolor)
  valids = find(~isnan(vals) & ~isnan(err));

  % First plot the SDs
  [x_err_poly, y_err_poly] = get_sem_poly(tvals(valids), vals(valids), err(valids));
  p=patch(x_err_poly,y_err_poly, pcolor, 'EdgeColor', 'none');
  alpha(p,0.5)
  hold on;
  plot(tvals,vals,'Color',color,'LineWidth',3);

  
%
% Returns your data as a polygon for bounding y at each x value with +/- y_off
%  for real nice SEM/SD plots
%
function [ret_x, ret_y] = get_sem_poly(x, y, y_off);
  ret_x = zeros(2*length(x),1);
  ret_y = zeros(2*length(x),1);
  
  l = length(x);
  
  for i=1:length(ret_x)
    if (i < l)
      ret_x(i) = x(i);
      ret_y(i) = y(i) + y_off(i);
    elseif (i == l || i == l+1)
      ret_x(i) = x(l);
      ret_y(i) = y(l) + y_off(l);
      if ( i == l + 1) ; ret_y(i) = y(l) - y_off(l); end
    else
      ret_x(i) = x(2*l - i + 1);
      ret_y(i) = y(2*l - i + 1) - y_off(2*l - i + 1);
    end
  end

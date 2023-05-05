%
% Corr that removes nan's and checks that you are using the right vector orientation
%
%  USAGE:
%
%    [r p] = nancorr(x,y,type)
%
%  ARGUMENTS:
%
%    x, y: vectors that must be of same length.  If entries are nan, they are 
%          ignored for correlation computation.  In the special case where
%          x is a matrix and y is not provided, the correlation of the columns
%          will be returned as a matrix.
%    type: Pearson, Spearmann etc.
%
function [r p] = nancorr(x,y,type)
	if (nargin < 3 || length(type) == 0) ; type = 'Pearson' ;end

    if (nargin > 1)
		if (size(x,2) ~=1 & size(x,1) == 1) ; x = x'; end
		if (size(y,2) ~=1 & size(y,1) == 1) ; y = y'; end
		vali = find(~isnan(x) & ~isnan(y)); 

		if (length(vali) > 0)
			[r p] = corr (x(vali),y(vali),'type',type);
		else
			r = nan;
            p = nan;
		end
    else 
	    r = ones(size(x,2));
        p = nan*ones(size(x,2));
        if (nargout == 1)
            for i1=1:size(x,2)
                tic;
                vali1 = find(~isnan(x(:,i1)));
                for i2=i1+1:size(x,2)
                    vali = intersect(vali1, find(~isnan(x(:,i2))));
                    if (length(vali) > 0)
                        rval = corr(x(vali,i1),x(vali,i2),'type',type);
                        r(i1,i2) = rval;
                        r(i2,i1) = rval;
                    end
                end
                t_elapsed = toc;
                disp(['Done with ' num2str(i1) ' out of ' num2str(size(x,2)) ' loop time: ' num2str(t_elapsed)]); 
            end
        else
            for i1=1:size(x,2)
                vali1 = find(~isnan(x(:,i1)));
                for i2=i1+1:size(x,2)
                    vali = intersect(vali1, find(~isnan(x(:,i2))));
                    if (length(vali) > 0)
                        [rval pval] = corr(x(vali,i1),x(vali,i2),'type',type);
                        r(i1,i2) = rval;
                        r(i2,i1) = rval;
                        p(i1,i2) = pval;
                        p(i2,i1) = pval;
                    end
                end
            end
        end
	end


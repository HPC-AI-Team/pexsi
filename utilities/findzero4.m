function [xHist, fHist, gHist] = findzero( fFunc, gFunc, x0, ...
	xmin, xmax, fmin, fmax, gmin, gmax, maxiter, fTol )
% FINDZERO4 uses the Newton's method with quadratic interpolation as
% safeguard.
%
% Lin Lin
% 12/28/2012

if( fmin >= 0 || fmax <= 0 )
	error('Problem with the initial interval.');
end

x     = x0;
xHist = [];
fHist = [];
gHist = [];

for iter = 1 : maxiter
	fVal = fFunc(x);
	gVal = gFunc(x);

	xHist = [xHist; x];
	fHist = [fHist; fVal];
	gHist = [gHist; gVal];

	if( abs( fVal ) < fTol )
		fprintf('Found zero in %d steps, xsol=%25.15f.\n', iter, x );
		break;
	end	

	% Compute the newton update
	xNewton = x - fVal / gVal;

	% Shrink interval
	if( fVal < 0 )
		xmin = x;
		fmin = fVal;
		gmin = gVal;
	else
		xmax = x;
		fmax = fVal;
		gmax = gVal;
	end

	if( xNewton < xmin || xNewton > xmax )
		% Quadratic interpolation
		coef = [0; gVal; fVal];
		if( fVal < 0 )
			d = xmax - xmin;
			coef(1) = ( fmax - fVal - gVal * d ) / d^2;
			x = roots(coef);
			x = x (find(x>0 & x < d)) + xmin;
		else
			d = xmin - xmax;
			coef(1) = ( fmin - fVal - gVal * d ) / d^2;
			x = roots(coef);
			x = x ( find(x > d & x < 0) ) + xmax;
		end
	else
		x = xNewton;
	end

end

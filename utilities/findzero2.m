function [xHist, fHist, gHist] = findzero( fFunc, gFunc, x0, ...
	xmin, xmax, fmin, fmax, gmin, gmax, maxiter, fTol )
% FINDZERO2 uses the Newton's method with bisection as safeguard.
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

	xold = x;

	% Compute the newton update
	xNewton = x - fVal / gVal;

	if( xNewton < xmin || xNewton > xmax )
		% Bisection
		if( fVal < 0 )
			x = (x + xmax ) / 2;
		else
			x = (x + xmin ) / 2;
		end
	else
		x = xNewton;
	end

	% Shrink interval
	if( fVal < 0 )
		xmin = xold;
		fmin = fVal;
		gmin = gVal;
	else
		xmax = xold;
		fmax = fVal;
		gmax = gVal;
	end

end

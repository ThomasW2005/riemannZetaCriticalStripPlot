#include "zeta.h"

const long double LOWER_THRESHOLD = 1.0e-6;
const long double UPPER_BOUND = 1.0e+4;
const int MAXNUM = 100;

std::complex<long double> cdgamma(std::complex<long double> x)
{
	std::complex<long double> y;
	long double xr, xi, wr, wi, ur, ui, vr, vi, yr, yi, t;

	xr = x.real();
	xi = x.imag();
	if (xr < 0) {
		wr = 1 - xr;
		wi = -xi;
	}
	else {
		wr = xr;
		wi = xi;
	}
	ur = wr + 6.00009857740312429;
	vr = ur * (wr + 4.99999857982434025) - wi * wi;
	vi = wi * (wr + 4.99999857982434025) + ur * wi;
	yr = ur * 13.2280130755055088 + vr * 66.2756400966213521 +
		0.293729529320536228;
	yi = wi * 13.2280130755055088 + vi * 66.2756400966213521;
	ur = vr * (wr + 4.00000003016801681) - vi * wi;
	ui = vi * (wr + 4.00000003016801681) + vr * wi;
	vr = ur * (wr + 2.99999999944915534) - ui * wi;
	vi = ui * (wr + 2.99999999944915534) + ur * wi;
	yr += ur * 91.1395751189899762 + vr * 47.3821439163096063;
	yi += ui * 91.1395751189899762 + vi * 47.3821439163096063;
	ur = vr * (wr + 2.00000000000603851) - vi * wi;
	ui = vi * (wr + 2.00000000000603851) + vr * wi;
	vr = ur * (wr + 0.999999999999975753) - ui * wi;
	vi = ui * (wr + 0.999999999999975753) + ur * wi;
	yr += ur * 10.5400280458730808 + vr;
	yi += ui * 10.5400280458730808 + vi;
	ur = vr * wr - vi * wi;
	ui = vi * wr + vr * wi;
	t = ur * ur + ui * ui;
	vr = yr * ur + yi * ui + t * 0.0327673720261526849;
	vi = yi * ur - yr * ui;
	yr = wr + 7.31790632447016203;
	ur = std::log(yr * yr + wi * wi) * 0.5 - 1;
	ui = std::atan2(wi, yr);
	yr = std::exp(ur * (wr - 0.5) - ui * wi - 3.48064577727581257) / t;
	yi = ui * (wr - 0.5) + ur * wi;
	ur = yr * std::cos(yi);
	ui = yr * std::sin(yi);
	yr = ur * vr - ui * vi;
	yi = ui * vr + ur * vi;
	if (xr < 0) {
		wr = xr * 3.14159265358979324;
		wi = std::exp(xi * 3.14159265358979324);
		vi = 1 / wi;
		ur = (vi + wi) * std::sin(wr);
		ui = (vi - wi) * std::cos(wr);
		vr = ur * yr + ui * yi;
		vi = ui * yr - ur * yi;
		ur = 6.2831853071795862 / (vr * vr + vi * vi);
		yr = ur * vr;
		yi = ur * vi;
	}

	return std::complex<long double>(yr, yi);
}

//analytic continuation of zeta function
std::complex<long double> zeta(const std::complex<long double>& s)
{
	std::complex<long double> a_arr[MAXNUM + 1];

	//constants
	std::complex<long double> half(0.5, 0.0);
	std::complex<long double> one(1.0, 0.0);
	std::complex<long double> two(2.0, 0.0);
	std::complex<long double> rev(-1.0, 0.0);

	std::complex<long double> prev(1.0e+20, 0.0);

	//initialize with a_0 = 0.5 / (1 - 2^(1-s))
	a_arr[0] = half / (one - std::pow(two, (one - s)));

	std::complex<long double> sum(0.0, 0.0);
	sum = sum + a_arr[0];


	for (int n = 1; n <= MAXNUM; n++)
	{
		std::complex<long double> n_c(n, 0.0);

		for (int k = 0; k < n; k++)
		{
			std::complex<long double> k_c(k, 0.0);

			a_arr[k] = half * a_arr[k] * (n_c / (n_c - k_c));
			sum = sum + a_arr[k];
		}

		a_arr[n] = (rev * a_arr[n - 1] * std::pow((n_c / (n_c + one)), s) / n_c);
		sum = sum + a_arr[n];

		//If the difference is less than or equal to the threshold value, it is considered to be convergent and the calculation is terminated.
		//if( (std::abs(prev - sum) / std::abs( sum )) < LOWER_THRESHOLD )
		if (std::abs(prev - sum) < LOWER_THRESHOLD)
			break;
		
		//doesn't work for large values, so it gets terminated when it exceeds UPPER_BOUND
		if (std::abs(sum) > UPPER_BOUND)
			break;

		prev = sum;
	}

	return sum;
}

//Complex zeta function streamlined using functional equations
std::complex<long double> complex_zeta(std::complex<long double> s)
{
	long double real = s.real();
	long double imag = s.imag();

	if (real < 0.0)
	{
		// Flipping left and right in a function equation
		std::complex<long double> inv_s(1 - real, -imag);
		std::complex<long double> inv_one(-1, 0);
		std::complex<long double> two(2.0, 0.0);
		std::complex<long double> pi_c(3.141592653589793238462643383279, 0.0);

		return std::pow(two, s) * std::pow(pi_c, inv_one * inv_s) * sin(pi_c * s / two) * cdgamma(inv_s) * zeta(inv_s);
	}
	else
		return zeta(s);

}
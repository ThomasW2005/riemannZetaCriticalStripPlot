#include "zeta.h"

const long double LOWER_THRESHOLD = 1.0e-6;
const long double UPPER_BOUND = 1.0e+4;
const int MAXNUM = 100;

std::complex<long double> Γ(std::complex<long double> x)
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

//zeta function
std::complex<long double> zeta(const std::complex<long double>& s)
{
	std::complex<long double> a_arr[MAXNUM + 1];
	std::complex<long double> half(0.5, 0.0);
	std::complex<long double> one(1.0, 0.0);
	std::complex<long double> two(2.0, 0.0);
	std::complex<long double> rev(-1.0, 0.0);
	std::complex<long double> sum(0.0, 0.0);
	std::complex<long double> prev(1.0e+20, 0.0);
	
	a_arr[0] = half / (one - std::pow(two, (one - s))); //initialize with a_0 = 0.5 / (1 - 2^(1-s))
	sum += a_arr[0];

	for (int n = 1; n <= MAXNUM; n++)
	{
		std::complex<long double> nCplx(n, 0.0); //complex index

		for (int k = 0; k < n; k++)
		{
			std::complex<long double> kCplx(k, 0.0); //complex index

			a_arr[k] *= half * (nCplx / (nCplx - kCplx));
			sum += a_arr[k];
		}

		a_arr[n] = (rev * a_arr[n - 1] * std::pow((nCplx / (nCplx + one)), s) / nCplx);
		sum += a_arr[n];


		if (std::abs(prev - sum) < LOWER_THRESHOLD)//If the difference is less than or equal to the threshold value, it is considered to be convergent and the calculation is terminated.
			break;

		if (std::abs(sum) > UPPER_BOUND)//doesn't work for large values, so it gets terminated when it exceeds UPPER_BOUND
			break;

		prev = sum;
	}
	return sum;
}

//analytic continuation of zeta function
std::complex<long double> ζ(std::complex<long double> s)
{
	if (s.real() < 1.0)
	{
		std::complex<long double> oneCplx(1.0, 0.0);
		std::complex<long double> twoCplx(2.0, 0.0);
		std::complex<long double> piCplx(std::numbers::pi, 0.0);

		return std::pow(twoCplx, s) * std::pow(piCplx, s - oneCplx) * sin(piCplx * s / twoCplx) * Γ(oneCplx - s) * zeta(oneCplx - s);
	}
	else
		return zeta(s);
}
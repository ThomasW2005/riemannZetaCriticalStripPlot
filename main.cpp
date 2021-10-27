/*
Thomas Weichhart
27.10.2021 v1.3
std::riemann_zeta ζ
Refs:
	https://github.com/junpeitsuji/iwannatouchzeta
	http://mathonline.wikidot.com/the-complex-cosine-and-sine-functions
#include <Windows.h> + SetConsoleOutputCP(CP_UTF8);
*/

#include <SDL.h>
#include <iostream>
#include <complex>
#include <deque>
#include <algorithm>
#include "zeta.h"

constexpr int WIDTH = 1280;
constexpr int HEIGHT = 720;
constexpr float AR = (float)HEIGHT / WIDTH;

long double map(long double x, long double in_min, long double in_max, long double out_min, long double out_max);
void translateComplex(std::complex<long double> n, int& x, int& y, long double xMin, long double xMax, long double iMin, long double iMax);
SDL_Color HSVtoRGB(float h, float s = 1., float v = 1.);

int main(int argc, char* argv[])
{
	SDL_Init(SDL_INIT_EVERYTHING);
	SDL_Renderer* renderer;
	SDL_Window* window;
	SDL_Event event;
	bool run = true;
	long double counter = 0;
	int x, y, lastX, lastY;
	std::deque<std::complex<long double>> points;
	SDL_CreateWindowAndRenderer(WIDTH, HEIGHT, 0, &window, &renderer);
	SDL_SetWindowTitle(window, "Riemann ζ(Zeta) Function Critical Strip Plot ");
	SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
	SDL_RenderClear(renderer);

	while (run)
	{
		SDL_SetRenderDrawColor(renderer, 50, 50, 50, 255);
		SDL_RenderDrawLine(renderer, 0, HEIGHT / 2, WIDTH, HEIGHT / 2);
		SDL_RenderDrawLine(renderer, WIDTH / 2, 0, WIDTH / 2, HEIGHT);

		SDL_Color c = HSVtoRGB(counter * 10);
		SDL_SetRenderDrawColor(renderer, c.r, c.g, c.b, 255);
		auto s = complex_zeta(std::complex<long double>(0.5, counter));

		points.push_back(s);
		if (points.size() > 1000)
			points.pop_front();

		for (int i = 0; i < points.size() - 1; i++)
		{
			long double adjustedCounter = counter - .010 * i;
			SDL_Color c = HSVtoRGB(adjustedCounter * 10, 1., std::clamp((float)map(adjustedCounter, counter - .010 * points.size(), counter, 3.0, 0.0), (float)0.0, (float)1.0));
			SDL_SetRenderDrawColor(renderer, c.r, c.g, c.b, 255);

			translateComplex(points[i], x, y, -4, 4, -4, 4);
			translateComplex(points[i + 1], lastX, lastY, -4, 4, -4, 4);
			SDL_RenderDrawLine(renderer, lastX, lastY, x, y);
		}
		SDL_RenderPresent(renderer);
		SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
		SDL_RenderClear(renderer);
		counter += .01;

		Uint32 timeout = SDL_GetTicks() + 1;
		do
		{
			while (SDL_PollEvent(&event)) //eventhandling
			{
				switch (event.type)
				{
				case SDL_MOUSEMOTION:
				case SDL_MOUSEBUTTONUP:
				case SDL_MOUSEBUTTONDOWN:
					break;
				case SDL_QUIT:
					run = false;
					break;
				case SDL_KEYDOWN:
					switch (event.key.keysym.sym)
					{
					case SDLK_ESCAPE:
						run = false;
						break;
					}
					break;
				}
			}
		} while (!SDL_TICKS_PASSED(SDL_GetTicks(), timeout));
	}
	SDL_DestroyRenderer(renderer);
	SDL_DestroyWindow(window);
	SDL_Quit();
	return 0;
}

long double map(long double x, long double in_min, long double in_max, long double out_min, long double out_max)
{
	return (x - in_min) * (out_max - out_min) / (in_max - in_min) + out_min;
}

void translateComplex(std::complex<long double> n, int& x, int& y, long double xMin, long double xMax, long double iMin, long double iMax)
{
	long double real = n.real();
	long double imag = n.imag();

	int offset = WIDTH / 2 - map(0, xMin, xMax, 0, WIDTH * AR);
	x = map(real, xMin, xMax, 0, WIDTH * AR) + offset;
	y = map(imag, iMin, iMax, HEIGHT, 0);

	//SDL_RenderDrawPoint(renderer, x, y);
	//std::cout << real << ":" << imag << "xy:";
	//std::cout << x << ":" << y << '\n';
}

SDL_Color HSVtoRGB(float h, float s, float v)
{
	if (h >= 359.9999f)
	{
		int a = h / 360.;
		h -= 360 * a;
	}
	float r, g, b;
	int i;
	float f, p, q, t;
	if (s == 0)
		return{ 0, 0, 0 };
	h /= 60;
	i = floor(h);
	f = h - i;
	p = v * (1 - s);
	q = v * (1 - s * f);
	t = v * (1 - s * (1 - f));
	switch (i) {
	case 0: r = v; g = t; b = p; break;
	case 1: r = q; g = v; b = p; break;
	case 2: r = p; g = v; b = t; break;
	case 3: r = p; g = q; b = v; break;
	case 4: r = t; g = p; b = v; break;
	default: r = v; g = p; b = q; break;
	}
	return { (unsigned char)(255 * r),(unsigned char)(255 * g),(unsigned char)(255 * b) };
}
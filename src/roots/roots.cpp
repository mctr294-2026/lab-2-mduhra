#include "roots.hpp"
#include <cmath>
#include <iostream>

const double TOLERANCE = 1e-6;
const int MAX_ITERATIONS = 1e6;

// Bisection Method
bool bisection(std::function<double(double)> f, double a, double b, double *root) {
    double fa = f(a);
    double fb = f(b);

    // Requires opposite signs at the boundaries poiints. There is a chance there isn't a root.
    if (fa * fb > 0) return false;

    for (int i = 0; i < MAX_ITERATIONS; ++i) {
        double mid = (a + b) / 2.0;  // Calculates the mid point
        double fmid = f(mid);

        if (std::abs(fmid) < TOLERANCE || (b - a) / 2.0 < TOLERANCE) {  // if the middle point is close enough to zero, the interval shrinks
            *root = mid;
            return true;
        }

        if (fa * fmid < 0) {  // if the root is on the left side, shift the right boundary to the middle, vice versa.
            b = mid;
            fb = fmid;
        } else {
            a = mid;
            fa = fmid;
        }
    }
    return false;
}

// Regula Falsi
bool regula_falsi(std::function<double(double)> f, double a, double b, double *root) {   // straight line between midpoints
    double fa = f(a);
    double fb = f(b);

    if (fa * fb > 0) return false;

    for (int i = 0; i < MAX_ITERATIONS; ++i) {

        double c = b - (fb * (b - a)) / (fb - fa);   // checks if C is close to the root
        double fc = f(c);

        if (std::abs(fc) < TOLERANCE) {
            *root = c;
            return true;
        }

        if (fa * fc < 0) {   // shrinks the interval based on where the sing change occurs
            b = c;
            fb = fc;
        } else {
            a = c;
            fa = fc;
        }
    }
    return false;
}

// Newton-Raphson
bool newton_raphson(std::function<double(double)> f, std::function<double(double)> g,  // calls secnd function to get the slope at x
                    double a, double b, double c, double *root) {
    double x = c;

    for (int i = 0; i < MAX_ITERATIONS; ++i) {
        double fx = f(x);
        double gx = g(x); 

        if (std::abs(gx) < 1e-12) return false; // Can't divide by 0, therefore test would fail

        double x_next = x - fx / gx;   // jumps to next guess using tangent line

        if (x_next < a || x_next > b) return false; // jump stay within the interval [a, b]

        if (std::abs(x_next - x) < TOLERANCE || std::abs(f(x_next)) < TOLERANCE) {
            *root = x_next;
            return true;
        }
        x = x_next;
    }
    return false;
}

// Secant Method
bool secant(std::function<double(double)> f, double a, double b, double c, double *root) {
    double x1 = c;

    double x0 = (x1 + TOLERANCE <= b) ? x1 + TOLERANCE : x1 - TOLERANCE;    // create another point near C as a secant line needs two points

    for (int i = 0; i < MAX_ITERATIONS; ++i) {
        double fx0 = f(x0);
        double fx1 = f(x1);

        if (std::abs(fx1 - fx0) < 1e-12) return false; // Cant divide by 0. This is because two functions had the same y value

        double x_next = x1 - fx1 * (x1 - x0) / (fx1 - fx0);

        if (x_next < a || x_next > b) return false;  // slope must stay within the range

        if (std::abs(x_next - x1) < TOLERANCE || std::abs(f(x_next)) < TOLERANCE) {
            *root = x_next;
            return true;
        }

        x0 = x1;
        x1 = x_next;
    }
    return false;
}
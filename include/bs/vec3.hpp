#pragma once

#include <cmath>
#include <algorithm>
#include <array>
#include <complex>
#include <limits>
#include <string>

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884
#endif

// ================================================================
// Vec3
// ================================================================
struct Vec3
{
    double x, y, z;
};

inline Vec3 operator+(const Vec3& a, const Vec3& b)
{
    return { a.x + b.x, a.y + b.y, a.z + b.z };
}

inline Vec3 operator-(const Vec3& a, const Vec3& b)
{
    return { a.x - b.x, a.y - b.y, a.z - b.z };
}

inline Vec3 operator*(double s, const Vec3& v)
{
    return { s * v.x, s * v.y, s * v.z };
}

inline Vec3& operator+=(Vec3& a, const Vec3& b)
{
    a.x += b.x;
    a.y += b.y;
    a.z += b.z;
    return a;
}

inline double dot(const Vec3& a, const Vec3& b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

inline double norm(const Vec3& v)
{
    return std::sqrt(dot(v, v));
}

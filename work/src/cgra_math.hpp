//---------------------------------------------------------------------------
//
// Copyright (c) 2016 Taehyun Rhee, Joshua Scott, Ben Allen
//
// This software is provided 'as-is' for assignment of COMP308 in ECS,
// Victoria University of Wellington, without any express or implied warranty. 
// In no event will the authors be held liable for any damages arising from
// the use of this software.
//
// The contents of this file may not be copied or duplicated in any form
// without the prior permission of its owner.
//
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
// 
// CGRA Math Library
//
// Vector math library that mimics api and behaviour of GLSL without being as
// complicated as GLM. Functionality includes vector2, vector3, vector4, matrix2,
// matrix3 and matrix4 and their respective data types (float, double, int,
// unsigned and bool) and functions.
//
// Important things to note, vectors and matrices are mutable. Matrices
// are column major (stored with columns in consectutive order) and are 
// accessed by mat[column][row].
//
// Features NOT avaliable include swizzling and some rarely used functions:
// - Common Functions
// - - floatBitsToInt
// - - intBitsToFloat
// - - fma
// - - frexp
// - - ldexp
// - Floating-Point Pack and Unpack Functions
// - Vector Relational Functions
// - Integer Functions
//
//
// See the official GLSL docs for relevant documentation
// - https://cvs.khronos.org/svn/repos/ogl/trunk/ecosystem/public/sdk/docs/man4/html/
// - https://www.opengl.org/registry/doc/GLSLangSpec.4.50.pdf
//
// Or these other links that may be more helpful
// - http://www.shaderific.com/glsl-functions/   (recommended)
// - https://web.eecs.umich.edu/~sugih/courses/eecs487/common/notes/APITables.xml
// - http://docs.gl/
//
//
// NOTE: This is a new library which hasn't be throughly tested. If you suspect
// any bugs that are reproducible, please let the creators know, and they'll fix it.
//
//
//----------------------------------------------------------------------------

#pragma once

#include <cassert>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <type_traits>

namespace cgra {

	template <typename> class vector2;
	using  vec2 = vector2<float>;
	using dvec2 = vector2<double>;
	using ivec2 = vector2<int>;
	using uvec2 = vector2<unsigned>;
	using bvec2 = vector2<bool>;

	template <typename> class vector3;
	using  vec3 = vector3<float>;
	using dvec3 = vector3<double>;
	using ivec3 = vector3<int>;
	using uvec3 = vector3<unsigned>;
	using bvec3 = vector3<bool>;

	template <typename> class vector4;
	using  vec4 = vector4<float>;
	using dvec4 = vector4<double>;
	using ivec4 = vector4<int>;
	using uvec4 = vector4<unsigned>;
	using bvec4 = vector4<bool>;

	template <typename> class matrix2;
	using  mat2 = matrix2<float>;
	using dmat2 = matrix2<double>;

	template <typename> class matrix3;
	using  mat3 = matrix3<float>;
	using dmat3 = matrix3<double>;

	template <typename> class matrix4;
	using  mat4 = matrix4<float>;
	using dmat4 = matrix4<double>;

	namespace math {

		// random
		template <typename T> inline T random(T lower = 0, T upper = 1) {
			static std::default_random_engine re { std::random_device()() };
			std::uniform_real_distribution<double> dist(lower, upper);
			return T(dist(re));
		}

		// pi
		inline double pi() {
			return 3.1415926535897932384626433832795;
		}

		// natural log base
		inline double e() {
			return 2.7182818284590452353602874713527;
		}

		// golden ratio
		inline double phi() {
			return 1.61803398874989484820458683436563811;
		}
	}

	template <typename T> inline T radians(T val) {
		return val * math::pi() / 180.0;
	}

	template <typename T> inline T degrees(T val) {
		return val / math::pi() * 180.0;
	}

	template <typename T> inline T log2(const T &a) {
		return std::log(a) * 1.4426950408889634073599246810019;
	}

	template <typename T> inline T exp2(const T &a) {
		return std::pow(2, a);
	}

	template <typename T> inline std::enable_if_t<std::is_arithmetic<T>::value, T> atan(const T &y, const T &x) {
		return std::atan2(y, x);
	}


	template <typename T> inline int sign(T val) {
		return (T(0) < val) - (val < T(0));
	}

	template <typename T> inline T inf() {
		// use like: inf<float>()
		// only for floating point types
		return std::numeric_limits<T>::infinity();
	}

	template <typename T> inline bool isinf(T a) {
		return std::numeric_limits<T>::max() < std::abs(a);
	}

	template <typename T> inline T nan() {
		// use like: nan<float>()
		// only for floating point types
		return T(0) / T(0);
	}

	template <typename T> inline bool isnan(T a) {
		return a != a;
	}


	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////    ___  _____   __      ________ _____ _______ ____  _____                                                        ////
	////   |__ \|  __ \  \ \    / /  ____/ ____|__   __/ __ \|  __ \                                                     ////
	////      ) | |  | |  \ \  / /| |__ | |       | | | |  | | |__) |                                                  ////
	////     / /| |  | |   \ \/ / |  __|| |       | | | |  | |  _  /                                                 ////
	////    / /_| |__| |    \  /  | |___| |____   | | | |__| | | \ \                                                   ////
	////   |____|_____/      \/   |______\_____|  |_|  \____/|_|  \_\                                                    ////
	////                                                                                                                   ////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	template <typename T>
	class vector2 {
	public:
		union{ T x; T r;};
		union{ T y; T g;};

		// T constructors
		vector2() : x(0), y(0) {}
		explicit vector2(T v) : x(v), y(v) {}
		vector2(T _x, T _y) : x(_x), y(_y) {}

		template <typename U>
		vector2(const vector2<U> &other) : x(other.x), y(other.y) { }

		static vector2 random(T lower = 0, T upper = 1) { 
			return vector2(math::random<T>(lower, upper), math::random<T>(lower, upper));
		}

		static vector2 i() {return vector2(1, 0);}
		static vector2 j() {return vector2(0, 1);}

		static vector2 checknan(const vector2 &v) {
			T sum = v.x + v.y;
			assert(sum == sum);
			return v;
		}

		explicit operator T *() {
			return &(x);
		}

		T * dataPointer() {
			return &(x);
		}

		const T * dataPointer() const{
			return &(x);
		}

		const T & operator[](size_t i) const {
			assert(i < 2);
			return *(&x + i);
		}

		T & operator[](size_t i) {
			assert(i < 2);
			return *(&x + i);
		}

		// stream insertion
		inline friend std::ostream & operator<<(std::ostream &out, const vector2 &v) {
			return out << '(' << v.x << ", " << v.y << ')';
		}


		// Operator overload - assign
		// 

		// assign
		template <typename U>
		vector2 & operator=(const vector2<U> &other) {
			x = other.x;
			y = other.y;
			return *this;
		}

		// add-assign
		template <typename U>
		vector2 & operator+=(const vector2<U> &rhs) {
			x += rhs.x;
			y += rhs.y;
			return *this;
		}
		
		// add-assign
		vector2 & operator+=(T rhs) {
			x += rhs;
			y += rhs;
			return *this;
		}

		// subtract-assign
		template <typename U>
		vector2 & operator-=(const vector2<U> &rhs) {
			x -= rhs.x;
			y -= rhs.y;
			return *this;
		}

		// subtract-assign
		vector2 & operator-=(T rhs) {
			x -= rhs;
			y -= rhs;
			return *this;
		}

		// mulitply-assign
		template <typename U>
		vector2 & operator*=(const vector2<U> &rhs) {
			x *= rhs.x;
			y *= rhs.y;
			return *this;
		}

		// mulitply-assign
		vector2 & operator*=(T rhs) {
			x *= rhs;
			y *= rhs;
			return *this;
		}

		// divide-assign
		template <typename U>
		vector2 & operator/=(const vector2<U> &rhs) {
			x /= rhs.x;
			y /= rhs.y;
			vector2::checknan(*this);
			return *this;
		}

		// divide-assign
		vector2 & operator/=(T rhs) {
			x /= rhs;
			y /= rhs;
			vector2::checknan(*this);
			return *this;
		}
	};



	// Vector / Vector Operator Overloads
	//

	// equality
	template <typename T>
	inline bool operator==(const vector2<T> &lhs, const vector2<T> &rhs) {
		return lhs.x == rhs.x && lhs.y == rhs.y;
	}

	// inequality
	template <typename T>
	inline bool operator!=(const vector2<T> &lhs, const vector2<T> &rhs) {
		return !(lhs == rhs);
	}

	// negate
	template <typename T>
	inline vector2<T> operator-(const vector2<T> &rhs) {
		return vector2<T>(-rhs.x, -rhs.y);
	}

	// add
	template <typename T1, typename T2>
	inline auto operator+(const vector2<T1> &lhs, const vector2<T2> &rhs) {
		vector2<std::common_type_t<T1, T2>> v = lhs;
		return v += rhs;
	}

	// subtract
	template <typename T1, typename T2>
	inline auto operator-(const vector2<T1> &lhs, const vector2<T2> &rhs) {
		vector2<std::common_type_t<T1, T2>> v = lhs;
		return v -= rhs;
	}

	// multiply
	template <typename T1, typename T2>
	inline auto operator*(const vector2<T1> &lhs, const vector2<T2> &rhs) {
		vector2<std::common_type_t<T1, T2>> v = lhs;
		return v *= rhs;
	}

	// divide
	template <typename T1, typename T2>
	inline auto operator/(const vector2<T1> &lhs, const vector2<T2> &rhs) {
		vector2<std::common_type_t<T1, T2>> v = lhs;
		return v /= rhs;
	}


	// Vector / Scalar Operator Overloads
	//

	// add right
	template <typename T1, typename T2>
	inline auto operator+(const vector2<T1> &lhs, T2 rhs) {
		vector2<std::common_type_t<T1, T2>> v = lhs;
		return v += rhs;
	}

	// add left
	template <typename T1, typename T2>
	inline auto operator+(T1 lhs, const vector2<T2> &rhs) {
		vector2<std::common_type_t<T1, T2>> v = rhs;
		return v += lhs;
	}

	// subtract right
	template <typename T1, typename T2>
	inline auto operator-(const vector2<T1> &lhs, T2 rhs) {
		vector2<std::common_type_t<T1, T2>> v = lhs;
		return v -= rhs;
	}

	// subtract left
	template <typename T1, typename T2>
	inline auto operator-(T1 lhs, const vector2<T2> &rhs) {
		vector2<std::common_type_t<T1, T2>> v(lhs);
		return v -= rhs;
	}

	// multiply right
	template <typename T1, typename T2>
	inline auto operator*(const vector2<T1> &lhs, T2 rhs) {
		vector2<std::common_type_t<T1, T2>> v = lhs;
		return v *= rhs;
	}

	// multiply left
	template <typename T1, typename T2>
	inline auto operator*(T1 lhs, const vector2<T2> &rhs) {
		vector2<std::common_type_t<T1, T2>> v = rhs;
		return v *= lhs;
	}

	// divide right
	template <typename T1, typename T2>
	inline auto operator/(const vector2<T1> &lhs, T2 rhs) {
		vector2<std::common_type_t<T1, T2>> v = lhs;
		return v /= rhs;
	}

	// divide left
	template <typename T1, typename T2>
	inline auto operator/(T1 lhs, const vector2<T2> &rhs) {
		vector2<std::common_type_t<T1, T2>> v(lhs);
		return v /= rhs;
	}


	// Angle and Trigonometry Functions
	// 

	// degrees to radians
	template <typename T>
	inline vector2<T> radians(const vector2<T> &d) {
		return vector2<T>(radians(d.x), radians(d.y));
	}

	// radians to degrees
	template <typename T>
	inline vector2<T> degrees(const vector2<T> &r) {
		return vector2<T>(degrees(r.x), degrees(r.y));
	}

	// sine
	template <typename T>
	inline vector2<T> sin(const vector2<T> &v) {
		return vector2<T>(std::sin(v.x), std::sin(v.y));
	}

	// cosine
	template <typename T>
	inline vector2<T> cos(const vector2<T> &v) {
		return vector2<T>(std::cos(v.x), std::cos(v.y));
	}

	// tangent
	template <typename T>
	inline vector2<T> tan(const vector2<T> &v) {
		return vector2<T>(std::tan(v.x), std::tan(v.y));
	}

	// arc sine
	template <typename T>
	inline vector2<T> asin(const vector2<T> &v) {
		return vector2<T>(std::asin(v.x), std::asin(v.y));
	}

	// arc cosine
	template <typename T>
	inline vector2<T> acos(const vector2<T> &v) {
		return vector2<T>(std::acos(v.x), std::acos(v.y));
	}

	// arc tangent of y/x
	template <typename T1, typename T2>
	inline auto atan(const vector2<T1> &y, const vector2<T2> &x) {
		using common_t = std::common_type_t<T1, T2>;
		return vector2<common_t>(std::atan2<common_t>(y.x, x.x), std::atan2<common_t>(y.y, x.y));
	}

	// arc tangent
	template <typename T>
	inline vector2<T> atan(const vector2<T> &v) {
		return vector2<T>(std::atan(v.x), std::atan(v.y));
	}


	// Exponential Functions
	// 

	// v raised to the power of e
	template <typename T1, typename T2>
	inline auto pow(const vector2<T1> &v, const vector2<T2> &e) {
		return vector2<std::common_type_t<T1, T2>>(std::pow(v.x, e.x), std::pow(v.y, e.y));
	}

	// natural exponentiation of vector
	template <typename T>
	inline vector2<T> exp(const vector2<T> &v) {
		return vector2<T>(std::exp(v.x), std::exp(v.y));
	}

	// natural logarithm of vector
	template <typename T>
	inline vector2<T> log(const vector2<T> &v) {
		return vector2<T>(std::log(v.x), std::log(v.y));
	}

	// base 2 exponentiation of vector
	template <typename T>
	inline vector2<T> exp2(const vector2<T> &v) {
		return vector2<T>(exp2(v.x), exp2(v.y));
	}

	// base 2 logarithm of vector
	template <typename T>
	inline vector2<T> log2(const vector2<T> &v) {
		return vector2<T>(log2(v.x), log2(v.y));
	}

	// square root of vector
	template <typename T>
	inline vector2<T> sqrt(const vector2<T> &v) {
		return vector2<T>(std::sqrt(v.x), std::sqrt(v.y));
	}

	// inverse of the square root of vector
	template <typename T>
	inline vector2<T> inversesqrt(const vector2<T> &v) {
		return vector2<T>::checknan(vector2<T>(T(1)/std::sqrt(v.x), T(1)/std::sqrt(v.y)));
	}


	// Common Functions
	// 

	// absolute value of vector
	template <typename T>
	inline vector2<T> abs(const vector2<T> &v) {
		return vector2<T>(std::abs(v.x), std::abs(v.y));
	}

	// sign (-1, 0, 1) of vector
	template <typename T>
	inline vector2<T> sign(const vector2<T> &v) {
		return vector2<T>(sign(v.x), sign(v.y));
	}

	// floor of vector
	template <typename T>
	inline vector2<T> floor(const vector2<T> &v) {
		return vector2<T>(std::floor(v.x), std::floor(v.y));
	}

	// ceil of vector
	template <typename T>
	inline vector2<T> ceil(const vector2<T> &v) {
		return vector2<T>(std::ceil(v.x), std::ceil(v.y));
	}

	// fractional part of vector : v-floor(v)
	template <typename T>
	inline vector2<T> fract(const vector2<T> &v) {
		return v-floor(v);
	}

	// modulas of vector : v-m*floor(v/m)
	template <typename T1, typename T2>
	inline auto mod(const vector2<T1> &v, T2 m) {
		return v-m*floor(v/m);
	}

	// modulas of vector : v-m*floor(v/m)
	template <typename T1, typename T2>
	inline auto mod(const vector2<T1> &v, const vector2<T2> &m) {
		return v-m*floor(v/m);
	}

	// minimum of vector components and T
	template <typename T1, typename T2, typename=std::enable_if_t<std::is_arithmetic<T2>::value>>
	inline auto min(const vector2<T1> &lhs, T2 rhs) {
		using common_t = std::common_type_t<T1, T2>;
		return vector2<std::common_type_t<T1, T2>>(std::min<common_t>(lhs.x, rhs), std::min<common_t>(lhs.y, rhs));
	}

	// minimum of vector components
	template <typename T1, typename T2>
	inline auto min(const vector2<T1> &lhs, const vector2<T2> &rhs) {
		using common_t = std::common_type_t<T1, T2>;
		return vector2<common_t>(std::min<common_t>(lhs.x, rhs.x), std::min<common_t>(lhs.y, rhs.y));
	}

	// maximum of vector components and T
	template <typename T1, typename T2, typename=std::enable_if_t<std::is_arithmetic<T2>::value>>
	inline auto max(const vector2<T1> &lhs, T2 rhs) {
		using common_t = std::common_type_t<T1, T2>;
		return vector2<common_t>(std::max<common_t>(lhs.x, rhs), std::max<common_t>(lhs.y, rhs));
	}

	// maximum of vector components
	template <typename T1, typename T2>
	inline auto max(const vector2<T1> &lhs, const vector2<T2> &rhs) {
		using common_t = std::common_type_t<T1, T2>;
		return vector2<common_t>(std::max<common_t>(lhs.x, rhs.x), std::max<common_t>(lhs.y, rhs.y));
	}

	// clamp components of vector between minVal and maxVal
	template <typename T1, typename T2, typename T3>
	inline auto clamp(const vector2<T1> &v, T2 minVal, T3 maxVal) {
		return min(max(v, minVal), maxVal);
	}

	// clamp components of vector between minVal and maxVal components
	template <typename T1, typename T2, typename T3>
	inline auto clamp(const vector2<T1> &v, const vector2<T2> &minVal, const vector2<T3> &maxVal) {
		return min(max(v, minVal), maxVal);
	}

	// linear blend of vectors : x*(1-a) + y*a
	template <typename T1, typename T2, typename T3>
	inline auto mix(const vector2<T1> &lhs, const vector2<T2> &rhs, T3 a) {
		using common_t = std::common_type_t<T1, T2, T3>;
		return lhs*(common_t(1)-a)+rhs*a;
	}

	// linear blend of vectors : x*(1-a) + y*a
	template <typename T1, typename T2, typename T3>
	inline auto mix(const vector2<T1> &lhs, const vector2<T2> &rhs, const vector2<T3> &a) {
		using common_t = std::common_type_t<T1, T2, T3>;
		return lhs*(common_t(1)-a)+rhs*a;
	}

	// 0.0 if edge < v, else 1.0
	template <typename T1, typename T2>
	inline auto step(const vector2<T1> &edge, const vector2<T2> &v) {
		using common_t = std::common_type_t<T1, T2>;
		return vector2<std::common_type_t<T1, T2>>((edge.x<v.x)? common_t(0) : common_t(1), (edge.y<v.y)? common_t(0) : common_t(1));
	}

	// 0.0 if edge < v, else 1.0
	template <typename T1, typename T2>
	inline auto step(T1 edge, const vector2<T2> &v) {
		using common_t = std::common_type_t<T1, T2>;
		return vector2<std::common_type_t<T1, T2>>((edge<v.x)? common_t(0) : common_t(1), (edge<v.y)? common_t(0) : common_t(1));
	}

	// smooth hermit interpolation
	template <typename T1, typename T2, typename T3>
	inline auto smoothstep(const vector2<T1> &edge0, const vector2<T2> &edge1, T3 x) {
		using common_t = std::common_type_t<T1, T2, T3>;
		auto t = clamp((x-edge0)/(edge1-edge0),0, 1);
		return t * t * (common_t(3)-common_t(2)*t);
	}

	// smooth hermit interpolation
	template <typename T1, typename T2, typename T3>
	inline auto smoothstep(const vector2<T1> &edge0, const vector2<T2> &edge1, const vector2<T3> &x) {
		using common_t = std::common_type_t<T1, T2, T3>;
		auto t = clamp((x-edge0)/(edge1-edge0),0, 1);
		return t * t * (common_t(3)-common_t(2)*t);
	}

	// boolean vector of component-wise isnan
	template <typename T>
	inline bvec2 isnan(const vector2<T> &v) {
		return bvec2(isnan(v.x), isnan(v.y));
	}

	// boolean vector of component-wise isinf
	template <typename T>
	inline bvec2 isinf(const vector2<T> &v) {
		return bvec2(isinf(v.x), isinf(v.y));
	}


	// Geometric Functions
	// 

	// length/magnitude of vector
	template <typename T>
	inline T length(const vector2<T> &v) {
		return std::sqrt(v.x * v.x + v.y * v.y);
	}

	// distance between vectors
	template <typename T1, typename T2>
	inline auto distance(const vector2<T1> &lhs, const vector2<T2> &rhs) {
		return length(lhs-rhs);
	}

	// dot product
	template <typename T1, typename T2>
	inline auto dot(const vector2<T1> &lhs, const vector2<T2> &rhs) {
		return lhs.x * rhs.x +  lhs.y * rhs.y;
	}

	// returns unit vector
	template <typename T>
	inline vector2<T> normalize(const vector2<T> &v) {
		return v / length(v);
	}

	// if dot(nref, i) < 0 return n else return -n
	template <typename T1, typename T2, typename T3>
	inline auto faceforward(const vector2<T1> &n, const vector2<T2> &i, const vector2<T3> &nref) {
		using common_t = std::common_type_t<T1, T2, T3>;
		return (dot(nref, i) < common_t(0)) ? n : -n ;
	}

	// for incident i and surface normal n, returns the reflection direction
	template <typename T1, typename T2>
	inline auto reflect(const vector2<T1> &i, const vector2<T2> &n) {
		using common_t = std::common_type_t<T1, T2>;
		return i - common_t(2) * dot(n, i) * n;
	}

	// for incident i, surface normal n, and refraction index eta, return refraction vector
	template <typename T1, typename T2, typename T3>
	inline auto refract(const vector2<T1> &i, const vector2<T2> &n, T3 eta) {
		using common_t = std::common_type_t<T1, T2, T3>;
		auto k = common_t(1) - eta * eta * (common_t(1) - dot(n, i) * dot(n, i));
		if (k < common_t(0)) {
			return vector2<common_t>();
		}
		return eta * i - (eta * dot(n, i) + std::sqrt(k)) * n;
	}


	// Vector Relational Functions
	//

	// component-wise compare of l<r
	template <typename T1, typename T2>
	inline bvec2 lessThan(const vector2<T1> &lhs, const vector2<T2> &rhs) {
		return bvec2(lhs.x < rhs.x, lhs.y < rhs.y);
	}

	// component-wise compare of l<=r
	template <typename T1, typename T2>
	inline bvec2 lessThanEqual(const vector2<T1> &lhs, const vector2<T2> &rhs) {
		return bvec2(lhs.x <= rhs.x, lhs.y <= rhs.y);
	}

	// component-wise compare of l>r
	template <typename T1, typename T2>
	inline bvec2 greaterThan(const vector2<T1> &lhs, const vector2<T2> &rhs) {
		return bvec2(lhs.x < rhs.x, lhs.y < rhs.y);
	}

	// component-wise compare of l>=r
	template <typename T1, typename T2>
	inline bvec2 greaterThanEqual(const vector2<T1> &lhs, const vector2<T2> &rhs) {
		return bvec2(lhs.x >= rhs.x, lhs.y >= rhs.y);
	}

	// component-wise compare of l==r
	template <typename T1, typename T2>
	inline bvec2 equal(const vector2<T1> &lhs, const vector2<T2> &rhs) {
		return bvec2(lhs.x == rhs.x, lhs.y == rhs.y);
	}

	// component-wise compare of l!=r
	template <typename T1, typename T2>
	inline bvec2 notEqual(const vector2<T1> &lhs, const vector2<T2> &rhs) {
		return bvec2(lhs.x != rhs.x, lhs.y != rhs.y);
	}

	// true if ANY of v is true
	inline bool any(const bvec2 &v) {
		return v.x || v.y;
	}

	// true if ANY of v is true
	inline bool all(const bvec2 &v) {
		return v.x && v.y;
	}




	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////    ____  _____   __      ________ _____ _______ ____  _____                                                       ////
	////   |___ \|  __ \  \ \    / /  ____/ ____|__   __/ __ \|  __ \                                                    ////
	////     __) | |  | |  \ \  / /| |__ | |       | | | |  | | |__) |                                                 ////
	////    |__ <| |  | |   \ \/ / |  __|| |       | | | |  | |  _  /                                                ////
	////    ___) | |__| |    \  /  | |___| |____   | | | |__| | | \ \                                                  ////
	////   |____/|_____/      \/   |______\_____|  |_|  \____/|_|  \_\                                                   ////
	////                                                                                                                   ////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	template <typename T>
	class vector3 {
	public:
		union{ T x; T r;};
		union{ T y; T g;};
		union{ T z; T b;};

		// T constructors
		vector3() : x(0), y(0), z(0) {}
		explicit vector3(T v) : x(v), y(v), z(v) {}
		vector3(T _x, T _y, T _z) : x(_x), y(_y), z(_z) {}

		template <typename U>
		vector3(const vector3<U> &other) : x(other.x), y(other.y), z(other.z) {}

		// vector2 constructors
		template <typename U>
		vector3(const vector2<U> &v, T _z) : x(v.x), y(v.y), z(_z) {}

		template <typename U>
		vector3(T _x, const vector2<U> &v) : x(_x), y(v.x), z(v.y) {}

		// vector2 down-cast consctructor
		explicit operator vector2<T>() const {return vector2<T>(x, y);}

		static vector3 random(T lower = 0, T upper = 1) { 
			return vector3(math::random<T>(lower, upper), math::random<T>(lower, upper), math::random<T>(lower, upper));
		}

		static vector3 i() {return vector3(1, 0, 0);}
		static vector3 j() {return vector3(0, 1, 0);}
		static vector3 k() {return vector3(0, 0, 1);}

		static vector3 checknan(const vector3 &v) {
			T sum = v.x + v.y + v.z;
			assert(sum == sum);
			return v;
		}

		explicit operator T *() {
			return &(x);
		}

		T * dataPointer() {
			return &(x);
		}

		const T * dataPointer() const{
			return &(x);
		}

		const T & operator[](size_t i) const {
			assert(i < 3);
			return *(&x + i);
		}

		T & operator[](size_t i) {
			assert(i < 3);
			return *(&x + i);
		}

		// stream insertion
		inline friend std::ostream & operator<<(std::ostream &out, const vector3 &v) {
			return out << '(' << v.x << ", " << v.y << ", " << v.z << ')';
		}

		// Operator overload - assign
		// 

		// assign
		template <typename U>
		vector3 & operator=(const vector3<U> &other) {
			x = other.x;
			y = other.y;
			z = other.z;
			return *this;
		}

		// add-assign
		template <typename U>
		vector3 & operator+=(const vector3<U> &rhs) {
			x += rhs.x;
			y += rhs.y;
			z += rhs.z;
			return *this;
		}

		// add-assign
		vector3 & operator+=(T rhs) {
			x += rhs;
			y += rhs;
			z += rhs;
			return *this;
		}

		// subtract-assign
		template <typename U>
		vector3 & operator-=(const vector3<U> &rhs) {
			x -= rhs.x;
			y -= rhs.y;
			z -= rhs.z;
			return *this;
		}

		// subtract-assign
		vector3 & operator-=(T rhs) {
			x -= rhs;
			y -= rhs;
			z -= rhs;
			return *this;
		}

		// mulitply-assign
		template <typename U>
		vector3 & operator*=(const vector3<U> &rhs) {
			x *= rhs.x;
			y *= rhs.y;
			z *= rhs.z;
			return *this;
		}

		// mulitply-assign
		vector3 & operator*=(T rhs) {
			x *= rhs;
			y *= rhs;
			z *= rhs;
			return *this;
		}

		// divide-assign
		template <typename U>
		vector3 & operator/=(const vector3<U> &rhs) {
			x /= rhs.x;
			y /= rhs.y;
			z /= rhs.z;
			vector3::checknan(*this);
			return *this;
		}

		// divide-assign
		vector3 & operator/=(T rhs) {
			x /= rhs;
			y /= rhs;
			z /= rhs;
			vector3::checknan(*this);
			return *this;
		}
	};


	// Vector / Vector Operator Overloads
	//

	// equality
	template <typename T>
	inline bool operator==(const vector3<T> &lhs, const vector3<T> &rhs) {
		return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z;
	}

	// inequality
	template <typename T>
	inline bool operator!=(const vector3<T> &lhs, const vector3<T> &rhs) {
		return !(lhs == rhs);
	}

	// negate
	template <typename T>
	inline vector3<T> operator-(const vector3<T> &rhs) {
		return vector3<T>(-rhs.x, -rhs.y, -rhs.z);
	}

	// add
	template <typename T1, typename T2>
	inline auto operator+(const vector3<T1> &lhs, const vector3<T2> &rhs) {
		vector3<std::common_type_t<T1, T2>> v = lhs;
		return v += rhs;
	}

	// subtract
	template <typename T1, typename T2>
	inline auto operator-(const vector3<T1> &lhs, const vector3<T2> &rhs) {
		vector3<std::common_type_t<T1, T2>> v = lhs;
		return v -= rhs;
	}

	// multiply
	template <typename T1, typename T2>
	inline auto operator*(const vector3<T1> &lhs, const vector3<T2> &rhs) {
		vector3<std::common_type_t<T1, T2>> v = lhs;
		return v *= rhs;
	}

	// divide
	template <typename T1, typename T2>
	inline auto operator/(const vector3<T1> &lhs, const vector3<T2> &rhs) {
		vector3<std::common_type_t<T1, T2>> v = lhs;
		return v /= rhs;
	}


	// Vector / Scalar Operator Overloads
	//

	// add right
	template <typename T1, typename T2>
	inline auto operator+(const vector3<T1> &lhs, T2 rhs) {
		vector3<std::common_type_t<T1, T2>> v = lhs;
		return v += rhs;
	}

	// add left
	template <typename T1, typename T2>
	inline auto operator+(T1 lhs, const vector3<T2> &rhs) {
		vector3<std::common_type_t<T1, T2>> v = rhs;
		return v += lhs;
	}

	// subtract right
	template <typename T1, typename T2>
	inline auto operator-(const vector3<T1> &lhs, T2 rhs) {
		vector3<std::common_type_t<T1, T2>> v = lhs;
		return v -= rhs;
	}

	// subtract left
	template <typename T1, typename T2>
	inline auto operator-(T1 lhs, const vector3<T2> &rhs) {
		vector3<std::common_type_t<T1, T2>> v(lhs);
		return v -= rhs;
	}

	// multiply right
	template <typename T1, typename T2>
	inline auto operator*(const vector3<T1> &lhs, T2 rhs) {
		vector3<std::common_type_t<T1, T2>> v = lhs;
		return v *= rhs;
	}

	// multiply left
	template <typename T1, typename T2>
	inline auto operator*(T1 lhs, const vector3<T2> &rhs) {
		vector3<std::common_type_t<T1, T2>> v = rhs;
		return v *= lhs;
	}

	// divide right
	template <typename T1, typename T2>
	inline auto operator/(const vector3<T1> &lhs, T2 rhs) {
		vector3<std::common_type_t<T1, T2>> v = lhs;
		return v /= rhs;
	}

	// divide left
	template <typename T1, typename T2>
	inline auto operator/(T1 lhs, const vector3<T2> &rhs) {
		vector3<std::common_type_t<T1, T2>> v(lhs);
		return v /= rhs;
	}


	// Angle and Trigonometry Functions
	// 

	// degrees to radians
	template <typename T>
	inline vector3<T> radians(const vector3<T> &d) {
		return vector3<T>(radians(d.x), radians(d.y), radians(d.z));
	}

	// radians to degrees
	template <typename T>
	inline vector3<T> degrees(const vector3<T> &r) {
		return vector3<T>(degrees(r.x), degrees(r.y), degrees(r.z));
	}

	// sine
	template <typename T>
	inline vector3<T> sin(const vector3<T> &v) {
		return vector3<T>(std::sin(v.x), std::sin(v.y), std::sin(v.z));
	}

	// cosine
	template <typename T>
	inline vector3<T> cos(const vector3<T> &v) {
		return vector3<T>(std::cos(v.x), std::cos(v.y), std::cos(v.z));
	}

	// tangent
	template <typename T>
	inline vector3<T> tan(const vector3<T> &v) {
		return vector3<T>(std::tan(v.x), std::tan(v.y), std::tan(v.z));
	}

	// arc sine
	template <typename T>
	inline vector3<T> asin(const vector3<T> &v) {
		return vector3<T>(std::asin(v.x), std::asin(v.y), std::asin(v.z));
	}

	// arc cosine
	template <typename T>
	inline vector3<T> acos(const vector3<T> &v) {
		return vector3<T>(std::acos(v.x), std::acos(v.y), std::acos(v.z));
	}

	// arc tangent of y/x
	template <typename T1, typename T2>
	inline auto atan(const vector3<T1> &y, const vector3<T2> &x) {
		using common_t = std::common_type_t<T1, T2>;
		return vector3<common_t>(std::atan2<common_t>(y.x, x.x), std::atan2<common_t>(y.y, x.y), std::atan2<common_t>(y.z, x.z));
	}

	// arc tangent
	template <typename T>
	inline vector3<T> atan(const vector3<T> &v) {
		return vector3<T>(std::atan(v.x), std::atan(v.y), std::atan(v.z));
	}


	// Exponential Functions
	// 

	// v raised to the power of e
	template <typename T1, typename T2>
	inline auto pow(const vector3<T1> &v, const vector3<T2> &e) {
		return vector3<std::common_type_t<T1, T2>>(std::pow(v.x, e.x), std::pow(v.y, e.y), std::pow(v.z, e.z));
	}

	// natural exponentiation of vector
	template <typename T>
	inline vector3<T> exp(const vector3<T> &v) {
		return vector3<T>(std::exp(v.x), std::exp(v.y), std::exp(v.z));
	}

	// natural logarithm of vector
	template <typename T>
	inline vector3<T> log(const vector3<T> &v) {
		return vector3<T>(std::log(v.x), std::log(v.y), std::log(v.z));
	}

	// base 2 exponentiation of vector
	template <typename T>
	inline vector3<T> exp2(const vector3<T> &v) {
		return vector3<T>(exp2(v.x), exp2(v.y), exp2(v.z));
	}

	// base 2 logarithm of vector
	template <typename T>
	inline vector3<T> log2(const vector3<T> &v) {
		return vector3<T>(log2(v.x), log2(v.y), log2(v.z));
	}

	// square root of vector
	template <typename T>
	inline vector3<T> sqrt(const vector3<T> &v) {
		return vector3<T>(std::sqrt(v.x), std::sqrt(v.y), std::sqrt(v.z));
	}

	// inverse of the square root of vector
	template <typename T>
	inline vector3<T> inversesqrt(const vector3<T> &v) {
		return vector3<T>::checknan(vector3<T>(T(1)/std::sqrt(v.x), T(1)/std::sqrt(v.y), T(1)/std::sqrt(v.z)));
	}


	// Common Functions
	// 

	// absolute value of vector
	template <typename T>
	inline vector3<T> abs(const vector3<T> &v) {
		return vector3<T>(std::abs(v.x), std::abs(v.y), std::abs(v.z));
	}

	// sign (-1, 0, 1) of vector
	template <typename T>
	inline vector3<T> sign(const vector3<T> &v) {
		return vector3<T>(sign(v.x), sign(v.y), sign(v.z));
	}

	// floor of vector
	template <typename T>
	inline vector3<T> floor(const vector3<T> &v) {
		return vector3<T>(std::floor(v.x), std::floor(v.y), std::floor(v.z));
	}

	// ceil of vector
	template <typename T>
	inline vector3<T> ceil(const vector3<T> &v) {
		return vector3<T>(std::ceil(v.x), std::ceil(v.y), std::ceil(v.z));
	}

	// fractional part of vector : v-floor(v)
	template <typename T>
	inline vector3<T> fract(const vector3<T> &v) {
		return v-floor(v);
	}

	// modulas of vector : v-m*floor(v/m)
	template <typename T>
	inline vector3<T> mod(const vector3<T> &v, T m) {
		return v-m*floor(v/m);
	}

	// modulas of vector : v-m*floor(v/m)
	template <typename T1, typename T2>
	inline auto mod(const vector3<T1> &v, const vector3<T2> &m) {
		return v-m*floor(v/m);
	}

	// minimum of vector components and T
	template <typename T1, typename T2, typename=std::enable_if_t<std::is_arithmetic<T2>::value>>
	inline auto min(const vector3<T1> &lhs, T2 rhs) {
		using common_t = std::common_type_t<T1, T2>;
		return vector3<common_t>(std::min<common_t>(lhs.x, rhs), std::min<common_t>(lhs.y, rhs), std::min<common_t>(lhs.z, rhs));
	}

	// minimum of vector components
	template <typename T1, typename T2>
	inline auto min(const vector3<T1> &lhs, const vector3<T2> &rhs) {
		using common_t = std::common_type_t<T1, T2>;
		return vector3<common_t>(std::min<common_t>(lhs.x, rhs.x), std::min<common_t>(lhs.y, rhs.y), std::min<common_t>(lhs.z, rhs.z));
	}

	// maximum of vector components and T
	template <typename T1, typename T2, typename=std::enable_if_t<std::is_arithmetic<T2>::value>>
	inline auto max(const vector3<T1> &lhs, T2 rhs) {
		using common_t = std::common_type_t<T1, T2>;
		return vector3<common_t>(std::max<common_t>(lhs.x, rhs), std::max<common_t>(lhs.y, rhs), std::max<common_t>(lhs.z, rhs));
	}

	// maximum of vector components
	template <typename T1, typename T2>
	inline auto max(const vector3<T1> &lhs, const vector3<T2> &rhs) {
		using common_t = std::common_type_t<T1, T2>;
		return vector3<common_t>(std::max<common_t>(lhs.x, rhs.x), std::max<common_t>(lhs.y, rhs.y), std::max<common_t>(lhs.z, rhs.z));
	}

	// clamp components of vector between minVal and maxVal
	template <typename T1, typename T2, typename T3>
	inline auto clamp(const vector3<T1> &v, T2 minVal, T3 maxVal) {
		return min(max(v, minVal), maxVal);
	}

	// clamp components of vector between minVal and maxVal components
	template <typename T1, typename T2, typename T3>
	inline auto clamp(const vector3<T1> &v, const vector3<T2> &minVal, const vector3<T3> &maxVal) {
		return min(max(v, minVal), maxVal);
	}

	// linear blend of vectors : x*(1-a) + y*a
	template <typename T1, typename T2, typename T3>
	inline auto mix(const vector3<T1> &lhs, const vector3<T2> &rhs, T3 a) {
		using common_t = std::common_type_t<T1, T2, T3>;
		return lhs*(common_t(1)-a)+rhs*a;
	}

	// linear blend of vectors : x*(1-a) + y*a
	template <typename T1, typename T2, typename T3>
	inline auto mix(const vector3<T1> &lhs, const vector3<T2> &rhs, const vector3<T3> &a) {
		using common_t = std::common_type_t<T1, T2, T3>;
		return lhs*(common_t(1)-a)+rhs*a;
	}

	// 0.0 if edge < v, else 1.0
	template <typename T1, typename T2>
	inline auto step(const vector3<T1> &edge, const vector3<T2> &v) {
		using common_t = std::common_type_t<T1, T2>;
		return vector3<std::common_type_t<T1, T2>>((edge.x<v.x)? common_t(0) : common_t(1), (edge.y<v.y)? common_t(0) : common_t(1), (edge.z<v.z)? common_t(0) : common_t(1));
	}

	// 0.0 if edge < v, else 1.0
	template <typename T1, typename T2>
	inline auto step(T1 edge, const vector3<T2> &v) {
		using common_t = std::common_type_t<T1, T2>;
		return vector3<std::common_type_t<T1, T2>>((edge<v.x)? common_t(0) : common_t(1), (edge<v.y)? common_t(0) : common_t(1), (edge<v.z)? common_t(0) : common_t(1));
	}

	// smooth hermit interpolation
	template <typename T1, typename T2, typename T3>
	inline auto smoothstep(const vector3<T1> &edge0, const vector3<T2> &edge1, T3 x) {
		using common_t = std::common_type_t<T1, T2, T3>;
		auto t = clamp((x-edge0)/(edge1-edge0),0, 1);
		return t * t * (common_t(3)-common_t(2)*t);
	}

	// smooth hermit interpolation
	template <typename T1, typename T2, typename T3>
	inline auto smoothstep(const vector3<T1> &edge0, const vector3<T2> &edge1, const vector3<T3> &x) {
		using common_t = std::common_type_t<T1, T2, T3>;
		auto t = clamp((x-edge0)/(edge1-edge0),0, 1);
		return t * t * (common_t(3)-common_t(2)*t);
	}

	// boolean vector of component-wise isnan
	template <typename T>
	inline bvec3 isnan(const vector3<T> &v) {
		return bvec3(isnan(v.x), isnan(v.y), isnan(v.z));
	}

	// boolean vector of component-wise isinf
	template <typename T>
	inline bvec3 isinf(const vector3<T> &v) {
		return bvec3(isinf(v.x), isinf(v.y), isinf(v.z));
	}


	// Geometric Functions
	// 

	// length/magnitude of vector
	template <typename T>
	inline T length(const vector3<T> &v) {
		return std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
	}

	// distance between vectors
	template <typename T1, typename T2>
	inline auto distance(const vector3<T1> &lhs, const vector3<T2> &rhs) {
		return length(lhs-rhs);
	}

	// dot product
	template <typename T1, typename T2>
	inline auto dot(const vector3<T1> &lhs, const vector3<T2> &rhs) {
		return lhs.x * rhs.x +  lhs.y * rhs.y + lhs.z * rhs.z;
	}

	// cross product
	template <typename T1, typename T2>
	inline  auto cross(const vector3<T1> &lhs, const vector3<T2> &rhs){
		return vector3<std::common_type_t<T1, T2>>(lhs.y * rhs.z - lhs.z * rhs.y, lhs.z * rhs.x - lhs.x * rhs.z, lhs.x * rhs.y - lhs.y * rhs.x);
	}

	// returns unit vector
	template <typename T>
	inline vector3<T> normalize(const vector3<T> &v) {
		return v / length(v);
	}

	// if dot(nref, i) < 0 return n else return -n
	template <typename T1, typename T2, typename T3>
	inline auto faceforward(const vector3<T1> &n, const vector3<T2> &i, const vector3<T3> &nref) {
		using common_t = std::common_type_t<T1, T2, T3>;
		return (dot(nref, i) < common_t(0)) ? n : -n ;
	}

	// for incident i and surface normal n, returns the reflection direction
	template <typename T1, typename T2>
	inline auto reflect(const vector3<T1> &i, const vector3<T2> &n) {
		using common_t = std::common_type_t<T1, T2>;
		return i - common_t(2) * dot(n, i) * n;
	}

	// for incident i, surface normal n, and refraction index eta, return refraction vector
	template <typename T1, typename T2, typename T3>
	inline auto refract(const vector3<T1> &i, const vector3<T2> &n, T3 eta) {
		using common_t = std::common_type_t<T1, T2, T3>;
		auto k = common_t(1) - eta * eta * (common_t(1) - dot(n, i) * dot(n, i));
		if (k < common_t(0)) {
			return vector3<common_t>();
		}
		return eta * i - (eta * dot(n, i) + std::sqrt(k)) * n;
	}


	// Vector Relational Functions
	//

	// component-wise compare of l<r
	template <typename T1, typename T2>
	inline bvec3 lessThan(const vector3<T1> &lhs, const vector3<T2> &rhs) {
		return bvec3(lhs.x < rhs.x, lhs.y < rhs.y, lhs.z < rhs.z);
	}

	// component-wise compare of l<=r
	template <typename T1, typename T2>
	inline bvec3 lessThanEqual(const vector3<T1> &lhs, const vector3<T2> &rhs) {
		return bvec3(lhs.x <= rhs.x, lhs.y <= rhs.y, lhs.z <= rhs.z);
	}

	// component-wise compare of l>r
	template <typename T1, typename T2>
	inline bvec3 greaterThan(const vector3<T1> &lhs, const vector3<T2> &rhs) {
		return bvec3(lhs.x < rhs.x, lhs.y < rhs.y, lhs.z < rhs.z);
	}

	// component-wise compare of l>=r
	template <typename T1, typename T2>
	inline bvec3 greaterThanEqual(const vector3<T1> &lhs, const vector3<T2> &rhs) {
		return bvec3(lhs.x >= rhs.x, lhs.y >= rhs.y, lhs.z >= rhs.z);
	}

	// component-wise compare of l==r
	template <typename T1, typename T2>
	inline bvec3 equal(const vector3<T1> &lhs, const vector3<T2> &rhs) {
		return bvec3(lhs.x == rhs.x, lhs.y == rhs.y, lhs.z == rhs.z);
	}

	// component-wise compare of l!=r
	template <typename T1, typename T2>
	inline bvec3 notEqual(const vector3<T1> &lhs, const vector3<T2> &rhs) {
		return bvec3(lhs.x != rhs.x, lhs.y != rhs.y, lhs.z != rhs.z);
	}

	// true if ANY of v is true
	inline bool any(const bvec3 &v) {
		return v.x || v.y || v.z;
	}

	// true if ANY of v is true
	inline bool all(const bvec3 &v) {
		return v.x && v.y && v.z;
	}



	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////    _  _   _____   __      ________ _____ _______ ____  _____                                                      ////
	////   | || | |  __ \  \ \    / /  ____/ ____|__   __/ __ \|  __ \                                                   ////
	////   | || |_| |  | |  \ \  / /| |__ | |       | | | |  | | |__) |                                                ////
	////   |__   _| |  | |   \ \/ / |  __|| |       | | | |  | |  _  /                                               ////
	////      | | | |__| |    \  /  | |___| |____   | | | |__| | | \ \                                                 ////
	////      |_| |_____/      \/   |______\_____|  |_|  \____/|_|  \_\                                                  ////
	////                                                                                                                   ////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	template <typename T>
	class vector4 {
	public:
		union{ T x; T r;};
		union{ T y; T g;};
		union{ T z; T b;};
		union{ T w; T a;};

		// T constructors
		vector4() : x(0), y(0), z(0), w(0) {}
		explicit vector4(T v) : x(v), y(v), z(v), w(v) {}
		vector4(T _x, T _y, T _z, T _w) : x(_x), y(_y), z(_z), w(_w) {}

		template <typename U>
		vector4(const vector2<U> &other) : x(other.x), y(other.y), z(other.z), w(other.w) {}

		// vector2 constructors
		template <typename U>
		vector4(const vector2<U> &v, T _z, T _w) : x(v.x), y(v.y), z(_z), w(_w) {}
		template <typename U>
		vector4(T _x, const vector2<U> &v, T _w) : x(_x), y(v.x), z(v.y), w(_w) {}
		template <typename U>
		vector4(T _x, T _y, const vector2<U> &v) : x(_x), y(_y), z(v.x), w(v.y) {}
		template <typename U1, typename U2>
		vector4(const vector2<U1> &v0, const vector2<U2> &v1) : x(v0.x), y(v0.y), z(v1.x), w(v1.y) {}

		// vector3 constructors
		template <typename U>
		vector4(const vector3<U> &v, T _w) : x(v.x), y(v.y), z(v.z), w(_w) {}
		template <typename U>
		vector4(T _x, const vector3<U> &v) : x(_x), y(v.x), z(v.y), w(v.z) {}

		// vector2, vector3 down-cast constructors
		explicit operator vector2<T>() const {return vector2<T>(x, y);}
		explicit operator vector3<T>() const {return vector3<T>(x, y, z);}

		static vector4 random(T lower = 0, T upper = 1) { 
			return vector4(math::random<T>(lower, upper), math::random<T>(lower, upper),
				math::random<T>(lower, upper), math::random<T>(lower, upper));
		}

		static vector4 i() {return vector4(1, 0, 0, 0);}
		static vector4 j() {return vector4(0, 1, 0, 0);}
		static vector4 k() {return vector4(0, 0, 1, 0);}
		static vector4 l() {return vector4(0, 0, 0, 1);}

		static vector4 checknan(const vector4 &v) {
			T sum = v.x + v.y + v.z + v.w;
			assert(sum == sum);
			return v;
		}

		explicit operator T *() {
			return &(x);
		}

		T * dataPointer() {
			return &(x);
		}

		const T * dataPointer() const{
			return &(x);
		}

		const T & operator[](size_t i) const {
			assert(i < 4);
			return *(&x + i);
		}

		T & operator[](size_t i) {
			assert(i < 4);
			return *(&x + i);
		}

		// stream insertion
		inline friend std::ostream & operator<<(std::ostream &out, const vector4 &v) {
			return out << '(' << v.x << ", " << v.y << ", " << v.z << ", " << v.w << ')';
		}

		// Operator overload - assign
		// 

		// assign
		template <typename U>
		vector4 & operator=(const vector4<U> &other) {
			x = other.x;
			y = other.y;
			z = other.z;
			w = other.w;
			return *this;
		}

		// add-assign
		template <typename U>
		vector4 & operator+=(const vector4<U> &rhs) {
			x += rhs.x;
			y += rhs.y;
			z += rhs.z;
			w += rhs.w;
			return *this;
		}

		// add-assign
		vector4 & operator+=(T rhs) {
			x += rhs;
			y += rhs;
			z += rhs;
			w += rhs;
			return *this;
		}

		// subtract-assign
		template <typename U>
		vector4 & operator-=(const vector4<U> &rhs) {
			x -= rhs.x;
			y -= rhs.y;
			z -= rhs.z;
			w -= rhs.w;
			return *this;
		}

		// subtract-assign
		vector4 & operator-=(T rhs) {
			x -= rhs;
			y -= rhs;
			z -= rhs;
			w -= rhs;
			return *this;
		}

		// mulitply-assign
		template <typename U>
		vector4 & operator*=(const vector4<U> &rhs) {
			x *= rhs.x;
			y *= rhs.y;
			z *= rhs.z;
			w *= rhs.w;
			return *this;
		}

		// mulitply-assign
		vector4 & operator*=(T rhs) {
			x *= rhs;
			y *= rhs;
			z *= rhs;
			w *= rhs;
			return *this;
		}

		// divide-assign
		template <typename U>
		vector4 & operator/=(const vector4<U> &rhs) {
			x /= rhs.x;
			y /= rhs.y;
			z /= rhs.z;
			w /= rhs.w;
			vector4::checknan(*this);
			return *this;
		}

		// divide-assign
		vector4 & operator/=(T rhs) {
			x /= rhs;
			y /= rhs;
			z /= rhs;
			w /= rhs;
			vector4::checknan(*this);
			return *this;
		}
	};



	// Vector / Vector Operator Overloads
	//

	// equality
	template <typename T>
	inline bool operator==(const vector4<T> &lhs, const vector4<T> &rhs) {
		return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z && lhs.w == rhs.w;
	}

	// inequality
	template <typename T>
	inline bool operator!=(const vector4<T> &lhs, const vector4<T> &rhs) {
		return !(lhs == rhs);
	}

	// negate
	template <typename T>
	inline vector4<T> operator-(const vector4<T> &rhs) {
		return vector4<T>(-rhs.x, -rhs.y, -rhs.z, -rhs.w);
	}

	// add
	template <typename T1, typename T2>
	inline auto operator+(const vector4<T1> &lhs, const vector4<T2> &rhs) {
		vector4<std::common_type_t<T1, T2>> v = lhs;
		return v += rhs;
	}

	// subtract
	template <typename T1, typename T2>
	inline auto operator-(const vector4<T1> &lhs, const vector4<T2> &rhs) {
		vector4<std::common_type_t<T1, T2>> v = lhs;
		return v -= rhs;
	}

	// multiply
	template <typename T1, typename T2>
	inline auto operator*(const vector4<T1> &lhs, const vector4<T2> &rhs) {
		vector4<std::common_type_t<T1, T2>> v = lhs;
		return v *= rhs;
	}

	// divide
	template <typename T1, typename T2>
	inline auto operator/(const vector4<T1> &lhs, const vector4<T2> &rhs) {
		vector4<std::common_type_t<T1, T2>> v = lhs;
		return v /= rhs;
	}


	// Vector / Scalar Operator Overloads
	//

	// add right
	template <typename T1, typename T2>
	inline auto operator+(const vector4<T1> &lhs, T2 rhs) {
		vector4<std::common_type_t<T1, T2>> v = lhs;
		return v += rhs;
	}

	// add left
	template <typename T1, typename T2>
	inline auto operator+(T1 lhs, const vector4<T2> &rhs) {
		vector4<std::common_type_t<T1, T2>> v = rhs;
		return v += lhs;
	}

	// subtract right
	template <typename T1, typename T2>
	inline auto operator-(const vector4<T1> &lhs, T2 rhs) {
		vector4<std::common_type_t<T1, T2>> v = lhs;
		return v -= rhs;
	}

	// subtract left
	template <typename T1, typename T2>
	inline auto operator-(T1 lhs, const vector4<T2> &rhs) {
		vector4<std::common_type_t<T1, T2>> v(lhs);
		return v -= rhs;
	}

	// multiply right
	template <typename T1, typename T2>
	inline auto operator*(const vector4<T1> &lhs, T2 rhs) {
		vector4<std::common_type_t<T1, T2>> v = lhs;
		return v *= rhs;
	}

	// multiply left
	template <typename T1, typename T2>
	inline auto operator*(T1 lhs, const vector4<T2> &rhs) {
		vector4<std::common_type_t<T1, T2>> v = rhs;
		return v *= lhs;
	}

	// divide right
	template <typename T1, typename T2>
	inline auto operator/(const vector4<T1> &lhs, T2 rhs) {
		vector4<std::common_type_t<T1, T2>> v = lhs;
		return v /= rhs;
	}

	// divide left
	template <typename T1, typename T2>
	inline auto operator/(T1 lhs, const vector4<T2> &rhs) {
		vector4<std::common_type_t<T1, T2>> v(lhs);
		return v /= rhs;
	}


	// Angle and Trigonometry Functions
	// 

	// degrees to radians
	template <typename T>
	inline vector4<T> radians(const vector4<T> &d) {
		return vector4<T>(radians(d.x), radians(d.y), radians(d.z), radians(d.w));
	}

	// radians to degrees
	template <typename T>
	inline vector4<T> degrees(const vector4<T> &r) {
		return vector4<T>(degrees(r.x), degrees(r.y), degrees(r.z), degrees(r.w));
	}

	// sine
	template <typename T>
	inline vector4<T> sin(const vector4<T> &v) {
		return vector4<T>(std::sin(v.x), std::sin(v.y), std::sin(v.z), std::sin(v.w));
	}

	// cosine
	template <typename T>
	inline vector4<T> cos(const vector4<T> &v) {
		return vector4<T>(std::cos(v.x), std::cos(v.y), std::cos(v.z), std::cos(v.w));
	}

	// tangent
	template <typename T>
	inline vector4<T> tan(const vector4<T> &v) {
		return vector4<T>(std::tan(v.x), std::tan(v.y), std::tan(v.z), std::tan(v.w));
	}

	// arc sine
	template <typename T>
	inline vector4<T> asin(const vector4<T> &v) {
		return vector4<T>(std::asin(v.x), std::asin(v.y), std::asin(v.z), std::asin(v.w));
	}

	// arc cosine
	template <typename T>
	inline vector4<T> acos(const vector4<T> &v) {
		return vector4<T>(std::acos(v.x), std::acos(v.y), std::acos(v.z), std::acos(v.w));
	}

	// arc tangent of y/x
	template <typename T1, typename T2>
	inline auto atan(const vector4<T1> &y, const vector4<T2> &x) {
		using common_t = std::common_type_t<T1, T2>;
		return vector4<common_t>(std::atan2<common_t>(y.x, x.x), std::atan2<common_t>(y.y, x.y), std::atan2<common_t>(y.z, x.z), std::atan2<common_t>(y.w, x.w));
	}

	// arc tangent
	template <typename T>
	inline vector4<T> atan(const vector4<T> &v) {
		return vector4<T>(std::atan(v.x), std::atan(v.y), std::atan(v.z), std::atan(v.w));
	}


	// Exponential Functions
	// 

	// v raised to the power of e
	template <typename T1, typename T2>
	inline auto pow(const vector4<T1> &v, const vector4<T2> &e) {
		return vector4<std::common_type_t<T1, T2>>(std::pow(v.x, e.x), std::pow(v.y, e.y), std::pow(v.z, e.z), std::pow(v.w, e.w));
	}

	// natural exponentiation of vector
	template <typename T>
	inline vector4<T> exp(const vector4<T> &v) {
		return vector4<T>(std::exp(v.x), std::exp(v.y), std::exp(v.z), std::exp(v.w));
	}

	// natural logarithm of vector
	template <typename T>
	inline vector4<T> log(const vector4<T> &v) {
		return vector4<T>(std::log(v.x), std::log(v.y), std::log(v.z), std::log(v.w));
	}

	// base 2 exponentiation of vector
	template <typename T>
	inline vector4<T> exp2(const vector4<T> &v) {
		return vector4<T>(exp2(v.x), exp2(v.y), exp2(v.z), exp2(v.w));
	}

	// base 2 logarithm of vector
	template <typename T>
	inline vector4<T> log2(const vector4<T> &v) {
		return vector4<T>(log2(v.x), log2(v.y), log2(v.z), log2(v.w));
	}

	// square root of vector
	template <typename T>
	inline vector4<T> sqrt(const vector4<T> &v) {
		return vector4<T>(std::sqrt(v.x), std::sqrt(v.y), std::sqrt(v.z), std::sqrt(v.w));
	}

	// inverse of the square root of vector
	template <typename T>
	inline vector4<T> inversesqrt(const vector4<T> &v) {
		return vector4<T>::checknan(vector4<T>(T(1)/std::sqrt(v.x), T(1)/std::sqrt(v.y), T(1)/std::sqrt(v.z), T(1)/std::sqrt(v.w)));
	}


	// Common Functions
	// 

	// absolute value of vector
	template <typename T>
	inline vector4<T> abs(const vector4<T> &v) {
		return vector4<T>(std::abs(v.x), std::abs(v.y), std::abs(v.z), std::abs(v.w));
	}

	// sign (-1, 0, 1) of vector
	template <typename T>
	inline vector4<T> sign(const vector4<T> &v) {
		return vector4<T>(sign(v.x), sign(v.y), sign(v.z), sign(v.w));
	}

	// floor of vector
	template <typename T>
	inline vector4<T> floor(const vector4<T> &v) {
		return vector4<T>(std::floor(v.x), std::floor(v.y), std::floor(v.z), std::floor(v.w));
	}

	// ceil of vector
	template <typename T>
	inline vector4<T> ceil(const vector4<T> &v) {
		return vector4<T>(std::ceil(v.x), std::ceil(v.y), std::ceil(v.z), std::ceil(v.w));
	}

	// fractional part of vector : v-floor(v)
	template <typename T>
	inline vector4<T> fract(const vector4<T> &v) {
		return v-floor(v);
	}

	// modulas of vector : v-m*floor(v/m)
	template <typename T>
	inline vector4<T> mod(const vector4<T> &v, T m) {
		return v-m*floor(v/m);
	}

	// modulas of vector : v-m*floor(v/m)
	template <typename T1, typename T2>
	inline auto mod(const vector4<T1> &v, const vector4<T2> &m) {
		return v-m*floor(v/m);
	}

	// minimum of vector components and T
	template <typename T1, typename T2, typename=std::enable_if_t<std::is_arithmetic<T2>::value>>
	inline auto min(const vector4<T1> &lhs, T2 rhs) {
		using common_t = std::common_type_t<T1, T2>;
		return vector4<common_t>(std::min<common_t>(lhs.x, rhs), std::min<common_t>(lhs.y, rhs), std::min<common_t>(lhs.z, rhs), std::min<common_t>(lhs.w, rhs));
	}

	// minimum of vector components
	template <typename T1, typename T2>
	inline auto min(const vector4<T1> &lhs, const vector4<T2> &rhs) {
		using common_t = std::common_type_t<T1, T2>;
		return vector4<common_t>(std::min<common_t>(lhs.x, rhs.x), std::min<common_t>(lhs.y, rhs.y), std::min<common_t>(lhs.z, rhs.z), std::min<common_t>(lhs.w, rhs.w));
	}

	// maximum of vector components and T
	template <typename T1, typename T2, typename=std::enable_if_t<std::is_arithmetic<T2>::value>>
	inline auto max(const vector4<T1> &lhs, T2 rhs) {
		using common_t = std::common_type_t<T1, T2>;
		return vector4<common_t>(std::max<common_t>(lhs.x, rhs), std::max<common_t>(lhs.y, rhs), std::max<common_t>(lhs.z, rhs), std::max<common_t>(lhs.w, rhs));
	}

	// maximum of vector components
	template <typename T1, typename T2>
	inline auto max(const vector4<T1> &lhs, const vector4<T2> &rhs) {
		using common_t = std::common_type_t<T1, T2>;
		return vector4<common_t>(std::max<common_t>(lhs.x, rhs.x), std::max<common_t>(lhs.y, rhs.y), std::max<common_t>(lhs.z, rhs.z), std::max<common_t>(lhs.w, rhs.w));
	}

	// clamp components of vector between minVal and maxVal
	template <typename T1, typename T2, typename T3>
	inline auto clamp(const vector4<T1> &v, T2 minVal, T3 maxVal) {
		return min(max(v, minVal), maxVal);
	}

	// clamp components of vector between minVal and maxVal components
	template <typename T1, typename T2, typename T3>
	inline auto clamp(const vector4<T1> &v, const vector4<T2> &minVal, const vector4<T3> &maxVal) {
		return min(max(v, minVal), maxVal);
	}

	// linear blend of vectors : x*(1-a) + y*a
	template <typename T1, typename T2, typename T3>
	inline auto mix(const vector4<T1> &lhs, const vector4<T2> &rhs, T3 a) {
		using common_t = std::common_type_t<T1, T2, T3>;
		return lhs*(common_t(1)-a)+rhs*a;
	}

	// linear blend of vectors : x*(1-a) + y*a
	template <typename T1, typename T2, typename T3>
	inline auto mix(const vector4<T1> &lhs, const vector4<T2> &rhs, const vector4<T3> &a) {
		using common_t = std::common_type_t<T1, T2, T3>;
		return lhs*(common_t(1)-a)+rhs*a;
	}

	// 0.0 if edge < v, else 1.0
	template <typename T1, typename T2>
	inline auto step(const vector4<T1> &edge, const vector4<T2> &v) {
		using common_t = std::common_type_t<T1, T2>;
		return vector4<std::common_type_t<T1, T2>>((edge.x<v.x)? common_t(0) : common_t(1), (edge.y<v.y)? common_t(0) : common_t(1), (edge.z<v.z)? common_t(0) : common_t(1), (edge.w<v.w)? common_t(0) : common_t(1));
	}

	// 0.0 if edge < v, else 1.0
	template <typename T1, typename T2>
	inline auto step(T1 edge, const vector4<T2> &v) {
		using common_t = std::common_type_t<T1, T2>;
		return vector4<std::common_type_t<T1, T2>>((edge<v.x)? common_t(0) : common_t(1), (edge<v.y)? common_t(0) : common_t(1), (edge<v.z)? common_t(0) : common_t(1), (edge<v.w)? common_t(0) : common_t(1));
	}

	// smooth hermit interpolation
	template <typename T1, typename T2, typename T3>
	inline auto smoothstep(const vector4<T1> &edge0, const vector4<T2> &edge1, T3 x) {
		using common_t = std::common_type_t<T1, T2, T3>;
		auto t = clamp((x-edge0)/(edge1-edge0),0, 1);
		return t * t * (common_t(3)-common_t(2)*t);
	}

	// smooth hermit interpolation
	template <typename T1, typename T2, typename T3>
	inline auto smoothstep(const vector4<T1> &edge0, const vector4<T2> &edge1, const vector4<T3> &x) {
		using common_t = std::common_type_t<T1, T2, T3>;
		auto t = clamp((x-edge0)/(edge1-edge0),0, 1);
		return t * t * (common_t(3)-common_t(2)*t);
	}

	// boolean vector of component-wise isnan
	template <typename T>
	inline bvec4 isnan(const vector4<T> &v) {
		return bvec4(isnan(v.x), isnan(v.y), isnan(v.z), isnan(v.w));
	}

	// boolean vector of component-wise isinf
	template <typename T>
	inline bvec4 isinf(const vector4<T> &v) {
		return bvec4(isinf(v.x), isinf(v.y), isinf(v.z), isinf(v.w));
	}


	// Geometric Functions
	// 

	// length/magnitude of vector
	template <typename T>
	inline T length(const vector4<T> &v) {
		return std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z + v.w * v.w);
	}

	// distance between vectors
	template <typename T1, typename T2>
	inline auto distance(const vector4<T1> &lhs, const vector4<T2> &rhs) {
		return length(lhs-rhs);
	}

	// dot product
	template <typename T1, typename T2>
	inline auto dot(const vector4<T1> &lhs, const vector4<T2> &rhs) {
		return lhs.x * rhs.x +  lhs.y * rhs.y + lhs.z * rhs.z + lhs.w * rhs.w;
	}

	// returns unit vector
	template <typename T>
	inline vector4<T> normalize(const vector4<T> &v) {
		return v / length(v);
	}

	// if dot(nref, i) < 0 return n else return -n
	template <typename T1, typename T2, typename T3>
	inline auto faceforward(const vector4<T1> &n, const vector4<T2> &i, const vector4<T3> &nref) {
		using common_t = std::common_type_t<T1, T2, T3>;
		return (dot(nref, i) < common_t(0)) ? n : -n ;
	}

	// for incident i and surface normal n, returns the reflection direction
	template <typename T1, typename T2>
	inline auto reflect(const vector4<T1> &i, const vector4<T2> &n) {
		using common_t = std::common_type_t<T1, T2>;
		return i - common_t(2) * dot(n, i) * n;
	}

	// for incident i, surface normal n, and refraction index eta, return refraction vector
	template <typename T1, typename T2, typename T3>
	inline auto refract(const vector4<T1> &i, const vector4<T2> &n, T3 eta) {
		using common_t = std::common_type_t<T1, T2, T3>;
		auto k = common_t(1) - eta * eta * (common_t(1) - dot(n, i) * dot(n, i));
		if (k < common_t(0)) {
			return vector4<common_t>();
		}
		return eta * i - (eta * dot(n, i) + std::sqrt(k)) * n;
	}


	// Vector Relational Functions
	//

	// component-wise compare of l<r
	template <typename T1, typename T2>
	inline bvec4 lessThan(const vector4<T1> &lhs, const vector4<T2> &rhs) {
		return bvec4(lhs.x < rhs.x, lhs.y < rhs.y, lhs.z < rhs.z, lhs.w < rhs.w);
	}

	// component-wise compare of l<=r
	template <typename T1, typename T2>
	inline bvec4 lessThanEqual(const vector4<T1> &lhs, const vector4<T2> &rhs) {
		return bvec4(lhs.x <= rhs.x, lhs.y <= rhs.y, lhs.z <= rhs.z, lhs.w <= rhs.w);
	}

	// component-wise compare of l>r
	template <typename T1, typename T2>
	inline bvec4 greaterThan(const vector4<T1> &lhs, const vector4<T2> &rhs) {
		return bvec4(lhs.x < rhs.x, lhs.y < rhs.y, lhs.z < rhs.z, lhs.w < rhs.w);
	}

	// component-wise compare of l>=r
	template <typename T1, typename T2>
	inline bvec4 greaterThanEqual(const vector4<T1> &lhs, const vector4<T2> &rhs) {
		return bvec4(lhs.x >= rhs.x, lhs.y >= rhs.y, lhs.z >= rhs.z, lhs.w >= rhs.w);
	}

	// component-wise compare of l==r
	template <typename T1, typename T2>
	inline bvec4 equal(const vector4<T1> &lhs, const vector4<T2> &rhs) {
		return bvec4(lhs.x == rhs.x, lhs.y == rhs.y, lhs.z == rhs.z, lhs.w == rhs.w);
	}

	// component-wise compare of l!=r
	template <typename T1, typename T2>
	inline bvec4 notEqual(const vector4<T1> &lhs, const vector4<T2> &rhs) {
		return bvec4(lhs.x != rhs.x, lhs.y != rhs.y, lhs.z != rhs.z, lhs.w != rhs.w);
	}

	// true if ANY of v is true
	inline bool any(const bvec4 &v) {
		return v.x || v.y || v.z || v.w;
	}

	// true if ANY of v is true
	inline bool all(const bvec4 &v) {
		return v.x && v.y && v.z && v.w;
	}




	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////  ___       ___    __  __       _______ _____  _______   __                                                        ////
	//// |__ \     |__ \  |  \/  |   /\|__   __|  __ \|_   _\ \ / /                                                      ////
	////    ) |_  __  ) | | \  / |  /  \  | |  | |__) | | |  \ V /                                                     ////
	////   / /\ \/ / / /  | |\/| | / /\ \ | |  |  _  /  | |   > <                                                    ////
	////  / /_ >  < / /_  | |  | |/ ____ \| |  | | \ \ _| |_ / . \                                                     ////
	//// |____/_/\_\____| |_|  |_/_/    \_\_|  |_|  \_\_____/_/ \_\                                                      ////
	////                                                                                                                   ////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	template <typename T>
	class matrix2 {
	private:
		vector2<T> data[2];

		void clear(T i=0) {
			data[0] = vector2<T>(i, 0);
			data[1] = vector2<T>(0, i);
		}

	public:
		matrix2() {
			clear();
		}

		explicit matrix2(T i) {
			clear(i);
		}

		template <typename U1, typename U2>
		matrix2(const vector2<U1> &a, const vector2<U2> &b) {
			data[0] = a;
			data[1] = b;
		}

		matrix2(
			T e00, T e10, 
			T e01, T e11
		) {
			data[0] = vector2<T>(e00, e10);
			data[1] = vector2<T>(e01, e11);
		}

		template <typename U>
		matrix2(const matrix2<U> &other) {
			data[0] = other.data[0];
			data[1] = other.data[1];
		}

		static matrix2 random(T lower = 0, T upper = 1) { 
			return matrix2(vector2<T>::random(lower, upper), vector2<T>::random(lower, upper));
		}

		static matrix2 identity() {return matrix2(1);}

		explicit operator T *() {
			return &(data[0].x);
		}

		T * dataPointer() {
			return &(data[0].x);
		}

		const T * dataPointer() const{
			return &(data[0].x);
		}

		const vector2<T> & operator[](size_t i) const {
			assert(i < 2);
			return data[i];
		}

		vector2<T> & operator[](size_t i) {
			assert(i < 2);
			return data[i];
		}

		// stream insertion
		inline friend std::ostream & operator<<(std::ostream &out, const matrix2 &m) {
			const size_t field_width = 10;
			std::ostringstream oss;
			oss << std::setprecision(4);
			oss << '[' << std::setw(field_width) << m[0][0] << ", " << std::setw(field_width) << m[1][0] << ']' << std::endl;
			oss << '[' << std::setw(field_width) << m[0][1] << ", " << std::setw(field_width) << m[1][1] << ']' << std::endl;
			return out << oss.str();
		}

		// Operator overload - assign
		// 

		// assign
		template <typename U>
		matrix2 & operator=(const matrix2<U> &other) {
			data[0] = other.data[0];
			data[1] = other.data[1];
			return *this;
		}

		// add-assign
		template <typename U>
		matrix2 & operator+=(const matrix2<U> &rhs) {
			data[0] += rhs[0];
			data[1] += rhs[1];
			return *this;
		}

		// add-assign
		matrix2 & operator+=(T rhs) {
			data[0] += rhs;
			data[1] += rhs;
			return *this;
		}

		// subtract-assign
		template <typename U>
		matrix2 & operator-=(const matrix2<U> &rhs) {
			data[0] -= rhs[0];
			data[1] -= rhs[1];
			return *this;
		}

		// subtract-assign
		matrix2 & operator-=(T rhs) {
			data[0] -= rhs;
			data[1] -= rhs;
			return *this;
		}

		// This is a special case where we use matrix product
		// if you want component-wise multiplication see matrixCompMult(matrix2, matrix2)
		// 
		// mulitply-assign
		template <typename U>
		matrix2 & operator*=(const matrix2<U> &rhs) {
			const matrix2 lhs = *this;
			for (int i = 0; i < 2; i++) {
				(*this)[i] = vector2<T>(0);
				for (int j = 0; j < 2; j++) {
					(*this)[i] += lhs[j] * rhs[i][j];
				}
			}
			return *this;
		}

		// mulitply-assign
		matrix2 & operator*=(T rhs) {
			data[0] *= rhs;
			data[1] *= rhs;
			return *this;
		}

		// divide-assign
		template <typename U>
		matrix2 & operator/=(const matrix2<U> &rhs) {
			data[0] /= rhs[0];
			data[1] /= rhs[1];
			return *this;
		}

		// divide-assign
		matrix2 & operator/=(T rhs) {
			data[0] /= rhs;
			data[1] /= rhs;
			return *this;
		}

	};

	// Matrix / Matrix Operator Overloads
	// Matrix / Vector Operator Overloads
	//

	// negate
	template <typename T>
	inline matrix2<T> operator-(const matrix2<T> &rhs) {
		matrix2<T> m;
		m[0] = -rhs[0];
		m[1] = -rhs[1];
		return m;
	}

	// add
	template <typename T1, typename T2>
	inline auto operator+(const matrix2<T1> &lhs, const matrix2<T2> &rhs) {
		matrix2<std::common_type_t<T1, T2>> m = lhs;
		return m += rhs;
	}

	// subtract
	template <typename T1, typename T2>
	inline auto operator-(const matrix2<T1> &lhs, const matrix2<T2> &rhs) {
		matrix2<std::common_type_t<T1, T2>> m = lhs;
		return m -= rhs;
	}

	// This is a special case where we use matrix product
	// if you want component-wise multiplication see matrixCompMult(matrix2<T>, matrix2<T>)
	// 
	// multiply
	template <typename T1, typename T2>
	inline auto operator*(const matrix2<T1> &lhs, const matrix2<T2> &rhs) {
		matrix2<std::common_type_t<T1, T2>> m = lhs;
		return m *= rhs;
	}

	// Left multiply matrix2<T> m with vector2<T> v
	// 
	// multiply
	template <typename T1, typename T2>
	inline auto operator*(const matrix2<T1> &lhs, const vector2<T2> &rhs) {
		vector2<std::common_type_t<T1, T2>> v(0);
		for (int j = 0; j < 2; j++) {
			v += lhs[j] * rhs[j];
		}
		return v;
	}

	// Right multiply vector2<T> v with matrix2<T> m
	// Equvilent to transpose(m) * v
	// 
	// mulitply-assign
	template <typename T1, typename T2>
	inline auto operator*=(vector2<T1> &lhs, const matrix2<T2> &rhs) {
		return lhs = vector2<std::common_type_t<T1, T2>>(dot(lhs, rhs[0]), dot(lhs, rhs[1]));
	}

	// multiply
	template <typename T1, typename T2>
	inline auto operator*(const vector2<T1> &lhs, const matrix2<T2> &rhs) {
		vector2<std::common_type_t<T1, T2>> v = lhs;
		return v *= rhs;
	}

	// divide
	template <typename T1, typename T2>
	inline auto operator/(const matrix2<T1> &lhs, const matrix2<T2> &rhs) {
		matrix2<std::common_type_t<T1, T2>> m = lhs;
		return m /= rhs;
	}


	// Vector / Scalar Operator Overloads
	//

	// add right
	template <typename T1, typename T2>
	inline auto operator+(const matrix2<T1> &lhs, T2 rhs) {
		matrix2<std::common_type_t<T1, T2>> m = lhs;
		return m += rhs;
	}

	// add left
	template <typename T1, typename T2>
	inline auto operator+(T1 lhs, const matrix2<T2> &rhs) {
		matrix2<std::common_type_t<T1, T2>> m = rhs;
		return m += lhs;
	}

	// subtract right
	template <typename T1, typename T2>
	inline auto operator-(const matrix2<T1> &lhs, T2 rhs) {
		matrix2<std::common_type_t<T1, T2>> m = lhs;
		return m -= rhs;
	}

	// subtract left
	template <typename T1, typename T2>
	inline auto operator-(T1 lhs, const matrix2<T2> &rhs) {
		vector2<std::common_type_t<T1, T2>> v(lhs);
		matrix2<std::common_type_t<T1, T2>> m(v, v);
		return m -= rhs;
	}

	// multiply right
	template <typename T1, typename T2>
	inline auto operator*(const matrix2<T1> &lhs, T2 rhs) {
		matrix2<std::common_type_t<T1, T2>> m = lhs;
		return m *= rhs;
	}

	// multiply left
	template <typename T1, typename T2>
	inline auto operator*(T1 lhs, const matrix2<T2> &rhs) {
		matrix2<std::common_type_t<T1, T2>> m = rhs;
		return m *= lhs;
	}

	// divide right
	template <typename T1, typename T2>
	inline auto operator/(const matrix2<T1> &lhs, T2 rhs) {
		matrix2<std::common_type_t<T1, T2>> m = lhs;
		return m /= rhs;
	}

	// divide left
	template <typename T1, typename T2>
	inline auto operator/(T1 lhs, const matrix2<T2> &rhs) {
		vector2<std::common_type_t<T1, T2>> v(lhs);
		matrix2<std::common_type_t<T1, T2>> m(v, v);
		return m /= rhs;
	}

	// Matrix functions
	// 

	// determinant of matrix
	template <typename T>
	inline T determinant(const matrix2<T> &m) {
		return m[0].x * m[1].y - m[1].x * m[0].y;
	}

	// inverse of matrix (error if not invertible)
	template <typename T>
	inline matrix2<T> inverse(const matrix2<T> &m) {
		T invdet = 1 / determinant(m);
		// FIXME proper detect infinite determinant
		assert(!isinf(invdet) && invdet == invdet && invdet != 0);
		return invdet * matrix2<T>(m[1].y, -m[0].y, -m[1].x, m[0].x);
	}

	// transpose of matrix
	template <typename T>
	inline matrix2<T> transpose(const matrix2<T> &m) {
		return matrix2<T>(
			m[0][0], m[1][0],
			m[0][1], m[1][1]
		);
	}

	// component-wise multiplication 
	// see (*) and (*=) operator overload for matrix product
	template <typename T1, typename T2>
	inline auto matrixCompMult(const matrix2<T1> &lhs, const matrix2<T2> &rhs) {
		matrix2<std::common_type_t<T1, T2>> m = lhs;
		m[0] *= rhs[0];
		m[1] *= rhs[1];
		return m;
	}

	template <typename T1, typename T2>
	inline auto outerProduct(const vector2<T1> &lhs, const vector2<T2> &rhs){
		matrix2<std::common_type_t<T1, T2>> m;
		m[0] = lhs * rhs[0];
		m[1] = lhs * rhs[1];
		return m;
	}





	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////  ____       ____    __  __       _______ _____  _______   __                                                      ////
	//// |___ \     |___ \  |  \/  |   /\|__   __|  __ \|_   _\ \ / /                                                    ////
	////   __) |_  __ __) | | \  / |  /  \  | |  | |__) | | |  \ V /                                                   ////
	////  |__ <\ \/ /|__ <  | |\/| | / /\ \ | |  |  _  /  | |   > <                                                  ////
	////  ___) |>  < ___) | | |  | |/ ____ \| |  | | \ \ _| |_ / . \                                                   ////
	//// |____//_/\_\____/  |_|  |_/_/    \_\_|  |_|  \_\_____/_/ \_\                                                    ////
	////                                                                                                                   ////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	template <typename T>
	class matrix3 {
	private:
		vector3<T> data[3];

		void clear(T i=0) {
			data[0] = vector3<T>(i, 0, 0);
			data[1] = vector3<T>(0, i, 0);
			data[2] = vector3<T>(0, 0, i);
		}

	public:
		matrix3() {
			clear();
		}

		explicit matrix3(T i) {
			clear(i);
		}

		template <typename U1, typename U2, typename U3>
		matrix3(const vector3<U1> &a, const vector3<U2> &b, const vector3<U3> &c) {
			data[0] = a;
			data[1] = b;
			data[2] = c;
		}

		matrix3(
			T e00, T e10, T e20,
			T e01, T e11, T e21,
			T e02, T e12, T e22
		) {
			data[0] = vector3<T>(e00, e10, e20);
			data[1] = vector3<T>(e01, e11, e21);
			data[2] = vector3<T>(e02, e12, e22);
		}

		template <typename U>
		matrix3(const matrix3<U> &other) {
			data[0] = other.data[0];
			data[1] = other.data[1];
			data[2] = other.data[2];
		}

		static matrix3 random(T lower = 0, T upper = 1) { 
			return matrix3(vector3<T>::random(lower, upper), vector3<T>::random(lower, upper), vector3<T>::random(lower, upper));
		}

		static matrix3 identity() {return matrix3(1);}

		explicit operator T *() {
			return &(data[0].x);
		}

		T * dataPointer() {
			return &(data[0].x);
		}

		const T * dataPointer() const{
			return &(data[0].x);
		}

		const vector3<T> & operator[](size_t i) const {
			assert(i < 3);
			return data[i];
		}

		vector3<T> & operator[](size_t i) {
			assert(i < 3);
			return data[i];
		}

		// stream insertion
		inline friend std::ostream & operator<<(std::ostream &out, const matrix3 &m) {
			const size_t field_width = 10;
			std::ostringstream oss;
			oss << std::setprecision(4);
			oss << '[' << std::setw(field_width) << m[0][0] << ", " << std::setw(field_width) << m[1][0] << ", "
				<< std::setw(field_width) << m[2][0] << ']' << std::endl;
			oss << '[' << std::setw(field_width) << m[0][1] << ", " << std::setw(field_width) << m[1][1] << ", "
				<< std::setw(field_width) << m[2][1] << ']' << std::endl;
			oss << '[' << std::setw(field_width) << m[0][2] << ", " << std::setw(field_width) << m[1][2] << ", "
				<< std::setw(field_width) << m[2][2] << ']' << std::endl;
			return out << oss.str();
		}

		static T det2x2(
			T e00, T e01,
			T e10, T e11
		) {
			return determinant(matrix2<T>(e00, e01, e10, e11));
		}

		// Operator overload - assign
		// 

		// assign
		template <typename U>
		matrix3 & operator=(const matrix3<U> &other) {
			data[0] = other.data[0];
			data[1] = other.data[1];
			data[2] = other.data[2];
			return *this;
		}

		// add-assign
		template <typename U>
		matrix3 & operator+=(const matrix3<U> &rhs) {
			data[0] += rhs[0];
			data[1] += rhs[1];
			data[2] += rhs[2];
			return *this;
		}

		// add-assign
		matrix3 & operator+=(T rhs) {
			data[0] += rhs;
			data[1] += rhs;
			data[2] += rhs;
			return *this;
		}

		// subtract-assign
		template <typename U>
		matrix3 & operator-=(const matrix3<U> &rhs) {
			data[0] -= rhs[0];
			data[1] -= rhs[1];
			data[2] -= rhs[2];
			return *this;
		}

		// subtract-assign
		matrix3 & operator-=(T rhs) {
			data[0] -= rhs;
			data[1] -= rhs;
			data[2] -= rhs;
			return *this;
		}

		// This is a special case where we use matrix product
		// if you want component-wise multiplication see matrixCompMult(matrix3, matrix3)
		// 
		// mulitply-assign
		template <typename U>
		matrix3 & operator*=(const matrix3<U> &rhs) {
			const matrix3 lhs = *this;
			for (int i = 0; i < 3; i++) {
				(*this)[i] = vector3<T>(0);
				for (int j = 0; j < 3; j++) {
					(*this)[i] += lhs[j] * rhs[i][j];
				}
			}
			return *this;
		}

		// mulitply-assign
		matrix3 & operator*=(T rhs) {
			data[0] *= rhs;
			data[1] *= rhs;
			data[2] *= rhs;
			return *this;
		}

		// divide-assign
		template <typename U>
		matrix3 & operator/=(const matrix3<U> &rhs) {
			data[0] /= rhs[0];
			data[1] /= rhs[1];
			data[2] /= rhs[2];
			return *this;
		}

		// divide-assign
		matrix3 & operator/=(T rhs) {
			data[0] /= rhs;
			data[1] /= rhs;
			data[2] /= rhs;
			return *this;
		}

	};

	// Matrix / Matrix Operator Overloads
	// Matrix / Vector Operator Overloads
	//

	// negate
	template <typename T>
	inline matrix3<T> operator-(const matrix3<T> &rhs) {
		matrix3<T> m;
		m[0] = -rhs[0];
		m[1] = -rhs[1];
		m[2] = -rhs[2];
		return m;
	}

	// add
	template <typename T1, typename T2>
	inline auto operator+(const matrix3<T1> &lhs, const matrix3<T2> &rhs) {
		matrix3<std::common_type_t<T1, T2>> m = lhs;
		return m += rhs;
	}

	// subtract
	template <typename T1, typename T2>
	inline auto operator-(const matrix3<T1> &lhs, const matrix3<T2> &rhs) {
		matrix3<std::common_type_t<T1, T2>> m = lhs;
		return m -= rhs;
	}

	// This is a special case where we use matrix product
	// if you want component-wise multiplication see matrixCompMult(matrix3<T>, matrix3<T>)
	// 
	// multiply
	template <typename T1, typename T2>
	inline auto operator*(const matrix3<T1> &lhs, const matrix3<T2> &rhs) {
		matrix3<std::common_type_t<T1, T2>> m = lhs;
		return m *= rhs;
	}

	// Left multiply matrix3<T> m with vector3<T> v
	// 
	// multiply
	template <typename T1, typename T2>
	inline auto operator*(const matrix3<T1> &lhs, const vector3<T2> &rhs) {
		vector3<std::common_type_t<T1, T2>> v(0);
		for (int j = 0; j < 3; j++) {
			v += lhs[j] * rhs[j];
		}
		return v;
	}

	// Right multiply vector3<T> v with matrix3<T> m
	// Equvilent to transpose(m) * v
	// 
	// mulitply-assign
	template <typename T1, typename T2>
	inline auto operator*=(vector3<T1> &lhs, const matrix3<T2> &rhs) {
		return vector3<std::common_type_t<T1, T2>>(dot(lhs, rhs[0]), dot(lhs, rhs[1]), dot(lhs, rhs[2]));
	}

	// multiply
	template <typename T1, typename T2>
	inline auto operator*(const vector3<T1> &lhs, const matrix3<T2> &rhs) {
		vector3<std::common_type_t<T1, T2>> v = lhs;
		return v*=rhs;
	}

	// divide
	template <typename T1, typename T2>
	inline auto operator/(const matrix3<T1> &lhs, const matrix3<T2> &rhs) {
		matrix3<std::common_type_t<T1, T2>> m = lhs;
		return m /= rhs;
	}


	// Vector / Scalar Operator Overloads
	//

	// add right
	template <typename T1, typename T2>
	inline auto operator+(const matrix3<T1> &lhs, T2 rhs) {
		matrix3<std::common_type_t<T1, T2>> m = lhs;
		return m += rhs;
	}

	// add left
	template <typename T1, typename T2>
	inline auto operator+(T1 lhs, const matrix3<T2> &rhs) {
		matrix3<std::common_type_t<T1, T2>> m = rhs;
		return m += lhs;
	}

	// subtract right
	template <typename T1, typename T2>
	inline auto operator-(const matrix3<T1> &lhs, T2 rhs) {
		matrix3<std::common_type_t<T1, T2>> m = lhs;
		return m -= rhs;
	}

	// subtract left
	template <typename T1, typename T2>
	inline auto operator-(T1 lhs, const matrix3<T2> &rhs) {
		vector3<std::common_type_t<T1, T2>> v(lhs);
		matrix3<std::common_type_t<T1, T2>> m(v, v, v);
		return m -= rhs;
	}

	// multiply right
	template <typename T1, typename T2>
	inline auto operator*(const matrix3<T1> &lhs, T2 rhs) {
		matrix3<std::common_type_t<T1, T2>> m = lhs;
		return m *= rhs;
	}

	// multiply left
	template <typename T1, typename T2>
	inline auto operator*(T1 lhs, const matrix3<T2> &rhs) {
		matrix3<std::common_type_t<T1, T2>> m = rhs;
		return m *= lhs;
	}

	// divide right
	template <typename T1, typename T2>
	inline auto operator/(const matrix3<T1> &lhs, T2 rhs) {
		matrix3<std::common_type_t<T1, T2>> m = lhs;
		return m /= rhs;
	}

	// divide left
	template <typename T1, typename T2>
	inline auto operator/(T1 lhs, const matrix3<T2> &rhs) {
		vector3<std::common_type_t<T1, T2>> v(lhs);
		matrix3<std::common_type_t<T1, T2>> m(v, v, v);
		return m /= rhs;
	}

	// Matrix functions
	// 

	// determinant of matrix
	template <typename T>
	inline T determinant(const matrix3<T> &m) {
		T d = 0;
		d += m[0][0] * m[1][1] * m[2][2];
		d += m[0][1] * m[1][2] * m[2][0];
		d += m[0][2] * m[1][0] * m[2][1];
		d -= m[0][0] * m[1][2] * m[2][1];
		d -= m[0][1] * m[1][0] * m[2][2];
		d -= m[0][2] * m[1][1] * m[2][0];
		return d;
	}

	// inverse of matrix (error if not invertible)
	template <typename T>
	inline matrix3<T> inverse(const matrix3<T> &m) {
		matrix3<T> mi;
		// first column of cofactors, can use for determinant
		T c00 = matrix3<T>::det2x2(m[1][1], m[1][2], m[2][1], m[2][2]);
		T c01 = -matrix3<T>::det2x2(m[1][0], m[1][2], m[2][0], m[2][2]);
		T c02 = matrix3<T>::det2x2(m[1][0], m[1][1], m[2][0], m[2][1]);
		// get determinant by expanding about first column
		T invdet = 1 / (m[0][0] * c00 + m[0][1] * c01 + m[0][2] * c02);
		// FIXME proper detect infinite determinant
		assert(!isinf(invdet) && invdet == invdet && invdet != 0);
		// transpose of cofactor matrix * (1 / det)
		mi[0][0] = c00 * invdet;
		mi[1][0] = c01 * invdet;
		mi[2][0] = c02 * invdet;
		mi[0][1] = -matrix3<T>::det2x2(m[0][1], m[0][2], m[2][1], m[2][2]) * invdet;
		mi[1][1] = matrix3<T>::det2x2(m[0][0], m[0][2], m[2][0], m[2][2]) * invdet;
		mi[2][1] = -matrix3<T>::det2x2(m[0][0], m[0][1], m[2][0], m[2][1]) * invdet;
		mi[0][2] = matrix3<T>::det2x2(m[0][1], m[0][2], m[1][1], m[1][2]) * invdet;
		mi[1][2] = -matrix3<T>::det2x2(m[0][0], m[0][2], m[1][0], m[1][2]) * invdet;
		mi[2][2] = matrix3<T>::det2x2(m[0][0], m[0][1], m[1][0], m[1][1]) * invdet;
		return mi;
	}

	// transpose of matrix
	template <typename T>
	inline matrix3<T> transpose(const matrix3<T> &m) {
		return matrix3<T>(
			m[0][0], m[1][0], m[2][0],
			m[0][1], m[1][1], m[2][1],
			m[0][2], m[1][2], m[2][2]
		);
	}

	// component-wise multiplication 
	// see (*) operator overload for matrix product
	template <typename T1, typename T2>
	inline auto matrixCompMult(const matrix3<T1> &lhs, const matrix3<T2> &rhs) {
		matrix3<std::common_type_t<T1, T2>> m = lhs;
		m[0] *= rhs[0];
		m[1] *= rhs[1];
		m[2] *= rhs[2];
		return m;
	}

	template <typename T1, typename T2>
	inline auto outerProduct(const vector3<T1> &lhs, const vector3<T2> &rhs){
		matrix3<std::common_type_t<T1, T2>> m;
		m[0] = lhs * rhs[0];
		m[1] = lhs * rhs[1];
		m[2] = lhs * rhs[2];
		return m;
	}




	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////    _  _        _  _     __  __       _______ _____  _______   __                                                 ////
	////   | || |      | || |   |  \/  |   /\|__   __|  __ \|_   _\ \ / /                                               ////
	////   | || |___  _| || |_  | \  / |  /  \  | |  | |__) | | |  \ V /                                              ////
	////   |__   _\ \/ /__   _| | |\/| | / /\ \ | |  |  _  /  | |   > <                                             ////
	////      | |  >  <   | |   | |  | |/ ____ \| |  | | \ \ _| |_ / . \                                              ////
	////      |_| /_/\_\  |_|   |_|  |_/_/    \_\_|  |_|  \_\_____/_/ \_\                                               ////
	////                                                                                                                  ////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	template <typename T>
	class matrix4 {
	private:
		vector4<T> data[4];

		void clear(T i=0) {
			data[0] = vector4<T>(i, 0, 0, 0);
			data[1] = vector4<T>(0, i, 0, 0);
			data[2] = vector4<T>(0, 0, i, 0);
			data[3] = vector4<T>(0, 0, 0, i);
		}

	public:
		matrix4() {
			clear();
		}

		explicit matrix4(T i) {
			clear(i);
		}

		template <typename U1, typename U2, typename U3, typename U4>
		matrix4(const vector4<U1> &a, const vector4<U2> &b, const vector4<U3> &c, const vector4<U4> &d) {
			data[0] = a;
			data[1] = b;
			data[2] = c;
			data[3] = d;
		}

		matrix4(
			T e00, T e10, T e20, T e30, 
			T e01, T e11, T e21, T e31, 
			T e02, T e12, T e22, T e32, 
			T e03, T e13, T e23, T e33
		) {
			data[0] = vector4<T>(e00, e10, e20, e30);
			data[1] = vector4<T>(e01, e11, e21, e31);
			data[2] = vector4<T>(e02, e12, e22, e32);
			data[3] = vector4<T>(e03, e13, e23, e33);
		}

		template <typename U>
		matrix4(const matrix4<U> &other) {
			data[0] = other.data[0];
			data[1] = other.data[1];
			data[2] = other.data[2];
			data[3] = other.data[3];
		}

		static matrix4 random(T lower = 0, T upper = 1) { 
			return matrix4(vector4<T>::random(lower, upper), vector4<T>::random(lower, upper),
				vector4<T>::random(lower, upper), vector4<T>::random(lower, upper));
		}

		static matrix4 identity() {return matrix4(1);}

		template <typename U1, typename U2, typename U3>
		static matrix4 lookAt(const vector3<U1> &eye,  const vector3<U2> &lookAt, const vector3<U3> &up) {
			vector3<T> vz = normalize(eye-lookAt);
			vector3<T> vx = normalize(cross(up, vz));
			vector3<T> vy = normalize(cross(vz, vx));
			matrix4 m = matrix4(
				vector4<T>(vx, 0),
				vector4<T>(vy, 0),
				vector4<T>(vz, 0),
				vector4<T>(eye, 1));
			return inverse(m);
		}

		static matrix4 lookAt(T ex, T ey, T ez, T lx, T ly, T lz, T ux, T uy, T uz) {
			return matrix4::lookAt(vector3<T>(ex, ey, ez), vector3<T>(lx, ly, lz), vector3<T>(ux, uy, uz));
		}

		// fovy in radians, aspect is w/h
		static matrix4 perspectiveProjection(T fovy, T aspect, T zNear, T zFar) {
			T f = T(1) / (fovy / T(2));

			matrix4 m;
			m[0][0] = f / aspect;
			m[1][1] = f;
			m[2][2] = (zFar + zNear) / (zNear - zFar);
			m[3][2] = (2 * zFar * zNear) / (zNear - zFar);
			m[2][3] = -1;
			return m;
		}

		static matrix4 orthographicProjection(T left, T right, T bottom, T top, T nearVal, T farVal) {
			matrix4 m;
			m[0][0] = T(2) / (right - left);
			m[3][0] = (right + left) / (right - left);
			m[1][1] = T(2) / (top - bottom);
			m[3][1] = (top + bottom) / (top - bottom);
			m[2][2] = -T(2) / (farVal - nearVal);
			m[3][2] = (farVal + nearVal) / (farVal - nearVal);
			m[3][3] = T(1);
			return m;
		}

		static matrix4 shear(int t_dim, int s_dim, T f) {
			matrix4 m;
			m[t_dim][s_dim] = f;
			return m;
		}

		static matrix4 translate(T dx, T dy, T dz) {
			matrix4 m(1);
			m[3][0] = dx;
			m[3][1] = dy;
			m[3][2] = dz;
			return m;
		}

		template <typename U>
		static matrix4 translate(const vector3<U> & d) {
			return matrix4::translate(d.x, d.y, d.z);
		}

		static matrix4 scale(T fx, T fy, T fz) {
			matrix4 m(1);
			m[0][0] = fx;
			m[1][1] = fy;
			m[2][2] = fz;
			return m;
		}
		
		template <typename U>
		static matrix4 scale(const vector3<U> & f) {
			return matrix4::scale(f.x(), f.y(), f.z());
		}
		
		static matrix4 scale(T f) {
			return matrix4::scale(f, f, f);
		}

		static matrix4 rotateX(T angle) {
			matrix4 m(1);
			m[1][1] = std::cos(angle);
			m[2][1] = -std::sin(angle);
			m[1][2] = std::sin(angle);
			m[2][2] = std::cos(angle);
			return m;
		}
		
		static matrix4 rotateY(T angle) {
			matrix4 m(1);
			m[0][0] = std::cos(angle);
			m[2][0] = std::sin(angle);
			m[0][2] = -std::sin(angle);
			m[2][2] = std::cos(angle);
			return m;
		}
		
		static matrix4 rotateZ(T angle) {
			matrix4 m(1);
			m[0][0] = std::cos(angle);
			m[1][0] = -std::sin(angle);
			m[0][1] = std::sin(angle);
			m[1][1] = std::cos(angle);
			return m;
		}

		explicit operator T *() {
			return &(data[0].x);
		}

		T * dataPointer() {
			return &(data[0].x);
		}

		const T * dataPointer() const{
			return &(data[0].x);
		}

		const vector4<T> & operator[](size_t i) const {
			assert(i < 4);
			return data[i];
		}

		vector4<T> & operator[](size_t i) {
			assert(i < 4);
			return data[i];
		}

		// stream insertion
		inline friend std::ostream & operator<<(std::ostream &out, const matrix4 &m) {
			const size_t field_width = 10;
			std::ostringstream oss;
			oss << std::setprecision(4);
			oss << '[' << std::setw(field_width) << m[0][0] << ", " << std::setw(field_width) << m[1][0] << ", "
				<< std::setw(field_width) << m[2][0] << ", " << std::setw(field_width) << m[3][0] << ']' << std::endl;
			oss << '[' << std::setw(field_width) << m[0][1] << ", " << std::setw(field_width) << m[1][1] << ", "
				<< std::setw(field_width) << m[2][1] << ", " << std::setw(field_width) << m[3][1] << ']' << std::endl;
			oss << '[' << std::setw(field_width) << m[0][2] << ", " << std::setw(field_width) << m[1][2] << ", "
				<< std::setw(field_width) << m[2][2] << ", " << std::setw(field_width) << m[3][2] << ']' << std::endl;
			oss << '[' << std::setw(field_width) << m[0][3] << ", " << std::setw(field_width) << m[1][3] << ", "
				<< std::setw(field_width) << m[2][3] << ", " << std::setw(field_width) << m[3][3] << ']';
			return out << oss.str();
		}

		static T det3x3(
			T e00, T e01, T e02,
			T e10, T e11, T e12,
			T e20, T e21, T e22
		) {
			T d = 0;
			d += e00 * e11 * e22;
			d += e01 * e12 * e20;
			d += e02 * e10 * e21;
			d -= e00 * e12 * e21;
			d -= e01 * e10 * e22;
			d -= e02 * e11 * e20;
			return d;
		}

		// Operator overload - assign
		// 

		// assign
		template <typename U>
		matrix4 & operator=(const matrix4<U> &other) {
			data[0] = other.data[0];
			data[1] = other.data[1];
			data[2] = other.data[2];
			data[3] = other.data[3];
			return *this;
		}

		// add-assign
		template <typename U>
		matrix4 & operator+=(const matrix4<U> &rhs) {
			data[0] += rhs[0];
			data[1] += rhs[1];
			data[2] += rhs[2];
			data[3] += rhs[3];
			return *this;
		}

		// add-assign
		matrix4 & operator+=(T rhs) {
			data[0] += rhs;
			data[1] += rhs;
			data[2] += rhs;
			data[3] += rhs;
			return *this;
		}

		// subtract-assign
		template <typename U>
		matrix4 & operator-=(const matrix4<U> &rhs) {
			data[0] -= rhs[0];
			data[1] -= rhs[1];
			data[2] -= rhs[2];
			data[3] -= rhs[3];
			return *this;
		}

		// subtract-assign
		matrix4 & operator-=(T rhs) {
			data[0] -= rhs;
			data[1] -= rhs;
			data[2] -= rhs;
			data[3] -= rhs;
			return *this;
		}

		// This is a special case where we use matrix product
		// if you want component-wise multiplication see matrixCompMult(matrix4, matrix4)
		// 
		// mulitply-assign
		template <typename U>
		matrix4 & operator*=(const matrix4<U> &rhs) {
			const matrix4 lhs = *this;
			for (int i = 0; i < 4; i++) {
				(*this)[i] = vector4<T>(0);
				for (int j = 0; j < 4; j++) {
					(*this)[i] += lhs[j] * rhs[i][j];
				}
			}
			return *this;
		}

		// mulitply-assign
		matrix4 & operator*=(T rhs) {
			data[0] *= rhs;
			data[1] *= rhs;
			data[2] *= rhs;
			data[3] *= rhs;
			return *this;
		}

		// divide-assign
		template <typename U>
		matrix4 & operator/=(const matrix4<U> &rhs) {
			data[0] /= rhs[0];
			data[1] /= rhs[1];
			data[2] /= rhs[2];
			data[3] /= rhs[3];
			return *this;
		}

		// divide-assign
		matrix4 & operator/=(T rhs) {
			data[0] /= rhs;
			data[1] /= rhs;
			data[2] /= rhs;
			data[3] /= rhs;
			return *this;
		}

	};

	// Matrix / Matrix Operator Overloads
	// Matrix / Vector Operator Overloads
	//

	// negate
	template <typename T>
	inline matrix4<T> operator-(const matrix4<T> &rhs) {
		matrix4<T> m;
		m[0] = -rhs[0];
		m[1] = -rhs[1];
		m[2] = -rhs[2];
		m[3] = -rhs[3];
		return m;
	}

	// add
	template <typename T1, typename T2>
	inline auto operator+(const matrix4<T1> &lhs, const matrix4<T2> &rhs) {
		matrix4<std::common_type_t<T1, T2>> m = lhs;
		return m += rhs;
	}

	// subtract
	template <typename T1, typename T2>
	inline auto operator-(const matrix4<T1> &lhs, const matrix4<T2> &rhs) {
		matrix4<std::common_type_t<T1, T2>> m = lhs;
		return m -= rhs;
	}

	// This is a special case where we use matrix product
	// if you want component-wise multiplication see matrixCompMult(matrix4<T>, matrix4<T>)
	// 
	// multiply
	template <typename T1, typename T2>
	inline auto operator*(const matrix4<T1> &lhs, const matrix4<T2> &rhs) {
		matrix4<std::common_type_t<T1, T2>> m = lhs;
		return m *= rhs;
	}

	// Left multiply matrix4<T> m with vector4<T> v
	// 
	// multiply
	template <typename T1, typename T2>
	inline auto operator*(const matrix4<T1> &lhs, const vector4<T2> &rhs) {
		vector4<std::common_type_t<T1, T2>> v(0);
		for (int j = 0; j < 4; j++) {
			v += lhs[j] * rhs[j];
		}
		return v;
	}

	// Right multiply vector4<T> v with matrix4<T> m
	// Equvilent to transpose(m) * v
	// 
	// mulitply-assign
	template <typename T1, typename T2>
	inline auto operator*=(vector4<T1> &lhs, const matrix4<T2> &rhs) {
		return vector4<std::common_type_t<T1, T2>>(dot(lhs, rhs[0]), dot(lhs, rhs[1]), dot(lhs, rhs[2]), dot(lhs, rhs[3]));
	}

	// multiply
	template <typename T1, typename T2>
	inline auto operator*(const vector4<T1> &lhs, const matrix4<T2> &rhs) {
		vector4<std::common_type_t<T1, T2>> v = lhs;
		return v*=rhs;
	}

	// divide
	template <typename T1, typename T2>
	inline auto operator/(const matrix4<T1> &lhs, const matrix4<T2> &rhs) {
		matrix4<std::common_type_t<T1, T2>> m = lhs;
		return m /= rhs;
	}


	// Vector / Scalar Operator Overloads
	//

	// add right
	template <typename T1, typename T2>
	inline auto operator+(const matrix4<T1> &lhs, T2 rhs) {
		matrix4<std::common_type_t<T1, T2>> m = lhs;
		return m += rhs;
	}

	// add left
	template <typename T1, typename T2>
	inline auto operator+(T1 lhs, const matrix4<T2> &rhs) {
		matrix4<std::common_type_t<T1, T2>> m = rhs;
		return m += lhs;
	}

	// subtract right
	template <typename T1, typename T2>
	inline auto operator-(const matrix4<T1> &lhs, T2 rhs) {
		matrix4<std::common_type_t<T1, T2>> m = lhs;
		return m -= rhs;
	}

	// subtract left
	template <typename T1, typename T2>
	inline auto operator-(T1 lhs, const matrix4<T2> &rhs) {
		vector4<std::common_type_t<T1, T2>> v(lhs);
		matrix4<std::common_type_t<T1, T2>> m(v, v, v, v);
		return m -= rhs;
	}

	// multiply right
	template <typename T1, typename T2>
	inline auto operator*(const matrix4<T1> &lhs, T2 rhs) {
		matrix4<std::common_type_t<T1, T2>> m = lhs;
		return m *= rhs;
	}

	// multiply left
	template <typename T1, typename T2>
	inline auto operator*(T1 lhs, const matrix4<T2> &rhs) {
		matrix4<std::common_type_t<T1, T2>> m = rhs;
		return m *= lhs;
	}

	// divide right
	template <typename T1, typename T2>
	inline auto operator/(const matrix4<T1> &lhs, T2 rhs) {
		matrix4<std::common_type_t<T1, T2>> m = lhs;
		return m /= rhs;
	}

	// divide left
	template <typename T1, typename T2>
	inline auto operator/(T1 lhs, const matrix4<T2> &rhs) {
		vector4<std::common_type_t<T1, T2>> v(lhs);
		matrix4<std::common_type_t<T1, T2>> m(v, v, v, v);
		return m /= rhs;
	}

	// Matrix functions
	// 

	// determinant of matrix
	template <typename T>
	inline T determinant(const matrix4<T> &m) {
		T d = 0;
		// expand about first column
		d += m[0][0] * matrix4<T>::det3x3(m[1][1], m[1][2], m[1][3], m[2][1], m[2][2], m[2][3], m[3][1], m[3][2], m[3][3]);
		d -= m[0][1] * matrix4<T>::det3x3(m[1][0], m[1][2], m[1][3], m[2][0], m[2][2], m[2][3], m[3][0], m[3][2], m[3][3]);
		d += m[0][2] * matrix4<T>::det3x3(m[1][0], m[1][1], m[1][3], m[2][0], m[2][1], m[2][3], m[3][0], m[3][1], m[3][3]);
		d -= m[0][3] * matrix4<T>::det3x3(m[1][0], m[1][1], m[1][2], m[2][0], m[2][1], m[2][2], m[3][0], m[3][1], m[3][2]);
		return d;
	}

	// inverse of matrix (error if not invertible)
	template <typename T>
	inline matrix4<T> inverse(const matrix4<T> &m) {
		matrix4<T> mi;
		// first column of cofactors, can use for determinant
		T c00 =  matrix4<T>::det3x3(m[1][1], m[1][2], m[1][3], m[2][1], m[2][2], m[2][3], m[3][1], m[3][2], m[3][3]);
		T c01 = -matrix4<T>::det3x3(m[1][0], m[1][2], m[1][3], m[2][0], m[2][2], m[2][3], m[3][0], m[3][2], m[3][3]);
		T c02 =  matrix4<T>::det3x3(m[1][0], m[1][1], m[1][3], m[2][0], m[2][1], m[2][3], m[3][0], m[3][1], m[3][3]);
		T c03 = -matrix4<T>::det3x3(m[1][0], m[1][1], m[1][2], m[2][0], m[2][1], m[2][2], m[3][0], m[3][1], m[3][2]);
		// get determinant by expanding about first column
		T invdet = 1 / (m[0][0] * c00 + m[0][1] * c01 + m[0][2] * c02 + m[0][3] * c03);
		// FIXME proper detect infinite determinant
		assert(!isinf(invdet) && invdet == invdet && invdet != 0);
		// transpose of cofactor matrix * (1 / det)
		mi[0][0] = c00 * invdet;
		mi[1][0] = c01 * invdet;
		mi[2][0] = c02 * invdet;
		mi[3][0] = c03 * invdet;
		mi[0][1] = -matrix4<T>::det3x3(m[0][1], m[0][2], m[0][3], m[2][1], m[2][2], m[2][3], m[3][1], m[3][2], m[3][3]) * invdet;
		mi[1][1] =  matrix4<T>::det3x3(m[0][0], m[0][2], m[0][3], m[2][0], m[2][2], m[2][3], m[3][0], m[3][2], m[3][3]) * invdet;
		mi[2][1] = -matrix4<T>::det3x3(m[0][0], m[0][1], m[0][3], m[2][0], m[2][1], m[2][3], m[3][0], m[3][1], m[3][3]) * invdet;
		mi[3][1] =  matrix4<T>::det3x3(m[0][0], m[0][1], m[0][2], m[2][0], m[2][1], m[2][2], m[3][0], m[3][1], m[3][2]) * invdet;
		mi[0][2] =  matrix4<T>::det3x3(m[0][1], m[0][2], m[0][3], m[1][1], m[1][2], m[1][3], m[3][1], m[3][2], m[3][3]) * invdet;
		mi[1][2] = -matrix4<T>::det3x3(m[0][0], m[0][2], m[0][3], m[1][0], m[1][2], m[1][3], m[3][0], m[3][2], m[3][3]) * invdet;
		mi[2][2] =  matrix4<T>::det3x3(m[0][0], m[0][1], m[0][3], m[1][0], m[1][1], m[1][3], m[3][0], m[3][1], m[3][3]) * invdet;
		mi[3][2] = -matrix4<T>::det3x3(m[0][0], m[0][1], m[0][2], m[1][0], m[1][1], m[1][2], m[3][0], m[3][1], m[3][2]) * invdet;
		mi[0][3] = -matrix4<T>::det3x3(m[0][1], m[0][2], m[0][3], m[1][1], m[1][2], m[1][3], m[2][1], m[2][2], m[2][3]) * invdet;
		mi[1][3] =  matrix4<T>::det3x3(m[0][0], m[0][2], m[0][3], m[1][0], m[1][2], m[1][3], m[2][0], m[2][2], m[2][3]) * invdet;
		mi[2][3] = -matrix4<T>::det3x3(m[0][0], m[0][1], m[0][3], m[1][0], m[1][1], m[1][3], m[2][0], m[2][1], m[2][3]) * invdet;
		mi[3][3] =  matrix4<T>::det3x3(m[0][0], m[0][1], m[0][2], m[1][0], m[1][1], m[1][2], m[2][0], m[2][1], m[2][2]) * invdet;
		return mi;
	}

	// transpose of matrix
	template <typename T>
	inline matrix4<T> transpose(const matrix4<T> &m) {
		return matrix4<T>(
			m[0][0], m[1][0], m[2][0], m[3][0],
			m[0][1], m[1][1], m[2][1], m[3][1],
			m[0][2], m[1][2], m[2][2], m[3][2],
			m[0][3], m[1][3], m[2][3], m[3][3]
		);
	}

	// component-wise multiplication 
	// see (*) operator overload for matrix product
	template <typename T1, typename T2>
	inline auto matrixCompMult(const matrix4<T1> &lhs, const matrix4<T2> &rhs) {
		matrix4<std::common_type_t<T1, T2>> m = lhs;
		m[0] *= rhs[0];
		m[1] *= rhs[1];
		m[2] *= rhs[2];
		m[3] *= rhs[3];
		return m;
	}

	template <typename T1, typename T2>
	inline auto outerProduct(const vector4<T1> &lhs, const vector4<T2> &rhs){
		matrix4<std::common_type_t<T1, T2>> m;
		m[0] = lhs * rhs[0];
		m[1] = lhs * rhs[1];
		m[2] = lhs * rhs[2];
		m[3] = lhs * rhs[3];
		return m;
	}

}

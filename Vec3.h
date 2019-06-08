#ifndef VEC3_H
#define VEC3_H
namespace Cu{
    template <typename T>
	class Vec3{
	public:
		Vec3();
		Vec3(T x, T y, T z);
		Vec3(const Vec3& in);
		const Vec3& operator=(const Vec3& in);
		const Vec3& operator+=(const Vec3& in);
		T x, y, z;
	};

    template <typename T>
	Vec3<T> operator+(const Vec3<T>& in1, const Vec3<T>& in2);
	template <typename T>
	Vec3<T> operator-(const Vec3<T>& in1, const Vec3<T>& in2);
	template <typename T>
	Vec3<T> operator*(const float& in1, const Vec3<T>& in2);
	template <typename T>
	Vec3<T> operator*(const Vec3<T>& in1, const float& in2);
	template <typename T>
	Vec3<T> operator/(const Vec3<T>& in1, const float& in2);
    
    #include "Vec3.hpp"
}

#endif
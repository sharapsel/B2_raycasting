#include <fstream>
#include <cmath>
#include <png.h>
#include <time.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <cstring>
using namespace std;

struct Vec3 
{
	double x, y, z;
	Vec3(){}
	Vec3(double x, double y, double z) : x(x), y(y), z(z) {}
	Vec3 operator + (Vec3& v) 
	{ 
		return Vec3(x + v.x, y + v.y, z + v.z); 
	}
	Vec3 operator - (Vec3& v) 
	{ 
		return Vec3(x - v.x, y - v.y, z - v.z); 
	}
	Vec3 operator * (double d)
	{ 
		return Vec3(x * d, y * d, z * d); 
	}
	Vec3 operator / (double d) 
	{
		return Vec3(x / d, y / d, z / d);
	}
	Vec3 normalize() 
	{
		double n = sqrt(x * x + y * y + z * z);
		return Vec3(x / n, y / n, z / n);
	}
};

double dot(Vec3& a, Vec3& b)
{
	return (a.x * b.x + a.y * b.y + a.z * b.z);
}

struct Ray
{
	Vec3 o, d;
	Ray(Vec3& o, Vec3& d) : o(o), d(d) {}
};

struct Sphere
{
	Vec3 c;
	double r;
	Vec3 col;
	Sphere() {}
	Sphere(Vec3& c, double r, Vec3& col) : c(c), r(r), col(col) {}
	Vec3 getNormal(Vec3& pi)
	{
		return (pi - c) / r;
	}
	bool intersect(Ray& ray, double &t)
	{
		Vec3 o = ray.o;
		Vec3 d = ray.d;
		Vec3 oc = o - c;
		double b = 2 * dot(oc, d);
		double c = dot(oc, oc) - r*r;
		double disc = b * b - 4 * c;
		if (disc < 1e-4)
			return false;
		disc = sqrt(disc);
		double t0 = -b - disc;
		double t1 = -b + disc;
		t = (t0 < t1) ? t0 : t1;
		return true;
	}
	bool sphereIntersect(Sphere& s2)
	{
		Vec3 delta = s2.c - c;
		double dist = sqrt(dot(delta, delta));
		return dist <= r + s2.r;
	}
};

struct Pixel
{
	uint8_t r;
	uint8_t g;
	uint8_t b;
	Pixel() : r(0), g(0), b(0) {}
	Pixel(uint8_t r, uint8_t g, uint8_t b) : r(r), g(g), b(b) {}
};

void vecToUINT8(Vec3& col, Pixel& p)
{
	p.r = (uint8_t)((col.x > 255) ? 255 : (col.x < 0) ? 0 : col.x);
	p.g = (uint8_t)((col.y > 255) ? 255 : (col.y < 0) ? 0 : col.y);
	p.b = (uint8_t)((col.z > 255) ? 255 : (col.z < 0) ? 0 : col.z);
}

struct Bitmap
{
	int h, w;
	Pixel* pixels;
	Bitmap(int h, int w) : h(h), w(w)
	{
		pixels = new Pixel[h * w];
	}
	void Clear()
	{
		for (int i = 0; i < h*w; i++)
			pixels[i] = Pixel();
	}
	~Bitmap()
	{
		delete(pixels);
	}
};

int saveBitmapToPNG(Bitmap& b, const char *path)
{
	/* create file */
	FILE * fp;
	png_structp png_ptr = NULL;
	png_infop info_ptr = NULL;
	int pixel_size = 3;
	int depth = 8;
	int code = 0;
	fp = fopen(path, "wb");
	if (!fp) 
	{
		code = 1;
		goto finalise;
	}
	/* initialize stuff */
	png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
	if (png_ptr == NULL) 
	{
		code = 2;
		goto finalise;
	}
	info_ptr = png_create_info_struct(png_ptr);
	if (info_ptr == NULL) 
	{
		code = 3;
		goto finalise;
	}
	if (setjmp(png_jmpbuf(png_ptr))) 
	{
		code = 4;
		goto finalise;
	}
	png_init_io(png_ptr, fp);
	/* write header */
	if (setjmp(png_jmpbuf(png_ptr)))
	{
		code = 5;
		goto finalise;
	}
	png_set_IHDR(png_ptr, info_ptr, b.w, b.h, depth, PNG_COLOR_TYPE_RGB,
		PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
	png_write_info(png_ptr, info_ptr);
	/* write bytes */
	if (setjmp(png_jmpbuf(png_ptr)))
	{
		code = 6;
		goto finalise;
	}
	png_bytep row = (png_bytep)malloc(pixel_size * b.w * sizeof(png_byte));
	int rev, revCoord;
	for (int y = 0; y < b.h; y++)
	{
		for (int x = 0; x < b.w; x++)
		{
			png_byte* ptr = &(row[x * 3]);
			rev = b.h - y - 1;
			revCoord = rev * b.w + x;
			ptr[0] = b.pixels[revCoord].r;
			ptr[1] = b.pixels[revCoord].g;
			ptr[2] = b.pixels[revCoord].b;
		}
		png_write_row(png_ptr, row);
	}
	/* end write */
	if (setjmp(png_jmpbuf(png_ptr)))
	{
		code = 7;
		goto finalise;
	}
	png_write_end(png_ptr, NULL);
finalise:
	if (fp != NULL) fclose(fp);
	if (info_ptr != NULL) png_free_data(png_ptr, info_ptr, PNG_FREE_ALL, -1);
	if (png_ptr != NULL) png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
	if (row != NULL) delete(row);
	return code;
}

void generateSpheres(Sphere **spheres, int h, int w, int sphereCount)
{
	srand(123456);
	*spheres = new Sphere[sphereCount];
	for (int i = 0; i < sphereCount; i++)
	{
		Sphere s(
			Vec3(rand() % w, rand() % h, rand() % h / 2), 
			h / 20 + (rand() % h / 10), 
			Vec3(rand() % 256, rand() % 256, rand() % 256));
		bool isOk = true;
		for (int j = 0; j < i; j++)
			if ((*spheres)[j].sphereIntersect(s))
			{
				isOk = false;
				break;
			}
		if (isOk)
			(*spheres)[i] = s;
		else i--;
	}
}

bool intersectArray(Sphere * spheres, Ray& ray, double &t, int &id, int spheresCount)
{
	bool intersects = false;
	double d = 0;
	t = 1e10;
	for (int i = 0; i < spheresCount; i++)
	{
		intersects = spheres[i].intersect(ray, d);
		if (intersects && d < t)
		{
			t = d;
			id = i;
		}
	}
	return t < 1e10;
}

void raytraceHost(Bitmap &b, Sphere *light, Sphere *spheres, int spheresCount)
{
	Vec3 black(0, 0, 0);
#pragma omp parallel for
	for (int y = 0; y < b.h; ++y)
	{
		for (int x = 0; x < b.w; ++x)
		{
			Vec3 pix_col(black);
			double t;
			int id;
			Ray ray(Vec3(x, y, -100), Vec3(0, 0, 1));
			bool intersects = intersectArray(spheres, ray, t, id, spheresCount);
			if (intersects)
			{
				Vec3 crossP = ray.o + ray.d*t;
				Vec3 lightDir = (*light).c;
				Vec3 normal = spheres[id].getNormal(crossP);
				double lambertian = dot(lightDir.normalize(), normal.normalize());
				pix_col = (spheres[id].col * lambertian) * 0.7;
				vecToUINT8(pix_col, b.pixels[y * b.w + x]);
			}
		}
	}
}

int main()
{
	clock_t tStart = clock();
	int h[] = { 720, 900, 1080, 1440 };
	int w[] = { 1280, 1600, 1920, 2560 };
	int sphereCount[] = { 100, 200, 300, 400 };
	Sphere lights[] = {
		Sphere(Vec3(200, 500, 150), 1, Vec3(255, 255, 255)),
		Sphere(Vec3(500, 600, 150), 1, Vec3(255, 255, 255)),
		Sphere(Vec3(1300, 200, 200), 1, Vec3(255, 255, 255)),
		Sphere(Vec3(2500, 1400, -100), 1, Vec3(255, 255, 255))
	};
	for (int i = 0; i < 4; i++)
	{
		Sphere *spheres = NULL;
		generateSpheres(&spheres, h[i], w[i], sphereCount[i]);
		Bitmap b(h[i], w[i]);

		tStart = clock();
		raytraceHost(b, &lights[i], spheres, sphereCount[i]);
		cout << "HOST - time elapsed: " << setprecision(5) << (double)(clock() - tStart) / CLOCKS_PER_SEC << " s" << endl;
		string hS = "file_host_";
		hS.operator+=(to_string(i)).operator+=(".png");
		saveBitmapToPNG(b, hS.c_str());
		delete(spheres);
	}
	getchar();
	return 0;
}
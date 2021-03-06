#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <thread>
#include <atomic>
#include <random>

struct RNG
{
    RNG() : gen(rd()), distribution(std::uniform_real_distribution<double>(0.0, 1.0)) { }

    double Next()
    {
        return distribution(gen);
    }

private:
    std::random_device rd;
    std::mt19937 gen;
    std::uniform_real_distribution<double> distribution;
};

struct Vec
{
    double x, y, z;
    Vec(double x_ = 0, double y_ = 0, double z_ = 0) { x = x_; y = y_; z = z_; }
    Vec operator+(const Vec &b) const { return {x + b.x,y + b.y,z + b.z}; }
    Vec operator-(const Vec &b) const { return {x - b.x,y - b.y,z - b.z}; }
    Vec operator-() const { return { -x, -y, -z}; }
    Vec operator*(double b) const { return {x * b,y * b,z * b}; }
    Vec Mult(const Vec &b) const { return {x * b.x, y * b.y, z * b.z}; }
    Vec& Norm() { return *this = *this * (1.0 / std::sqrt(x * x + y * y + z * z)); }
    double Dot(const Vec &b) const { return x * b.x + y * b.y + z * b.z; } // cross:
    Vec operator%(Vec& b) const { return {y * b.z - z * b.y,z * b.x - x * b.z,x * b.y - y * b.x }; }
};

struct Ray { Vec o, d; Ray(Vec o_, Vec d_) : o(o_), d(d_) { } };

struct Camera
{
    Vec eye, u, v, n;
    double w, h;
    double res_x, res_y;

    Camera(Vec o_, Vec dir_, Vec up_, double angle_, int res_x_, int res_y_) : eye(o_), n(dir_), res_x(res_x_), res_y(res_y_)
    {
        Vec nn = -n;
        u = up_ % nn;
        u.Norm();
        v = n % u;
        h = 2.0 * std::tan((M_PI * angle_ / 180.0) / 2.0);
        w = (res_x / res_y) * h;
    }

    Ray GenerateRay(double x, double y) const
    {
        Vec ps;
        ps.x = w * (x / float(res_x) - 0.5);
        ps.y = h * (y / float(res_y) - 0.5);
        ps.z = 1;

        Vec vX = u * ps.x;
        Vec vY = v * ps.y;
        Vec vZ = n * ps.z;
        Vec ray_dir = (vX + vY + vZ).Norm();
        return { eye, ray_dir };
    }
};

enum Interaction { DIFF, SPEC, REFR };

struct Sphere
{
    double rad; // radius
    Vec p, e, c; // position, emission, color
    Interaction interaction;

    Sphere(double rad_, Vec p_, Vec e_, Vec c_, Interaction inter_) :
        rad(rad_), p(p_), e(e_), c(c_), interaction(inter_) { }

    double Intersect(Ray const& r) const
    {
        Vec op = p - r.o; // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
        double t, eps = 1e-4, b = op.Dot(r.d), det = b * b - op.Dot(op) + rad * rad;
        if (det < 0) return 0;
        else det = std::sqrt(det);
        return (t = b - det) > eps ? t : ((t = b + det) > eps ? t : 0);
    }
};

Sphere spheres[] = { //Scene: radius, position, emission, color, material
  Sphere(1e5, Vec( 1e5+1,40.8,81.6), Vec(),Vec(.15,.803,.17),DIFF),//Left
  Sphere(1e5, Vec(-1e5+99,40.8,81.6),Vec(),Vec(.8,.15,.15),DIFF),//Right
  Sphere(1e5, Vec(50,40.8, 1e5),     Vec(),Vec(.8,.8,.8),DIFF),//Back
  Sphere(1e5, Vec(50,40.8,-1e5 + 170), Vec(),Vec(),           DIFF),//Front
  Sphere(1e5, Vec(50, 1e5, 81.6),    Vec(),Vec(.75,.75,.75),DIFF),//Bottom
  Sphere(1e5, Vec(50,-1e5+81.6,81.6),Vec(),Vec(.75,.75,.75),DIFF),//Top
  Sphere(16.5,Vec(27,16.5,47),       Vec(),Vec(1,1,1)*.999, SPEC),//Mirror
  Sphere(16.5,Vec(73,16.5,78),       Vec(),Vec(1,1,1)*.999, REFR),//Glass
  Sphere(600, Vec(50,681.6-.27,81.6),Vec(0.7,0.4,0.2) * 22,  Vec(), DIFF) //Lite
};

inline double clamp(double x) { return x < 0 ? 0 : x > 1 ? 1 : x; }
inline int toInt(double x) { return int(std::pow(clamp(x), 1 / 2.2) * 255 + .5); }

inline Vec Reflect(Vec const& d, Vec const& n)
{
    return d - n * 2 * n.Dot(d);
}

Vec SampleHemisphere(RNG& rng)
{
    double r1 = rng.Next();
    double r2 = rng.Next();

    double phi = 2 * M_PI * r1;
    double theta = std::sqrt(1 - r2);
    return Vec(std::cos(phi) * theta, std::sin(phi) * theta, std::sqrt(r2));
}

inline bool Intersect(const Ray &r, double &t, int &id)
{
    double d, inf = t = 1e20;
    for (int i = 9; i--;)
        if ((d = spheres[i].Intersect(r)) > 0 && d < t)
        {
            t = d;
            id = i;
        }

    return t < inf;
}

Vec Radiance(const Ray &r, int depth, RNG& rng)
{
    double t; // distance to intersection
    int id = 0; // id of intersected object
    if (!Intersect(r, t, id))
        return { 0, 0, 0 }; // if miss, return black

    const Sphere &obj = spheres[id]; // the hit object
    Vec x = r.o + r.d * t;
    Vec n = (x - obj.p).Norm();
    Vec nl = n.Dot(r.d) < 0 ? n : n * -1;
    Vec f = obj.c;
    // russian roulette
    double p = f.x > f.y && f.x > f.z ? f.x : f.y > f.z ? f.y : f.z;
    if (++depth > 5)
    {
        if (rng.Next() < p)
            f = f * (1 / p);
        else
            return obj.e;
    }

    if (obj.interaction == DIFF)
    {
        Vec w = nl, u = ((std::fabs(w.x) > .1 ? Vec(0, 1) : Vec(1)) % w).Norm(), v = w % u;
        Vec local_dir = SampleHemisphere(rng);
        Vec d = (u * local_dir.x + v * local_dir.y + w * local_dir.z).Norm();
        return obj.e + f.Mult(Radiance(Ray(x, d), depth, rng));
    } else if (obj.interaction == SPEC)
        return obj.e + f.Mult(Radiance(Ray(x, Reflect(r.d, n)), depth, rng));

    // Fresnel equation for Refraction
    Ray reflRay(x, Reflect(r.d, n));
    bool into = n.Dot(nl) > 0.0f;    // Ray from outside going in?
    double nc = 1, nt = 1.5, nnt = into ? nc / nt : nt / nc;
    double ddn = r.d.Dot(nl), cos2t;

    // (1 - ddn * ddn) = cos^2
    if ((cos2t = 1 - nnt * nnt * (1 - ddn * ddn)) < 0) // Total internal reflection
        return obj.e + f.Mult(Radiance(reflRay, depth, rng));

    Vec tdir = (r.d * nnt - n * ((into ? 1 : -1) * (ddn * nnt + std::sqrt(cos2t)))).Norm();
    double a = nt - nc, b = nt + nc, R0 = a * a / (b * b), c = 1 - (into ? - ddn : tdir.Dot(n));
    double Re = R0 + (1 - R0) * c * c * c * c * c, Tr = 1 - Re, P = .25 + .5 * Re, RP = Re / P, TP = Tr / (1 - P);
    return obj.e + f.Mult(depth > 2 ? (rng.Next() < P ? // Russian roulette
                          Radiance(reflRay, depth, rng) * RP : Radiance(Ray(x, tdir), depth, rng) * TP) :
                          Radiance(reflRay, depth, rng) * Re + Radiance(Ray(x, tdir), depth, rng) * Tr);
}

int main(int argc, char** argv)
{
    int w = 512, h = 512, samples = argc == 2 ? std::stoi(argv[1]) : 1; // # samples
    Vec* c = new Vec[w * h];
    Camera cam(Vec(50,52,295.6), Vec(0, -0.042612, -1).Norm(), Vec(0, 1, 0), 29, 512, 512);
    std::vector<std::thread> workers; std::atomic<int> work_done;
    unsigned worker_count = std::thread::hardware_concurrency();
    unsigned worker_modulo = h % worker_count;
    unsigned worker_height = h / worker_count;
    unsigned worker_final_inc = worker_modulo ? h - worker_height * worker_count : 0;
    unsigned min = 0; int total_paths = w * h;

    for (unsigned short i = 0; i < worker_count; i++, min += worker_height)
    {
        auto max = i == worker_count - 1 ? min + worker_height + worker_final_inc : min + worker_height;
        workers.emplace_back([min, max, w, samples, cam, c, &work_done]() {
            unsigned pathCount = w * (max - min);
            RNG rng;
            for (unsigned pathId = 0; pathId < pathCount; pathId++)
            {
                unsigned x = pathId % w;
                unsigned y = min + pathId / w;
                Vec r(0, 0, 0);
                for (int s = 0; s < samples; s++)
                {
                    double r1 = rng.Next() - 0.5;
                    double r2 = rng.Next() - 0.5;
                    Ray ray = cam.GenerateRay(x + r1, y + r2);
                    ray.o = ray.o + ray.d * 140; // hack camera origin
                    r = r + Radiance(ray, 0, rng) * (1.0 / samples);
                }

                auto i = static_cast<unsigned>(x + y * w);
                c[i] = Vec(clamp(r.x), clamp(r.y), clamp(r.z));
                work_done++;
            }
        });
    }

    while (work_done < total_paths)
    {
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
        fprintf(stdout, "Rendering: %5.2f %%\r", 100.0 * work_done / double(total_paths));
        fflush(stdout);
    }

    for (auto& item : workers)
        item.join();

    std::cout << std::endl; std::ofstream str("image.ppm");
    str << "P3\n" << w << " " << h << std::endl << 255 << std::endl;
    for (int i = 0;  i< w * h; i++)
        str << toInt(c[i].x) << " " << toInt(c[i].y) << " " << toInt(c[i].z) << " ";
    str.close();
    delete [] c;
    return 0;
}

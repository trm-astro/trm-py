cdef extern from "trm/vec3.h" namespace "Subs":
    cdef cppclass Vec3:
        Vec3()
        Vec3(double, double, double)
        Vec3(double*)
        double& x()
        double& y()
        double& z()
        void set(double, double, double)
        void set(double*)
        void get(double*)
        double length()
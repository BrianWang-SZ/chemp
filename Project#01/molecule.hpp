#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/Core"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;

class Atom {
public:
    int zval;
    double x;
    double y;
    double z;
};

class Molecule {
public:
    int natom;
    Atom *atoms;
    double **length;
    double ***angle;
    
    Molecule(const char*);
    ~Molecule();

    void coord(const char*);
    void bond_length();
    double*** unitvec();
    void bond_angle(double ***);
    void oop_angle(double ***);
    void tor_angle(double ***);
    void massctr(double*);
    void moi(double*);
    void rotation();
    void cleanup(double***);
    
    void print_coord();
    void print_len();
    void print_angle();
    void print_vec(double ***vec);
    void print_moi(Matrix, Matrix);

private:
    void oop_helper(double ***vec, int i, int j, int k, int l);
};

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
    int size;
    Atom *atoms;
    
    Molecule(const char *, const char *);
    void coord(const char *);
    double ** hessian(const char *);
    void weight(double **);
    void compute(double **);
    void cleanup(double **);
    void print_matrix(double **);
};

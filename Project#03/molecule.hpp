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
    int norb = 0, natom;
    Atom *atoms;
    double enuc;
    int *ioff;
    double **s, **t, **v, *eri, **ham;
    double **mux, **muy, **muz;
    
    void coord(const char*);
    void read_one_electron(const char *dir);
    void read_two_electron(const char *dir);
    double** readMatrix(const char *path);
    Molecule(const char *dir);
    //~Molecule();
    void print_matrix(double **mat);
    void print_matrix(Matrix mat);
    void hamiltonian();
    int calc_norb(const char *path);
    void compute_HF();
    void updateFock(Matrix &F, Matrix H, Matrix D);
    void updateDensity(Matrix &D, Matrix isqrt_S, Matrix F, Matrix C);
    void lookupTable();
    void cleanNan(Matrix &M);
    double calc_energy(Matrix D, Matrix H, Matrix F);
    double calc_rms(Matrix D, Matrix new_D);
    void mobasis(Matrix C, Matrix F);
    void read_dipole(const char *dir);
    void compute_dipole(Matrix D);
};

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
    int norb = 0, natom;
    Atom *atoms;
    double enuc;
    int *ioff;
    double **s, **t, **v, *eri, **ham;
    double **mux, **muy, **muz;
    
    void read_coord(const char*);
    void read_one_electron(const char *dir);
    void read_two_electron(const char *dir);
    double** readMatrix(const char *path);
    Molecule(const char *dir);
    ~Molecule();
    void print_matrix(double **mat, int size);
    void print_matrix(Matrix mat);
    void hamiltonian();
    int calc_norb(const char *path);
    void compute();
    double compute_hf(Matrix isqrt_S, Matrix &C, Matrix &D, Matrix &F);
    void initialize(int toprint, Matrix &isqrt_S, Matrix &C, Matrix &D, Matrix &F);
    void updateFock(Matrix &F, Matrix D);
    void updateDensity(Matrix &D, Matrix isqrt_S, Matrix &C, Matrix F);
    double compute_diis(Matrix D, Matrix F, Matrix isqrt_S);
    void lookupTable();
    double calc_hf_energy(Matrix D, Matrix F);
    double calc_ccsd_energy(Matrix Fs, double ****mospin, double ****t_ijab, double **t_ia);
    double calc_rms(Matrix D, Matrix new_D);
    double calc_rms(double **mat, double **newmat, int size);
    void mobasis(Matrix C, Matrix F);
    void read_dipole(const char *dir);
    void compute_dipole(Matrix D);
    double**** create4dmat(int size);
    void free4dmat(double ****mat, int size);
    double compute_mp2(Matrix C, Matrix F, Matrix isqrt_S);
    void massctr(double *ctr);
    double**** spatial_to_spin(double *moeri);
    double* spatial_atom(double *eri, Matrix C);
    int calc_nomo();
    double compute_ccsd(Matrix C, Matrix F, Matrix isqrt_S);
    void free2dvec(double **vec, int size);
    double** create2dvec(int size);
    double tau(double ****t_ijab, double **t_ia, int i, int j, int a, int b);
    double taut(double ****t_ijab, double **t_ia, int i, int j, int a, int b);
    void update_interm(double ****mospin, double ****t_ijab, double **t_ia, 
                       Matrix Fs, Matrix &Fae, Matrix &Fmi, Matrix &Fme,
                       double ****Wmnij, double ****Wabef, double ****Wmbej);
    void updateT(double ****mospin, double ****t_ijab, double **t_ia, Matrix Fs, Matrix Fae, 
                       Matrix Fmi, Matrix Fme, double ****Wmnij, double ****Wabef, double ****Wmbej,
                       double ****D_ijab, double **D_ia);
};

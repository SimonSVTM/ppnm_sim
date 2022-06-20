using System;
using System.Diagnostics;

namespace homework
{
    public partial class matrix
    {
        private delegate T unitary_grid_operation<T>(int x);
        private delegate T binary_grid_operation<T>(int x, int y);
        public readonly int size1, size2;
        private double[] data; // keep matrix elements in one big array
        public matrix(int m, int n)
        {// constructor

            size1 = m; size2 = n;
            data = new double[size1 * size2];
        }

        public void setData(double[] data)
        {
            Debug.Assert(data.Length == size1 * size2);
            this.data = data;
        } 

        new public string ToString()
        {
            String res = "[";
            for (int i = 0; i < size1; i++)
            {
                res += "[";
                for (int j = 0; j < size2; j++)
                {
                    res += this[i, j];
                    if (j < size2 - 1)
                    {
                        res += ", ";
                    }
                }
                if (i == size1 - 1) res += "]]";
                else res += "],";
                res += "\n";

            }
            
            return res;
        }
        public double this[int i, int j]
        { // indexer
            get => data[i + j * size1];      // getter: for example, aij=A[i,j];
            set => data[i + j * size1] = value;// setter: for example, A[i,j]=5;
        }

        private static void mapto(matrix A, binary_grid_operation<double> b)
        {
            for (int i = 0; i < A.size1; i++) for (int j = 0; j < A.size2; j++) A[i, j] = b(i, j);
        }


        private static bool check_condition(int m, int n, binary_grid_operation<bool> b)
        {
            for (int i = 0; i < m; i++) {
                for (int j = 0; j < n; j++)
                {
                    bool cond = b(i, j);
                    if (!cond) return false;
                }

            }
            return true;

        }

        

        public matrix resize(int m, int n)
        {
            matrix this_new = new matrix(m, n);
            mapto(this_new, (i, j) => { if (i < size1 && j < size2) return this[i, j]; else return 0;});
            return this_new;
        }

        public bool isUpperTriangulary()
        {

            return check_condition(this.size1, this.size2, (i, j) => { if (i > j) return approxEqaul(this[i, j], 0); else return true; });
        }

        public bool isEqualGlobal(matrix B)
        {
            Debug.Assert(this.size1 == B.size1 && this.size2 == B.size2, "Dimension mismatch");
            return check_condition(this.size1, this.size2, (i, j) => approxEqaul(this[i, j], B[i, j]));
        }

        private bool approxEqaul(double x, double y)
        {
            return Math.Abs(x - y) < 0.01;
        }

        public matrix getCopy()
        {
            matrix B = new matrix(this.size1, this.size2);
            mapto(B, (i, j) => this[i, j]);
            return B;
        }

        public static matrix getrandomsymmetric(int m, int n, Random rand)
        {
            matrix A = new matrix(m, n);

            int scale = 10000;
            mapto(A, (i, j) =>
            {
                if (j >= i) return rand.NextDouble() * scale;
                return A[j, i];
            });
            return A;
        }

        public matrix transpose()
        {
            matrix At = new matrix(this.size2, this.size1);
            mapto(At, (i, j) => this[j, i]);
            return At;
        }

        public matrix times(matrix B)
        {
            int m = this.size2;
            Debug.Assert(m == B.size1, "dimension mismatch");

            
            matrix P = new matrix(this.size1, B.size2);
            mapto(P, (i, j) => {
                double res = 0;
                for (int k = 0; k < m; k++) res += this[i, k] * B[k, j];
                return res;
            });
            return P;
        }

        public matrix add(matrix B)
        {
            matrix S = new matrix(size1, size2);
            mapto(S, (i, j) => this[i, j] + B[i, j]);
            return S;
        }

        public vector times(vector b)
        {
            return vector.matrixToVector(times(b.getMatrix()));
        }

        public matrix multiply_scalar(double c)
        {
            matrix P = new matrix(this.size1, this.size2);
            mapto(P, (i, j) => this[i, j] * c);
            return P;
        }

        public vector get_column(int j)
        {
            int m = size1;
            vector a = new vector(m);
            for (int i = 0; i < m; i++) a[i] = this[i, j];
            return a;
        }

        public matrix get_row(int i)
        {
            int n = size2;
            matrix a = new matrix(1, n);
            for (int j = 0; j < n; j++) a[0, j] = this[i, j];
            return a;
        }

        public static matrix getIdentity(int dim)
        {
            matrix I = new matrix(dim, dim);
            for (int i = 0; i < dim; i++) I[i, i] = 1;
            return I;
        }

        // QR -decompositioon by the Gram–Schmidt process
        public Tuple<matrix, matrix> QR_decomposition_Gram_Schmidt()
        {
            int m = this.size1, n = this.size2;
            matrix Q_matrix = new matrix(m, n);
            for (int j = 0; j < n; j++)
            {
                vector aj = this.get_column(j);
                for (int i = 0; i < m; i++) Q_matrix[i, j] = aj[i];

                for (int ju = 0; ju < j; ju++)
                {
                    vector proj = aj.project(Q_matrix.get_column(ju));
                    
                    for (int i = 0; i < m; i++) Q_matrix[i, j] -= proj[i];

                }
            }

            for (int j = 0; j < n; j++)
            {
                vector uj = Q_matrix.get_column(j);
                vector uj_normalized = uj.normalize();
                for (int i = 0; i < m; i++) Q_matrix[i, j] = uj_normalized[i];
            }

            return new Tuple<matrix, matrix> (Q_matrix, triangulate(Q_matrix));
        }

        public Tuple<matrix, matrix> QR_decomposition_Givens()
        {
            matrix R = this.getCopy();
            matrix Q = matrix.getIdentity(size1);
            int n, m, plim, qlim;

            for (int p= 0; p < Math.Min(size1, size2); p++)
            {
                for (int q = p + 1; q < size1; q++)
                {
                    if (q == size2 - 1 && q == size1 - 1 ) break;

                    double app = R[p, p];
                    double aqp = R[q, p];
                    double signpp = Math.Sign(app), signqp = Math.Sign(aqp);
                    double absx = signpp * app;
                    double z = Math.Max(absx, signqp * aqp);
                    double t, c, s;

                    if (z == absx)
                    {
                        t = aqp / app;
                        c = signpp / Math.Sqrt(1 + t * t);
                        s = c * t;
                    }
                    else
                    {
                        t = app / aqp;
                        s = signqp / Math.Sqrt(1 + t * t);
                        c = s * t;
                    }
                    JMD.Jtimes(R, p, q, c, s);
                    JMD.timesJ(Q, p, q, c, -s);
                }
            }

            return new Tuple<matrix, matrix> (Q.resize(size1, size2), R.resize(size2, size2));

        }

        // QR - decomposition method.
        public matrix inverse()
        {
            Debug.Assert(size1 == size2, "Attempted to invert non-square matrix.");
            matrix this_inv = getIdentity(size1);
            Tuple<matrix, matrix> res = QR_decomposition_Givens();
            matrix Qt = res.Item1.transpose();
            matrix R = res.Item2;
            for (int j = 0; j < size2; j ++)
            {
                vector b = this_inv.get_column(j);
                vector c = solve_uppertriangular_equation(R, Qt.times(b));
                for (int i = 0; i < size1; i++) this_inv[i, j] = c[i]; 
            }
            return this_inv;
            
        }

        private matrix triangulate(matrix Q_norm)
        {
            matrix R = new matrix(size2, size2);
            mapto(R, (i, j) =>
            {
                if (j >= i) return Q_norm.get_column(i).inner_product(get_column(j));
                return 0;
            });
            return R;
        }

        private static void singular_value_decomposition(double[,] A)
        {
            Debug.Assert(A.GetLength(0) >= A.GetLength(1), "dimension mismatch");

        }

        public static vector solve_uppertriangular_equation(matrix R, vector b)
        {
            int n = b.size;
            for (int j = n - 1; j >= 0; j--)
            {
                b[j] = b[j] / R[j, j];
                R[j, j] = 1;
                for (int k = 0; k < j; k++)
                {
                    b[k] -= R[k, j] * b[j];
                    R[k, j] = 0;
                }
            }
            return b;
        }

        public static void leastsq_matrix(matrix A, Func<double, double>[] fs, double[] x_vals, double[] y_errors)
        {
            mapto(A, (i, j) => fs[j](x_vals[i]) / y_errors[i]);
        }


        public static void calculate_covariance(matrix dcsdb_matrix, matrix covariance_matrix)
        {
            mapto(covariance_matrix, (i, j) => dcsdb_matrix.get_column(i).inner_product(dcsdb_matrix.get_column(j)));
        }

        public static matrix linspaces(double[] alphas, double[] betas, int N)
        {
            Debug.Assert(alphas.Length == betas.Length);
            matrix m = new matrix(N, alphas.Length);
            vector hs = new vector(alphas.Length);

            mapto(m, (i, j) => {
                if (i == 0)
                {
                    hs[j] = (betas[j] - alphas[j]) / (N - 1);
                    return alphas[j];
                }
                else if (i == N - 1) return betas[j];

                return m[i - 1, j] + hs[j];
            });
            return m;
        }

        public matrix extend(vector v)
        {
            int d = Math.Max(v.size, size1);
            matrix res = new matrix(d, size2 + 1);
            mapto(res, (i, j) => {
                if (j < size2 && i < size1) return this[i, j];
                else if (j == size2 && i < v.size) return v[i];
                return double.NaN;
            });
            return res;
        }
    }


    public class vector 
    {
        private matrix vecm;
        public readonly int size;
        public vector(int n)
        {// constructor
            vecm = new matrix(n, 1);
            size = n;
        }

        public vector(double[] xs)
        {// constructor
            int n = xs.Length;
            vecm = new matrix(n, 1);
            size = n;
            setData(xs);
        }

        public void setData(double[] data)
        {// constructor
            vecm.setData(data);
        }

        public double[] getData()
        {// constructor
            double[] data = new double[size];
            for (int i = 0; i < size; i++) data[i] = vecm[i, 0];
            return data;
        }

        public double this[int i]
        { // indexer
            get => vecm[i, 0];      // getter: for example, aij=A[i,j];
            set => vecm[i, 0] = value;// setter: for example, A[i,j]=5;
        }

        new public string ToString()
        {
            return vecm.ToString();
        }

        public bool isEqualGlobal(vector b)
        {
            return this.getMatrix().isEqualGlobal(b.getMatrix());
        }

        public double inner_product(vector v2)
        {
            return getMatrix().transpose().times(v2.getMatrix())[0, 0];
        }
        public matrix outer_product(vector v2)
        {
            return getMatrix().times(v2.getMatrix().transpose());
        }

        public vector project(vector vecU)
        {
            return vecU.multiply_scalar(vecU.inner_product(this) / vecU.inner_product(vecU));
        }

        public matrix transpose()
        {
            return vecm.transpose();
        }

        public vector multiply_scalar(double c)
        {
            return matrixToVector(this.getMatrix().multiply_scalar(c));
        }

        public vector normalize()
        {
            return this.multiply_scalar(1 / norm());
        }

        public double norm()
        {
            return Math.Sqrt(inner_product(this));
        }

        public matrix getMatrix()
        {
            return vecm.getCopy();
        }

        public vector getCopy()
        {
            return matrixToVector(this.getMatrix().getCopy());
        }

        public vector add(vector b)
        {
            return matrixToVector(getMatrix().add(b.getMatrix()));
        }

        public static vector linspace(double a, double b, int N)
        {
            vector v = new vector(N);
            v.vecm = matrix.linspaces(new double[]{ a}, new double[]{ b}, N);
            return v;
        }

        public static vector matrixToVector(matrix vecm)
        {
            Debug.Assert(vecm.size2 == 1, "not a vector");
            vector v_new = new vector(vecm.size1);
            v_new.vecm = vecm;
            return v_new;
        }
    }
}

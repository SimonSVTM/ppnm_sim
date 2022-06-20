using System;
using System.Diagnostics;
using System.Linq;

namespace homework
{
    public class leastsq
    {
        public static Tuple<vector, matrix> QR_leastsq_fit(Func<double, double>[] fs, double[] x_vals, double[] y_vals, double[] y_errors)
        {

            int npoints = x_vals.Length, no_funcs = fs.Length;
            double db = 0.0001;
            matrix dcsdb_matrix = new matrix(no_funcs, no_funcs);
            matrix covariance_matrix = new matrix(no_funcs, no_funcs);
            vector b = new vector(npoints);
            vector c_array_nochange = new vector(no_funcs);

            for (int i = 0; i < npoints; i++) b[i] = y_vals[i] / y_errors[i];

            for (int k = 0; k < no_funcs + 1; k++)
            {
                
                matrix A = new matrix(npoints, no_funcs);
                matrix.leastsq_matrix(A, fs, x_vals, y_errors);
                Tuple < matrix, matrix > res = A.QR_decomposition_Givens();
                //Tuple<matrix, matrix> res = A.QR_decomposition_Gram_Schmidt();
                matrix Q = res.Item1;
                matrix R = res.Item2;

                Debug.Assert(A.isEqualGlobal(Q.times(R)), "QR-decomposition unsuccesful");
                Debug.Assert(A.isEqualGlobal(Q.times(R)), "QR-decomposition unsuccesful");
                Debug.Assert(matrix.getIdentity(Q.size2).isEqualGlobal(Q.transpose().times(Q)), "Q is not unitary");

                Debug.Assert(R.isUpperTriangulary(), "R is not upper triangulary");
                

                matrix R_copy = R.getCopy();
                vector b_copy = b.getCopy();
                vector c = matrix.solve_uppertriangular_equation(R, Q.transpose().times(b));

                b = b_copy;
                Debug.Assert(R_copy.times(c).isEqualGlobal(Q.transpose().times(b)), "Linear equation not solved");

                if (k != 0) {
                    for (int i = 0; i < no_funcs; i++) dcsdb_matrix[i, k - 1] = (c[i] - c_array_nochange[i]) / db;
                    b[k - 1] -= db;
                }
                else c_array_nochange = c.getCopy();
                b[k] += db;
            }


            matrix.calculate_covariance(dcsdb_matrix, covariance_matrix);

            return new Tuple<vector, matrix>(c_array_nochange, covariance_matrix);

        }


    }
}

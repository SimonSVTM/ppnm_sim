using System;
using System.Windows.Forms;
using System.Diagnostics;
using System.Linq;


namespace homework
{
    static class Program
    {
        private static Random rand;

        [STAThread]
        static void Main(string[] args)
        {
            rand = new Random();
            Application.EnableVisualStyles();
            Application.SetCompatibleTextRenderingDefault(false);
            Application.Run(new Form());

        }

        public static void test_jmd()
        {

            int dim = 15;

            for (int iterator = 0; iterator < 100; iterator++)
            {
                matrix A = matrix.getrandomsymmetric(dim, dim, rand);

                matrix A_init = A.getCopy();

                matrix V = matrix.getIdentity(dim);
                JMD.jacobi_cyclic(A, V);
                matrix Vt = V.transpose();

                matrix D = Vt.times(A_init.times(V));
                Debug.Assert(V.times(D.times(Vt)).isEqualGlobal(A_init), "Error VDV^T != A");

                for (int i = 0; i < dim; i++) D[i, i] = 1;

                //Console.WriteLine(String.Join(" ", D.Cast<double>()));
                Debug.Assert(D.isEqualGlobal(matrix.getIdentity(dim)), "D is not diagonal");
                Debug.Assert(Vt.times(V).isEqualGlobal(matrix.getIdentity(dim)), "V not unitary");
                Debug.Assert(V.times(Vt).isEqualGlobal(matrix.getIdentity(dim)), "V not unitary");
            }
            Debug.WriteLine("JMD Test Completed!");
        }

        public static Tuple<double[], matrix> hydrogen_jmd(double dr)
        {
            double rmax = 20;
            int npoints = (int)(rmax / dr) - 1;
            double[] r = new double[npoints];
            for (int i = 0; i < npoints; i++) r[i] = dr * (i + 1);
            matrix H = new matrix(npoints, npoints);
            for (int i = 0; i < npoints - 1; i++)
            {
                H[i, i] = -2;
                H[i, i + 1] = 1;
                H[i + 1, i] = 1;
            }
            H[npoints - 1, npoints - 1] = -2;
            H = H.multiply_scalar(-0.5 / dr / dr);
            for (int i = 0; i < npoints; i++) H[i, i] += -1 / r[i];

            matrix V = matrix.getIdentity(npoints);
            matrix H_init = H.getCopy();
            JMD.jacobi_cyclic(H, V);
            matrix Vt = V.transpose();
            matrix D = Vt.times(H_init.times(V));
            double[] eigenvalues = new double[npoints];
            for (int i = 0; i < npoints; i++) eigenvalues[i] = D[i, i];
            return new Tuple<double[], matrix>(eigenvalues, V);

        }



    }




}




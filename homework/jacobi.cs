using System;
using System.Diagnostics;
using System.Collections;
namespace homework
{
    /*
    JACOBI MATRIX DIAGONALISATION
    */


    public partial class JMD
    {

        public static void ROT2D(matrix A, int[] indeces, double c, double s)
        {
            int i1 = indeces[0], i2 = indeces[1], i3 = indeces[2], i4 = indeces[3];
            double x = A[i1, i2], y = A[i3, i4];
            A[i1, i2] = c * x -s * y;
            A[i3, i4] = s * x + c * y;
        }

        public static void Jtimes(matrix A, int p, int q, double c, double s)
        {
            for (int j = 0; j < A.size2; j++) ROT2D(A, new int[] { p, j, q, j}, c, -s);

        }

        public static void timesJ(matrix A, int p, int q, double c, double s)
        {

            for (int i = 0; i < A.size1; i++) ROT2D(A, new int[] { i, p, i, q }, c, s);
        }



        public static void timesJ(matrix A, int p, int q, double theta)
        {
            timesJ(A, p, q, Math.Cos(theta), Math.Sin(theta));
        }

        public static void Jtimes(matrix A, int p, int q, double theta)
        {
            Jtimes(A, p, q, Math.Cos(theta), Math.Sin(theta));
        }


        public static void jacobi_cyclic(matrix A, matrix V)
        {
            int m = A.size1, n = A.size2;

            bool changed;
            do
            {
                changed = false;
                for (int p = 0; p < m - 1; p++)
                    for (int q = p + 1; q < n; q++)
                    {
                        double apq = A[p, q];
                        double app = A[p, p];
                        double aqq = A[q, q]; 
                        double theta = 0.5 * Math.Atan2(2 * apq, aqq - app);
                        double c = Math.Cos(theta), s = Math.Sin(theta);
                        double new_app = c * c * app - 2 * s * c * apq + s * s * aqq;
                        double new_aqq = s * s * app + 2 * s * c * apq + c * c * aqq;
                        if (new_app != app || new_aqq != aqq) // do rotation
                        {
                            changed = true;
                            timesJ(A, p, q, theta);
                            Jtimes(A, p, q, -theta); // A←J^T*A*J 
                            timesJ(V, p, q, theta); // V←V*J
                        }
                    }
            } while (changed);
        }
    }
}

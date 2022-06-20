using System;
using System.Diagnostics;

namespace homework
{
    class montecarlo
    {

        private static Random rand;
        private static vector x_rand;
        private static int n;
        private static void inititateRandom(string discrepancy, int dim)
        {
            switch (discrepancy)
            {
                case "low":
                    n = 0;
                    x_rand = new vector(dim); 
                    halton(n, x_rand); break;
                case "pure": rand = new Random(); break;
            }
        }

        private static double nextDouble(string discrepancy, int k)
        {
            switch (discrepancy)
            {
                case "low":

                    if (k == 0)
                    {
                        halton(n, x_rand);
                        n++;
                    }
                    return x_rand[k];
                case "pure": return rand.Next() / (double) Int32.MaxValue;
                default: return double.NaN;
            }
        }

        private static double corput(int n, int b )
        {
            double q = 0, bk = (double)1 / b;
            while (n > 0) { q += (n % b ) * bk; n /= b; bk /= b; }
            return q;
        }

        public static void halton(int n, vector x)
        {
            int[] b = new int[] { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67 };
            int maxd = b.Length; 
            int bas;
            for (int i = 0; i < x.size; i++)
            {
                if (i >= maxd) bas = b[i % maxd] * b[maxd - i % maxd - 1];
                else bas = b[i];
                x[i] = corput(n, bas);
            }
        }


        public static (double, double) plainmc(Func<vector, double> f, vector a, vector b, int N = 100000, string discrepancy = "pure")
        {
            inititateRandom(discrepancy, a.size);
            int dim = a.size; double V = 1; for (int i = 0; i < dim; i++) V *= b[i] - a[i];
            double sum = 0, sum2 = 0;
            var x = new vector(dim);
            for (int i = 0; i < N; i++)
            {
               // Debug.WriteLine(nextDouble());
                for (int k = 0; k < dim; k++) x[k] = a[k] + nextDouble(discrepancy, k) * (b[k] - a[k]);
                double fx = f(x); sum += fx; sum2 += fx * fx;
            }
            double mean = sum / N, sigma = Math.Sqrt(sum2 / N - mean * mean);
            var result = (mean * V, sigma * V / Math.Sqrt(N));
            return result;
        }
    }


}

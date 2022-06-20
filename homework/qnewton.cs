using System;
using System.Diagnostics;


namespace homework
{
    class qnewton
    {

        public static vector gradiant(Func<vector, double> f, vector point, double scale = 1)
        {
            double current_value = f(point);
            double dx = scale * 0.0001;
            int n = point.size;
            vector gradf = new vector(n);
            for (int i = 0; i < n; i ++)
            {
                point[i] += dx;
                gradf[i] = (f(point) - current_value) / dx;
                point[i] -= dx;
            }
            return gradf;
        }

        private static (double, bool) linesearch(Func<vector, double> f, vector x, vector dx, vector neggradf)
        {
            double lambda = 1f;
            double fx = f(x);
            int maxiter = 100000;
            int endit = maxiter - 1;
            for (int i = 0; i < maxiter; i++)
            {
                lambda /= 2f;
                if (f(x.add(dx.multiply_scalar(lambda))) < fx - 0.001 * lambda * dx.inner_product(neggradf) ) { endit = i; break; } //Armijo condition
            }
           
            return (lambda, endit == maxiter - 1); 
        }

        private static (double, bool) linesearch_root(Func<vector, vector> f, vector x, vector dx, double fnorm)
        {
            double lambda = 1f;
            int maxiter = 100000;
            int endit = maxiter - 1;
            for (int i = 0; i < maxiter;i++)
            {
                lambda /= 2f;
                if (f(x.add(dx.multiply_scalar(lambda))).norm() < (1 - lambda / 2f) * fnorm || lambda < 1f / 1028f) { endit = i; break; };
            }

            return (lambda, endit == maxiter - 1);
        }

        public static vector root(
            Func<vector, vector> f, /* objective function */
            vector start, /* starting point */
            double acc = 0.00001 /* accuracy goal, on exit |f(x)| should be < acc */
        )
        {
            int n = start.size;

            vector x = start.getCopy(), fcurrent = f(x), dx = fcurrent.multiply_scalar(-1);
            matrix J = matrix.getIdentity(n), Jinv;
            double fnorm = fcurrent.norm();
            (double lambda, _) = linesearch_root(f, x, dx, fnorm);
            dx = dx.multiply_scalar(lambda);
            x = x.add(dx);
            double eps = Math.Pow(2, -26);
            bool failed;
            while (fnorm > acc && dx.norm() > x.norm() * eps)
            {
                vector val = f(x.add(dx)).add(f(x).multiply_scalar(-1)).add(J.times(dx.multiply_scalar(-1)));
                matrix dJ = val.outer_product(val).multiply_scalar(1f / val.inner_product(dx));
                J = J.add(dJ);
                Jinv = J.inverse();
                dx = Jinv.times(fcurrent.multiply_scalar(-1));
                (lambda, failed) = linesearch_root(f, x, dx, fnorm);
                dx = dx.multiply_scalar(lambda);
                double dxnorm = dx.norm();
                if (Double.IsInfinity(dxnorm) || Double.IsNaN(dxnorm) || failed)
                {
                    dx = fcurrent.multiply_scalar(-1);
                    J = matrix.getIdentity(n);
                    (lambda, _) = linesearch_root(f, x, dx, fnorm);
                    dx = dx.multiply_scalar(lambda);
                }
                x = x.add(dx);
                fcurrent = f(x);
                fnorm = fcurrent.norm();
                
            }

            return x;
        }

        public static vector minimize(
            Func<vector, double> f, /* objective function */
            vector start, /* starting point */
            double acc = 0.00001 /* accuracy goal, on exit |gradient| should be < acc */
        )
        {
            
            int n = start.size;
            vector x = start.getCopy(), dx, neggradf;
            neggradf = gradiant(f, x).multiply_scalar(-1);
            dx = neggradf;

            double lambda;
            bool failed;
            matrix H_inv = matrix.getIdentity(n);
            (lambda, _) = linesearch(f, x, dx, neggradf);
            dx = dx.multiply_scalar(lambda);
            x = x.add(dx);

            int nosteps = 0;
            while (neggradf.norm() > acc && lambda > acc * 0.00001)
            {
                neggradf = gradiant(f, x).multiply_scalar(-1);
                

                vector y = gradiant(f, x.add(dx)).add(neggradf);
                vector v = dx.add(H_inv.times(y).multiply_scalar(-1));
                //SR1 update
                H_inv = H_inv.add(v.outer_product(v).multiply_scalar(1f / (v.inner_product(y))));

                dx = H_inv.times(neggradf).normalize();

                (lambda, failed) = linesearch(f, x, dx, neggradf);

                dx = dx.multiply_scalar(lambda);
                double dxnorm = dx.norm();
                if (Double.IsInfinity(dxnorm) || Double.IsNaN(dxnorm) || failed)
                {
                    dx = neggradf;
                    H_inv = matrix.getIdentity(n);
                    (lambda, _) = linesearch(f, x, dx, neggradf);
                    dx = dx.multiply_scalar(lambda);
                }

                x = x.add(dx);
                nosteps++;
            }

            Debug.WriteLine("Converged in " + nosteps.ToString() + " steps");
            return x;

        }
    }


}

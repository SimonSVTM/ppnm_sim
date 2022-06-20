using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;


namespace homework
{
	
    class ode_rke
    {
		private static int i = 1;
		private static Dictionary<int, (vector, matrix, vector)> butcher_tableau;

		private static readonly double[] A45 = 
			new double[] { 
			0, 0, 0, 0, 0, 0,
			(double)1/4, 0, 0, 0, 0, 0 ,
			(double)3/32, (double)9/32, 0, 0, 0, 0,
			(double)1932/2197, -(double)7200/2197, (double)7296/2197, 0, 0, 0,
			(double)439/216, -(double)8, (double)3680/513, -(double)845/4104, 0, 0,
			-(double)8/27, 2, -(double)3544/2565, (double)1859/4104, -(double)11/40, 0};
		private static readonly double[] c45 = new double[] {0, (double)1 /4, (double)3 /8, (double)12 /13, 1, (double)1 /2};
		private static readonly double[] b5 = new double[] { (double)16 /135, 0, (double)6656 /12825, (double)28561 /56430, -(double)9 /50, (double)2 /55 };
		private static double[] b4 = new double[] { (double)24 /216, 0, (double)1408 /2565, (double)2197 /4104, -(double)1 /5, 0 };

		private double frac(int y, int x)
		{
			return y / x;
		}

		public static (vector, matrix) getCA(int l)
        {
			int order = l * (l + 1) / 2 + 1;
			vector c = new vector(order);
			matrix A = new matrix(order, order);
			c[0] = 0;
			int p = 1;
			int u = 0;
			double dl = 1 / (double) l;
 			for (int i = 1; i < order; i ++)
            {
				
				if (c[i - 1] * l >= l - 1)
                {
					c[i] = p * dl;
					p++;
					u = i - 1;
				}
				else c[i] = c[i-1] + dl;
				A[i, u] = c[i];
				for (int j = u + 1; j < i; j++) A[i, j] = dl;
			}

			return (c, A);
        }

		public static (vector, matrix, vector) estimateBvector(double h, int l)
        {
			int order = l * (l + 1) / 2 + 1;
			Func<double, vector> f_integral = (x) =>
			{
				vector v = new vector(1);
				v.setData(new double[] { Math.Exp(- Math.Pow(x, 2) / 2)});
				return v;
			};
			Func<double, vector, vector> f = (x, y) =>
			{
				vector v = new vector(1);
				v.setData(new double[] {- x * y[0] });
				return v;
			};
			(vector c, matrix A) = getCA(l);
			int added = 0;
			vector ydiff = new vector(order + added + 1);
			matrix k_values = new matrix(order + added + 1, order);
			vector y = new vector(1);
			for (int t = 0; t < order + added; t++)
            {
				double x = t * h;
				y = f_integral(x);
				
				matrix Ks = calcKs(f, x, y, h, c, A);
				for (int j = 0; j < order; j++) k_values[t, j] = Ks[0, j];
				ydiff[t] = f_integral(x + h)[0] - f_integral(x)[0];
				
			}
			for (int j = 0; j < order; j++) k_values[order + added, j] = 1;
			ydiff[order + added] = 1;
			//Debug.WriteLine(k_values.ToString());
			Tuple<matrix, matrix> res = k_values.QR_decomposition_Givens();
			matrix Q = res.Item1;
			matrix R = res.Item2;
			//Debug.WriteLine(ydiff.ToString());
			//Debug.WriteLine(k_values.inverse().ToString());
			vector b = matrix.solve_uppertriangular_equation(R, Q.transpose().times(ydiff));
			
			return (c, A, b); 
		}

		public static matrix calcKs(Func<double, vector, vector> f, double x, vector y, double h, vector c, matrix A)
        {
			int order = c.size;
			int ysize = y.size;
			matrix ks = new matrix(ysize, order);
			for (int o = 0; o < order; o++)
            { 
				vector y_new = y.getCopy();
				for (int j = 0; j < o; j++) y_new = y_new.add(ks.get_column(j).multiply_scalar(A[o, j]));
				
				vector Ko = f(x + c[o] * h, y_new).multiply_scalar(h);
				for (int i = 0; i < ysize; i++) ks[i, o] = Ko[i];
            }

			//Debug.WriteLine(ks.ToString());
			return ks;
        }

		public static List<vector> ode(Func<double, vector, vector> f, 
			Func<Func<double, vector, vector>, double, vector, double, (vector, vector)> nummethod, 
			double a, vector ya, double b, double dx)
        {
			List<vector> ys = new List<vector>();
			while (a < b)
            {
				ya = driver(f, nummethod, a, ya, a + dx, dx / 2);
				ys.Add(ya);
				a += dx;
            }
			return ys;

		}

		public static vector driver(
		Func<double, vector, vector> f, /* the f from dy/dx=f(x,y) */
		Func<Func<double, vector, vector>, double, vector, double, (vector, vector)> nummethod,
		double a,                     /* the start-point a */
		vector ya,                    /* y(a) */
		double b,                     /* the end-point of the integration */
		double h = 0.01,                  /* initial step-size */
		double acc = 0.01,              /* absolute accuracy goal */
		double eps = 0.01               /* relative accuracy goal */
	)
		{
			if (a > b) throw new Exception("driver: a>b");
			double x = a; vector y = ya;
			do
			{
				if (x >= b) return y; /* job done */
				if (x + h > b) h = b - x;   /* last step should end at b */
				var (yh, erv) = nummethod(f, x, y, h);
				double tol = Math.Max(acc, yh.norm() * eps) * Math.Sqrt(h / (b - a));
				double err = erv.norm();
				if (err <= tol) { x += h; y = yh; } // accept step
				h *= Math.Min(Math.Pow(tol / err, 0.25) * 0.95, 2); // reajust stepsize
			} while (true);
		}//driver

		public static void initiate_tableau() 
		{
			butcher_tableau = new Dictionary<int, (vector, matrix, vector)>();
		}

		public static (vector, vector) rkstepLK(
		Func<double, vector, vector> f, double x, vector y, double h, int l1, int l2)
		{
			vector c1, c2, b1, b2;
			matrix A1, A2;
			(vector, matrix, vector) res1, res2;
			if (!butcher_tableau.TryGetValue(l1, out res1))
			{
				(c1, A1, b1) = estimateBvector(h, l1);
				butcher_tableau[l1] = (c1, A1, b1);
			}
			else (c1, A1, b1) = res1;

			if (!butcher_tableau.TryGetValue(l2, out res2))
			{
				(c2, A2, b2) = estimateBvector(h, l2);
				butcher_tableau[l2] = (c2, A2, b2);

			}
			else (c2, A2, b2) = res2;


			matrix Ks1 = calcKs(f, x, y, h, c1, A1);
			matrix Ks2 = calcKs(f, x, y, h, c2, A2);

			return rkstep_estimates(y, Ks1, Ks2, b1.getData(), b2.getData());
		}

		public static (vector, vector) rkstep45(
		Func<double, vector, vector> f, double x, vector y, double h)
		{

			vector c = new vector(6);
			c.setData(c45);
			matrix A = new matrix(6, 6);

			A.setData(ode_rke.A45);

			matrix Ks = calcKs(f, x, y, h, c, A);
			
			return rkstep_estimates(y, Ks, Ks, b4, b5);
		}


		private static (vector, vector) rkstep_estimates(vector y, matrix Ks1, matrix Ks2, double[] b1, double[] b2)
        {
			vector y_new1 = y.getCopy();
			vector y_new2 = y.getCopy();
			for (int j = 0; j < b1.Length; j++)
			{
				vector Kj1 = Ks1.get_column(j);
				y_new1 = y_new1.add(Kj1.multiply_scalar(b1[j]));
				
			}
			for (int j = 0; j < b2.Length; j++)
			{
				vector Kj2 = Ks2.get_column(j);
				y_new2 = y_new2.add(Kj2.multiply_scalar(b2[j]));
			}
			vector er = y_new2.add(y_new1.multiply_scalar(-1));

			return (y_new2, er);
		}


		public static (vector, vector) rkstep12(
		Func<double, vector, vector> f, /* the f from dy/dx=f(x,y) */
		double x,   /* the current value of the variable */
		vector y,   /* the current value y(x) of the sought function */
		double h    /* the step to be taken */
		)
		{ // Runge-Kutta Euler/Midpoint method (probably not the most effective)
			vector k0 = f(x, y); /* embedded lower order formula (Euler) */
			vector k1 = f(x + h / 2, y.add(k0.multiply_scalar(h / 2))); /* higher order formula (midpoint) */
			vector yh = y.add(k1.multiply_scalar(h));     /* y(x+h) estimate */
			vector er = k1.add(k0.multiply_scalar(-1)).multiply_scalar(h);  /* error estimate */
			return (yh, er);
		}
	}


}

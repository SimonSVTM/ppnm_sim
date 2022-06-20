using System;
using System.Diagnostics;

namespace homework
{
    class splines
    {
		public static double linterp(vector x, vector y, double z)
		{
			int i = binsearch(x, z);
			double dx = x[i + 1] - x[i]; if (!(dx > 0)) throw new Exception("uups...");
			double dy = y[i + 1] - y[i];
			return y[i] + dy / dx * (z - x[i]);
		}

		public static int binsearch(vector x, double z)
		{/* locates the interval for z by bisection */
			if (!(x[0] <= z && z <= x[x.size - 1])) throw new Exception("binsearch: bad z");
			int i = 0, j = x.size - 1;
			while (j - i > 1)
			{
				int mid = (i + j) / 2;
				if (z > x[mid]) i = mid; else j = mid;
			}
			return i;
		}

		public static double linterpInteg(vector x, vector y, double z)
        {
			if (!(x[0] <= z && z <= x[x.size - 1])) throw new Exception("linterpInteg: bad z");
			double res = 0; 
			int i = 0;
			while (true)
			{
				double dx, dy;
				if (x[i + 1] >= z)
				{
					dx = x[i + 1] - x[i]; if (!(dx > 0)) throw new Exception("uups...");
					dy = y[i + 1] - y[i];
					double val = z - x[i];
					res += (val * dy / (2 * dx) + y[i]) * val;
					break;
				}
				dx = x[i + 1] - x[i]; if (!(dx > 0)) throw new Exception("uups...");
				dy = y[i + 1] - y[i];
				
				res += (dy / 2 + y[i]) * dx; 
				i++;
			}

			return res;
        }
	}

	class qspline
	{
		vector x, y, b, c;
		public qspline(vector xs, vector ys)
		{
			
			x = xs.getCopy();
			y = ys.getCopy();
			int n = x.size;
			Debug.Assert(n == ys.size, "different lengths of arrays.");
			vector p = new vector(n - 1);
			vector dxs = new vector(n - 1);
			for (int i = 0; i < n - 1; i++)
            {
				double dx = x[i + 1] - x[i]; if (!(dx > 0)) throw new Exception("uups...");
				double dy = y[i + 1] - y[i];
				dxs[i] = dx;
				p[i] = dy / dx;
            }

			b = new vector(n - 1);
			c = new vector(n - 1);
			c[0] = 0;
			for (int i = 1; i < n - 1; i++) c[i] = (p[i] - p[i - 1] -c[i - 1] * dxs[i - 1]) / dxs[i];
			vector c_backp = c.getCopy();
			c_backp[n - 2] *= 1 / 2;
			for (int i = n - 3; i >= 0; i--) c_backp[i] = (p[i + 1] - p[i] - c_backp[i + 1] * dxs[i + 1]) / dxs[i];
			for (int i = 1; i < n - 1; i++) c[i] = (c[i] + c_backp[i]) / 2;
			for (int i = 1; i < n - 1; i++) b[i] = p[i] - c[i] * dxs[i];

			/* store xs and ys; calculate b and c */
		}
		public double spline(double z)
		{/* evaluate the spline */
			int i = splines.binsearch(x, z);
			double val = z - x[i];
			return y[i] + b[i] * val + c[i] * val * val;
		}
		public double derivative(double z) {
			int i = splines.binsearch(x, z);
			double val = z - x[i];
			return  b[i]  + 2 * c[i] * val;
		}
		public double integral(double z)
		{
			double res = 0, val;
			int i = 0;
			while (true)
			{
				if (x[i + 1] >= z)
                {
					val = z - x[i];
					res += (y[i] + b[i] * val / 2 + c[i] * val * val / 3) * val;
					break;
				}
				val = x[i + 1] - x[i];
				res += (y[i]  + b[i] * val / 2  + c[i] * val * val / 3) * val;
				i++;
			}
			
			return res;
			/* evaluate the integral */
		}
	}


	class Akima
	{
		vector x, y, b, c, d;
		public Akima(vector xs, vector ys)
		{

			x = xs.getCopy();
			y = ys.getCopy();
			int n = x.size;
			Debug.Assert(n == ys.size, "different lengths of arrays.");
			vector p = new vector(n - 1);
			vector dxs = new vector(n - 1);
			for (int i = 0; i < n - 1; i++)
			{
				double dx = x[i + 1] - x[i]; if (!(dx > 0)) throw new Exception("uups...");
				double dy = y[i + 1] - y[i];
				dxs[i] = dx;
				p[i] = dy / dx;
			}
			vector weigths = new vector(n - 1);
			for (int i = 1; i < n - 1; i++) weigths[i] = Math.Abs(p[i] - p[i - 1]);
			vector As = new vector(n);

			As[0] = p[0];
			As[1] = (p[0] + p[1]) / 2;
			As[n - 2] = (p[n - 2] + p[n - 3]) / 2;
			As[n - 1] = p[n - 2];
			for (int i = 2; i < n - 2; i++)
            {
				double w = weigths[i + 1] + weigths[i - 1];
				if (w == 0) As[i] = (p[i - 1] + p[i]) / 2;
				else As[i] = (weigths[i + 1] * p[i - 1] + weigths[i - 1] * p[i]) / w;
			}

			b = As.getCopy();
			c = new vector(n - 1);
			d = new vector(n - 1);

			for (int i = 0; i < n - 1; i++) c[i] = (3 * p[i] - 2 * As[i] - As[i + 1]) / dxs[i];
			for (int i = 0; i < n - 1; i++) d[i] = (As[i] + As[i + 1] - 2 * p[i]) / Math.Pow(dxs[i], 2);

			/* store xs and ys; calculate b and c */
		}
		public double spline(double z)
		{/* evaluate the spline */
			int i = splines.binsearch(x, z);
			double val = z - x[i];
			return y[i] + b[i] * val + c[i] * Math.Pow(val, 2) + d[i] * Math.Pow(val, 3);
		}
		public double derivative(double z)
		{
			int i = splines.binsearch(x, z);
			double val = z - x[i];
			return b[i] + (2 * c[i] + 3 * d[i] * val) * val;
		}
		public double integral(double z)
		{
			double res = 0, val;
			int i = 0;
			while (true)
			{
				if (x[i + 1] >= z)
				{
					val = z - x[i];
					res += (y[i] + b[i] * val / 2 + c[i] * Math.Pow(val, 2) / 3 + d[i] * Math.Pow(val, 3) / 4) * val;
					break;
				}
				val = x[i + 1] - x[i];
				res += (y[i] + b[i] * val / 2 + c[i] * Math.Pow(val, 2) / 3 + d[i] * Math.Pow(val, 3) / 4) * val;
				i++;
			}

			return res;
			/* evaluate the integral */
		}
	}
}

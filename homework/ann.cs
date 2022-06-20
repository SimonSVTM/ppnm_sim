using System;
using System.Diagnostics;

namespace homework
{
	class ann
	{
		int n; /* number of hidden neurons */
		Func<double, double> f; /* activation function */
		vector p; /* network parameters */
		public ann(int n, Func<double, double> f) 
		{ 
			this.f = f; this.n = n;
			p = new vector(3 * n);
		}
		public double response(double x)
		{
			double res = 0;
			for (int i = 0; i < n; i++) res += f((x - p[i]) / p[n + i]) * p[2 * n + 1];
			return res;
			/* return the response of the network to the input signal x */
		}

		public void train(vector x, vector y)
		{
			Func<vector, double> func = (vector p) =>
			{
				double cost = 0;
				this.p = p;
				for (int i = 0; i < x.size; i++) cost += Math.Pow(response(x[i]) - y[i], 2);
				Debug.WriteLine(cost);
				return cost;
			};
			vector start = new vector(3 * n);
			montecarlo.halton(10, start);
			p = qnewton.root((vector v) => qnewton.gradiant(func, v), start, acc : 0.1);
			/* train the network to interpolate the given table {x,y}*/
		}


	}
}

using System;
using System.Drawing;
using System.Windows.Forms;
using System.Linq;
using System.Diagnostics;
using OxyPlot;
using OxyPlot.Series;
using OxyPlot.WindowsForms;
using OxyPlot.Legends;
using OxyPlot.Axes;
using System.Collections.Generic;
using System.IO;

namespace homework
{

    public partial class Form : System.Windows.Forms.Form
    {
		public Form()
        {
            InitializeComponent();
        }

        private void Form1_Load(object sender, EventArgs e)
        {
			System.Globalization.CultureInfo customCulture = (System.Globalization.CultureInfo)System.Threading.Thread.CurrentThread.CurrentCulture.Clone();
			customCulture.NumberFormat.NumberDecimalSeparator = ".";
			System.Threading.Thread.CurrentThread.CurrentCulture = customCulture;

			//Program.test_jmd();
			//plot_hydrogen_eigens();
			//radioactivity_leastsq();
			//runge_kutta();
			//main2.test();
			//Main("test_genlist.txt");
			//minimum();
			//higgs();
			//montecarlotest();
			//splines_test();
			//qnewtonroot_test();
			//ann_test();
			quad_test();
		}

		private void plot_hydrogen_eigens()
        {
			double conversion_factor = 27.2114;
			double dr = 0.05;
			Tuple<double[], matrix> res = Program.hydrogen_jmd(dr);
			double[] eigenvalues = res.Item1;
			matrix eigenfunctions = res.Item2;
			Console.WriteLine(String.Join(" ", eigenvalues)); 

			int no_functions = 7;
			Func<int, double>[] fs = new Func<int, double>[no_functions];
			for (int i = 0; i < no_functions; i++) {
				vector eigenfunc = eigenfunctions.get_column(i);
				vector eigenfunc_normalized = eigenfunc.normalize();
				fs[i] = (j) => Math.Pow(eigenfunc_normalized[j], 2);
			}
			matrix xss = matrix.linspaces(Enumerable.Repeat((double ) 0, no_functions).ToArray(), Enumerable.Repeat(21 * dr, no_functions).ToArray(), 21);
			string[] flabels = new string[no_functions];
			for (int i = 0; i < no_functions; i++) flabels[i] = $"E = " + Math.Round(eigenvalues[i] * conversion_factor, 2).ToString() + "eV";

			plotf(fs, xss,  flabels, new int[] { 255, 170, 125 }, "distance", "Normalized Squared Eigenvector Values");

			Func<int, double>[] fs2 = new Func<int, double>[] { (i) => conversion_factor * eigenvalues[i]};

			plotf(fs2, xss, flabels, new int[] { 255, 170, 125 }, "distance", "Energy (eV)");

		}

		public void radioactivity_leastsq()
		{
			double[] time = { 1, 2, 3, 4, 6, 9, 10, 13, 15 };
			double[] activity = { 117, 100, 88, 72, 53, 29.5, 25.2, 15.2, 11.1 };
			double[] dys = { 5, 5, 5, 5, 5, 5, 1, 1, 1, 1 };

			double[] ln_activity = activity.Select(x => Math.Log(x)).ToArray();
			double[] dlny = new double[time.Length];
			for (int i = 0; i < time.Length; i++) dlny[i] = dys[i] / activity[i];

			Func<double, double>[] fs = { (t) => 1, (t) => -t };

			Tuple<vector, matrix> fit_values = leastsq.QR_leastsq_fit(fs, time, ln_activity, dlny);
			vector c = fit_values.Item1;
			matrix covariance_matrix = fit_values.Item2;
			double[] uncertainties = new double[fs.Length];
			for (int i = 0; i < fs.Length; i++) uncertainties[i] = Math.Sqrt(covariance_matrix[i, i]);
			double c01, c11, c02, c12;
			c01 = c[0] - uncertainties[0] / 2;
			c11 = c[1] - uncertainties[1] / 2;
			c02 = c01 + uncertainties[0];
			c12 = c11 + uncertainties[1];


			matrix xss = matrix.linspaces(new double[] { 0, 0, 0 }, new double[] { 15, 15, 15 }, (int) (15 / 0.1));


			for (int i = 0; i < time.Length; i++)
			{
				vector xs = new vector(3);
				double s = time[i];
				xs.setData(new double[] { s, s, s});
				xss = xss.extend(xs);
			}

			Func<double, double>[] funcs = new Func<double, double>[] { (x) => Math.Exp(c[0] - c[1] * x), (x) => Math.Exp(c01 - c11 * x), (x) => Math.Exp(c02 - c12 * x) };
			int tot = funcs.Length + time.Length;
			Func<int, double>[] fss = new Func<int, double>[tot];

			for (int i = 0; i < tot; i++)
			{
				if (i < funcs.Length)
				{
					Func<double, double> f = funcs[i];
					fss[i] = (j) => f(xss[j, 0]);
				}
				else
				{
					double a = activity[i - funcs.Length], b = dys[i - funcs.Length];
					fss[i] = (j) => a + (j - 1) * b / 2;
				}
			}
			string[] captions = new string[tot];
			captions[0] = "Fitted";
			captions[1] = "Error Bands";
			for (int i = 2; i < tot; i++) captions[i] = "";
			captions[6] = "Error Bars";
			plotf(fss, xss, captions, new int[] { 255, 50, 210 }, "Time", "Activity", markersize : 0);


		}


		private void runge_kutta()
		{
			vector c, b_test;
			matrix A;

			int l = 2;
			//(c, A) = ode_rke.getCA(l);
			//Debug.Write("c = \n" + c.ToString());
			//Debug.Write("A = \n" + A.ToString());
			//(c, A, b_test) = ode_rke.estimateBvector(0.1, l);
			//Debug.WriteLine("b = \n" + b_test.ToString());
			Func<double, vector, vector> f = (x, y) => {
				vector dydx = new vector(2);
				dydx.setData(new double[] { y[1], -y[0] });
				return dydx;
			};

			ode_rke.initiate_tableau();
			(vector, vector) rkstep32(Func<double, vector, vector> f, double x, vector y, double h) { return ode_rke.rkstepLK(f, x, y, h, 4, 5); }

			double h = 0.1;
			double a = 0.1, b = 4;

			vector ya = new vector(2);
			ya.setData(new double[] { 0, 1 });
			Func<Func<double, vector, vector>, double, vector, double, (vector, vector)>
				optmethodLK = (f, x, y, h) => ode_rke.rkstepLK(f, x, y, h, 1, 2);

			List<vector> ys = ode_rke.ode(f, ode_rke.rkstep45, a, ya, b, h);
			List<vector> ys2 = ode_rke.ode(f, ode_rke.rkstep12, a, ya, b, h);
			List<vector> ys3 = ode_rke.ode(f, rkstep32, a, ya, b, h);


			Func<int, double>[] fs = new Func<int, double>[] { (i) => ys[i][1], (i) => ys2[i][1], (i) => ys3[i][1] };
			matrix xss = matrix.linspaces(new double[] { a, a, a }, new double[] { b, b, b }, (int)((b - a) / h));
			plotf(fs, xss, new string[] { "RK45", "RK12", "RK45 - numerical tableau" }, new int[]{ 255, 170, 125 }, "x", "y(x)", true);


		}

		private void minimum()
        {
			Debug.WriteLine("STARTED");
			vector start = new vector(2);
			start[0] = 1.2; start[1] = 1.3;
			Debug.WriteLine(qnewton.minimize((vector y) => Math.Pow(1 - y[0], 2) + 100 * Math.Pow(y[1] - y[0] * y[0], 2), start).ToString());
			Debug.WriteLine(qnewton.minimize((vector y) => Math.Pow(y[0] * y[0] + y[1] - 11, 2) + 100 * Math.Pow(y[0] + y[1] * y[1] - 7, 2), start).ToString());
			Debug.WriteLine("DONE");

		}

		public void higgs()
        {
			double Bret_Wigner(double E, double m, double gamma, double A) { return A / (Math.Pow(E - m, 2) + gamma * gamma / 4); }

			list<vector> veclist = list_methods.file_to_list("higgs.txt");
			vector E = veclist.get(0), cross_section = veclist.get(1), dc = veclist.get(2);
			string res = "";
			for (veclist.start(); veclist.current != null; veclist.next()) res += veclist.current.item.ToString();

			Func<vector, double> deviation_func = (vector v) =>
			{
				vector deviation_vec = new vector(E.size);
				for (int i = 0; i < E.size; i++) deviation_vec[i] = (Bret_Wigner(E[i], v[0], v[1], v[2]) - cross_section[i]) / dc[i];
				return deviation_vec.inner_product(deviation_vec);
			};
			vector v = qnewton.minimize(deviation_func, new vector(new double[] { 125, 10, 100}), acc : 0.00000000000000000000001);
			Debug.WriteLine(v.ToString());

			matrix xss = vector.linspace(100, 160, 150).getMatrix();

			for (int i = 0; i < E.size; i++)
			{
				vector xs = new vector(3);
				double s = E[i];
				xs.setData(new double[] { s, s, s });
				xss = xss.extend(xs);
			}

			Func<double, double>[] funcs = new Func<double, double>[] { (E) => Bret_Wigner(E, v[0], v[1], v[2])};
			int tot = funcs.Length + E.size;
			Func<int, double>[] fss = new Func<int, double>[tot];

			for (int i = 0; i < tot; i++)
			{
				if (i < funcs.Length)
				{
					Func<double, double> f = funcs[i];
					fss[i] = (j) => f(xss[j, 0]);
				}
				else
				{
					double a = cross_section[i - funcs.Length], b = dc[i - funcs.Length];
					fss[i] = (j) => a + (j - 1) * b / 2;
				}
			}
			string[] captions = new string[tot];
			captions[0] = "Fitted";
			for (int i = 1; i < tot; i++) captions[i] = "";
			captions[15] = "Error Bars";
			plotf(fss, xss, captions, new int[] { 255, 50, 210 }, "Time", "Activity", markersize: 0);

		}


		private void montecarlotest()
		{

			Func<vector, double> f = (vector v) => {
				if (new vector(new double[] { v[0], v[1] }).norm() <= 1) return 1;
				return 0;
			};

			int count = 100000000;
			(double res, double err) = montecarlo.plainmc(f, new vector(new double[] { -1, -1 }), new vector(new double[] { 1, 1 }), N: count);
			Debug.WriteLine(res);

			(res, err) = montecarlo.plainmc(f, new vector(new double[] { -1, -1 }), new vector(new double[] { 1, 1 }), discrepancy: "low", N: count);
			Debug.WriteLine(res);

			f = (vector v) => 1 / (1 - Math.Cos(v[0]) * Math.Cos(v[1]) * Math.Cos(v[2]));

			(res, err) = montecarlo.plainmc(f, new vector(new double[] { 0, 0, 0 }), new vector(new double[] { Math.PI, Math.PI, Math.PI }));
			Debug.WriteLine(res / Math.Pow(Math.PI, 3));
		}

		private void splines_test()
		{
			Func<double, double> f = (double x) => Math.Sqrt(1 - Math.Pow(x, 2));
			int N = 100;
			vector xs = vector.linspace(-1, 1, N);

			vector ys = new vector(N);
			for (int i = 0; i < xs.size; i++) ys[i] = f(xs[i]);

			Debug.WriteLine(xs[xs.size - 1]);
			double res = splines.linterpInteg(xs, ys, 1);
			Debug.WriteLine(2 * res);

			qspline qspl = new qspline(xs, ys);
			res = qspl.integral(1);
			Debug.WriteLine(2 * res);

			Akima akima = new Akima(xs, ys);
			res = akima.integral(1);
			Debug.WriteLine(2 * res);


			vector xs_new = vector.linspace(-1, 1, N * 100);
			Func<int, double>[] fs = new Func<int, double>[] {(i) => splines.linterp(xs, ys, xs_new[i]), (i) => qspl.spline(xs_new[i]), (i) => akima.spline(xs_new[i]) };
			string[] captions = new string[] { "Linear", "Quadratic", "Akima" };
			matrix xsm = xs_new.getMatrix();
			xsm = xsm.extend(xs_new).extend(xs_new);
			plotf(fs, xsm, captions, new int[] { 255, 50, 210 }, "x", "y", markersize: 0);
		}

		private void qnewtonroot_test()
		{
			Func<vector, double> f = (vector v) => Math.Pow(1 - v[0], 2) + 100 * Math.Pow(v[1] - v[0] * v[0], 2);
			vector x = qnewton.root((vector v) => qnewton.gradiant(f, v), new vector(new double[] { 1.2, 1.3 }));
			Debug.WriteLine(x.ToString());
		}

		private void ann_test()
        {
			int N = 25;
			vector xs = vector.linspace(-1, 1, N);

			ann nn = new ann(N, (x) =>  Math.Exp(-x * x));
			
			Func<double, double> f = (x) => Math.Cos(5 * x - 1) * Math.Exp(-x * x) ;
			vector ys = new vector(N);
			for (int i = 0; i < xs.size; i++) ys[i] = f(xs[i]);
			nn.train(xs, ys);

			Func<int, double>[] fs = new Func<int, double>[] { (i) => ys[i], (i) => nn.response(xs[i]) };
			string[] captions = new string[] { "Real", "Neural Network"};
			matrix xsm = xs.getMatrix();
			xsm = xsm.extend(xs);
			plotf(fs, xsm, captions, new int[] { 255, 50, 210 }, "x", "y", markersize: 0);

		}

		private void quad_test()
        {
			Debug.WriteLine(quadratures.integrate((x) => Math.Sqrt(x), 0, 1));
			Debug.WriteLine(quadratures.integrate((x) => 1f /Math.Sqrt(x), 0, 1));
			Debug.WriteLine(quadratures.integrate((x) => 4f * Math.Sqrt(1 - x * x), 0, 1));
			Debug.WriteLine(quadratures.integrate((x) => Math.Log(x) / Math.Sqrt(x), 0, 1));
			Debug.WriteLine(quadratures.integrate_polar((x) => 1f / Math.Sqrt(x), 0, 1));
			Debug.WriteLine(quadratures.integrate_polar((x) => Math.Log(x) / Math.Sqrt(x), 0, 1));
		}

		public void Main(string location)
		{
			var list = list_methods.file_to_list(location);
			string res = "";
			for (list.start(); list.current != null; list.next()) res += list.current.item.ToString();
			Debug.Write(res);
		}

		private void plotf(Func<int, double>[] fs, matrix xs, string[] flabels, int[] rgb_basecolor, string xlabel = "x", string ylabel = "y", bool scattered = false, int markersize = 4)
        {
			PlotView pv = makePlotview();
			var x_axis1 = new LinearAxis { Position = AxisPosition.Bottom, Title = xlabel, MajorGridlineStyle = LineStyle.Solid, MinorGridlineStyle = LineStyle.Dot};
			x_axis1.Key = xlabel;
			pv.Model.Axes.Add(x_axis1);

			var y_axis1 = new LinearAxis { Position = AxisPosition.Left, Title = ylabel, MajorGridlineStyle = LineStyle.Solid, MinorGridlineStyle = LineStyle.Dot };
			y_axis1.Key = ylabel;
			pv.Model.Axes.Add(y_axis1);

			

			double xmax = double.MinValue, fmax = double.MinValue, xmin = double.MaxValue, fmin = double.MaxValue, x, f;

			for (int i = 0; i < fs.Length; i++)
			{
				OxyColor col = OxyColor.FromRgb(Convert.ToByte(i * rgb_basecolor[0] / fs.Length), Convert.ToByte(i * rgb_basecolor[1] / fs.Length), Convert.ToByte(i * rgb_basecolor[2] / fs.Length));
				LineSeries line1 = getLine(flabels[i], col, scattered, markersize);
				for (int j = 0; j < xs.size1; j++)
				{
					x = xs[j, i];
					if (!double.IsNaN(x))
                    {
						f = fs[i](j);
						if (x < xmin) xmin = x;
						if (x > xmax) xmax = x;
						if (f < fmin) fmin = f;
						if (f > fmax) fmax = f;


						line1.Points.Add(new DataPoint(x, f));
					}

					
				}
				line1.YAxisKey = y_axis1.Key;
				pv.Model.Series.Add(line1);

			}
			
			x_axis1.AbsoluteMinimum = Math.Min(xmin, 0f);
			x_axis1.AbsoluteMaximum = Math.Max(xmax, 0f);

			double dy = Math.Abs(fmax - fmin) * 0.005;
			y_axis1.AbsoluteMinimum = Math.Min(fmin, 0) - dy;
			y_axis1.AbsoluteMaximum = Math.Max(fmax, 0) + dy;

		}





		private LineSeries getLine(string title, OxyColor color, bool scattered, int markersize = 4)
        {
			int thickness = 3;
			if (scattered) thickness = 0;
			return new LineSeries()
			{
				Title = title,
				Color = color,
				StrokeThickness = thickness,
				MarkerSize = markersize,
				MarkerType = MarkerType.Circle

			};
		}

		private PlotView makePlotview()
        {
			PlotView pv = new PlotView();
			pv.Model = new PlotModel();
			pv.Model.Legends.Add(new Legend()
			{
				LegendTitle = "Legend",
				LegendPosition = LegendPosition.RightTop,
			});


			pv.Location = new Point(0, 0);
			pv.Size = new Size(800, 500);
			this.Controls.Add(pv);

			return pv;
		}

    }
}

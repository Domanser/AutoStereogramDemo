using System;
using System.Collections.Generic;
using System.Drawing;
using System.Drawing.Imaging;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Text;

namespace AutoStereogramDemo
{
	class Program
	{
		private class Circle3D
		{
			public Point3D Center { get; set; }
			public double Radius { get; set; }
		}

		private static string GetPathAtAssemblyLocation(string fileName)
		{
			return string.Format("{0}\\{1}", Path.GetDirectoryName(Assembly.GetExecutingAssembly().Location), fileName);
		}

		private static void CreateSpheres(int subpixelsCount)
		{
			var asb = new AutoStereogramBuilder(1024, 768, subpixelsCount);

			asb.AddHorizontalPlane(2);
			asb.AddSphere(asb.Width / 2, asb.Height / 2, 1.7, 0.25);
			asb.AddSphere(asb.Width / 2 - 0.25, asb.Height / 2 - 0.25, 1.85, 0.15);
			asb.AddSphere(asb.Width / 2 + 0.4, asb.Height / 2, 1.7, 0.1);
			asb.AddSphere(asb.Width / 2 + 0.1, asb.Height / 2 - 0.05, 1.3, 0.05);

			asb.GenerateBitmap().Save(GetPathAtAssemblyLocation(string.Format("output\\spheres{0}.gif", subpixelsCount)), ImageFormat.Gif);
		}

		private static void CreateFigures()
		{
			var asb = new AutoStereogramBuilder(1024, 768, 3);

			asb.AddHorizontalPlane(2);

			Vector3D v1 = new Vector3D { X = 1, Y = 0, Z = 0 };

			for (int n = 0; n < 4; n++)
			{
				Vector3D v2 = new Vector3D { X = 0, Y = Math.Cos(n * Math.PI / 8), Z = -Math.Sin(n * Math.PI / 8) };
				Vector3D v3 = new Vector3D { X = 0, Y = Math.Sin(n * Math.PI / 8), Z = Math.Cos(n * Math.PI / 8) };

				// hexagonal prism

				Point3D[] prismUpperPoints = new Point3D[6], prismLowerPoints = new Point3D[6];

				for (int m = 0; m < 6; m++)
				{
					prismUpperPoints[m] = new Point3D { X = -0.2 + n * 0.25, Y = -0.2, Z = 1.7 } + v1 * Math.Sin(m / 3.0 * Math.PI) * 0.1 + v2 * Math.Cos(m / 3.0 * Math.PI) * 0.1 - v3 * 0.05;
					prismLowerPoints[m] = new Point3D { X = -0.2 + n * 0.25, Y = -0.2, Z = 1.7 } + v1 * Math.Sin(m / 3.0 * Math.PI) * 0.1 + v2 * Math.Cos(m / 3.0 * Math.PI) * 0.1 + v3 * 0.05;
				}

				asb.AddPolygon(prismUpperPoints);
				asb.AddPolygon(prismLowerPoints);

				for (int m = 0; m < 6; m++)
					asb.AddPolygon(new Point3D[] { prismLowerPoints[m], prismUpperPoints[m], prismUpperPoints[(m + 1) % 6], prismLowerPoints[(m + 1) % 6] });

				// cylinder

				asb.AddCylinder(new Point3D { X = -0.2 + n * 0.25, Y = 0.1, Z = 1.7 } - v3 * 0.1, new Point3D { X = -0.2 + n * 0.25, Y = 0.1, Z = 1.7 } + v3 * 0.1, 0.08);

				// octahedron

				Point3D octahedronCenter = new Point3D { X = -0.2 + n * 0.25, Y = 0.4, Z = 1.7 };
				Point3D octahedronUpperPoint = octahedronCenter - v3 * 0.1 * Math.Sqrt(2), octahedronLowerPoint = octahedronCenter + v3 * 0.1 * Math.Sqrt(2);
				Point3D[] octahedronMiddlePoints = new Point3D[] { octahedronCenter + v1 * 0.1 + v2 * 0.1, octahedronCenter - v1 * 0.1 + v2 * 0.1, 
				 octahedronCenter - v1 * 0.1 - v2 * 0.1, octahedronCenter + v1 * 0.1 - v2 * 0.1 };

				for (int m = 0; m < 4; m++)
				{
					asb.AddPolygon(new Point3D[] { octahedronUpperPoint, octahedronMiddlePoints[m], octahedronMiddlePoints[(m + 1) % 4] });
					asb.AddPolygon(new Point3D[] { octahedronLowerPoint, octahedronMiddlePoints[m], octahedronMiddlePoints[(m + 1) % 4] });
				}
			}

			asb.GenerateBitmap().Save(GetPathAtAssemblyLocation("output\\figures.gif"), ImageFormat.Gif);
		}

		private static double? WavesSurfaceFunc(double x, double y)
		{
			double r = Math.Sqrt((x + 0.2) * (x + 0.2) + (y + 0.2) * (y + 0.2));
			return 2 - Math.Cos(r * 25) * 0.2 / (1 + r * 3);
		}

		private static List<Point3D> GetCirclePoints(Point3D center, double radius, Vector3D v1, Vector3D v2)
		{
			List<Point3D> points = new List<Point3D>(16);

			for (int n = 0; n < 32; n++)
			{
				double x = Math.Sin(n / 32.0 * Math.PI * 2) * radius;
				double y = Math.Cos(n / 32.0 * Math.PI * 2) * radius;

				points.Add(center + v1 * x + v2 * y);
			}

			return points;
		}

		private static void GetSomeOrthonormalizedVectors(Vector3D ortVec, out Vector3D v1, out Vector3D v2)
		{
			v1 = ortVec.GetSomeOrtVector().Normalize();
			v2 = Vector3D.CrossProduct(ortVec, v1).Normalize();
		}

		private static void Swap<T>(ref T value1, ref T value2)
		{
			T temp = value1;
			value1 = value2;
			value2 = temp;
		}

		private static void ExpressVector(Vector3D vector, Vector3D baseVector1, Vector3D baseVector2, out double coef1, out double coef2)
		{
			double[][] v = new double[][] { new double[] { baseVector1.X, baseVector2.X, vector.X }, new double[] { baseVector1.Y, baseVector2.Y, vector.Y },
			 new double[] { baseVector1.Z, baseVector2.Z, vector.Z } };

			if (Math.Abs(v[1][0]) > Math.Abs(v[0][0]) && Math.Abs(v[1][0]) > Math.Abs(v[2][0]))
				Swap(ref v[1], ref v[0]);
			else if (Math.Abs(v[2][0]) > Math.Abs(v[0][0]) && Math.Abs(v[2][0]) > Math.Abs(v[1][0]))
				Swap(ref v[2], ref v[0]);

			double k = v[1][0] / v[0][0];
			v[1][0] = 0;
			v[1][1] -= v[0][1] * k;
			v[1][2] -= v[0][2] * k;

			k = v[2][0] / v[0][0];
			v[2][0] = 0;
			v[2][1] -= v[0][1] * k;
			v[2][2] -= v[0][2] * k;

			if (Math.Abs(v[2][1]) > Math.Abs(v[1][1]))
				Swap(ref v[2], ref v[1]);

			coef2 = v[1][2] / v[1][1];
			coef1 = (v[0][2] - coef2 * v[0][1]) / v[0][0];
		}

		private static void GetNextOrthonormalizedVectors(Vector3D circle1OrtVec, Vector3D circle2OrtVec, Vector3D circle1V1, Vector3D circle1V2, 
		 out Vector3D circle2V1, out Vector3D circle2V2)
		{
			Vector3D commonVec = Vector3D.CrossProduct(circle1OrtVec, circle2OrtVec);
			if (commonVec.Abs() < 1e-8)
			{
				circle2V1 = circle1V1;
				circle2V2 = circle1V2;
				return;
			}

			commonVec.Normalize();
			Vector3D circle1AuxVec = Vector3D.CrossProduct(commonVec, circle1OrtVec).Normalize();
			Vector3D circle2AuxVec = Vector3D.CrossProduct(commonVec, circle2OrtVec).Normalize();

			double coef1, coef2;
			ExpressVector(circle1V1, commonVec, circle1AuxVec, out coef1, out coef2);
			circle2V1 = (commonVec * coef1 + circle2AuxVec * coef2).Normalize();
			ExpressVector(circle1V2, commonVec, circle1AuxVec, out coef1, out coef2);
			circle2V2 = (commonVec * coef1 + circle2AuxVec * coef2).Normalize();
		}

		private static void CreatePipe(AutoStereogramBuilder asb, List<Circle3D> circles, Vector3D firstDirection = null)
		{
			List<Point3D> prevPoints = null;
			Vector3D v1 = null, v2 = null;

			if (firstDirection == null)
				firstDirection = circles[1].Center - circles[0].Center;

			for (int n = 0; n < circles.Count; n++)
			{
				if (n == 0)
					GetSomeOrthonormalizedVectors(firstDirection, out v1, out v2);
				else
				{
					GetNextOrthonormalizedVectors((n == 1 ? firstDirection : circles[n - 1].Center - circles[n - 2].Center).Normalize(), 
					 (circles[n].Center - circles[n - 1].Center).Normalize(), v1, v2, out v1, out v2);
				}

				List<Point3D> points = GetCirclePoints(circles[n].Center, circles[n].Radius, v1, v2);

				if (prevPoints != null)
				{
					for (int m = 0; m < points.Count; m++)
					{
						Point3D p1 = points[m], p2 = prevPoints[m], p3 = points[(m + 1) % points.Count], p4 = prevPoints[(m + 1) % points.Count];

						asb.AddPolygon(new Point3D[] { p1, p2, p3 });
						asb.AddPolygon(new Point3D[] { p2, p3, p4 });
					}
				}

				prevPoints = points;
			}
		}

		private static void CreatePipe()
		{
			var asb = new AutoStereogramBuilder(1024, 768, 3);

			asb.AddHorizontalPlane(2);

			List<Circle3D> circles = new List<Circle3D>();

			for (double x = -0.18; x < 0.55; x += 0.003)
				circles.Add(new Circle3D { Center = new Point3D { X = x, Y = 0.1 + Math.Sin(x * 50) * 0.2, Z = 1.5 + Math.Cos(x * 50) * 0.2 }, Radius = 0.04 + x * 0.03 });

			CreatePipe(asb, circles);

			asb.GenerateBitmap().Save(GetPathAtAssemblyLocation("output\\pipe.gif"), ImageFormat.Gif);
		}

		private static void CreateObjectsAtTwoSides()
		{
			var asb = new AutoStereogramBuilder(1680, 1050, 3);

			asb.AddHorizontalPlane(2);
			asb.AddHorizontalPlane(1.8, (x, y, z) => (x - asb.Width / 2) * (x - asb.Width / 2) / 0.8 + (y - asb.Height / 2) * (y - asb.Height / 2) / 0.4 <= 1);
			asb.AddHorizontalPlane(asb.Width * 0.35, asb.Width * 0.65, asb.Height * 0.35, asb.Height * 0.65, -0.2);
			asb.AddSphere(asb.Width / 2, asb.Height / 2, -0.25, 0.03);

			asb.GenerateBitmap().Save(GetPathAtAssemblyLocation("output\\twosides.gif"), ImageFormat.Gif);
		}

		private static void CreateSurface()
		{
			var asb = new AutoStereogramBuilder(1024, 768, 3);

			asb.AddSurface((x, y) => 2 + Math.Sin(Math.Sqrt(x * x + y * y) * 40 + Math.Atan2(y, x)) * 0.25 / (4 * Math.Sqrt(x * x + y * y) + 1), 2.25, 1.75);
			asb.GenerateBitmap().Save(GetPathAtAssemblyLocation("output\\surface.gif"), ImageFormat.Gif);
		}

		private static void CreateWavesWithImage()
		{
			var asb = new AutoStereogramBuilder(1024, 768, 3);

			asb.AddSphere(-0.2, -0.2, 1.6, 0.025);
			asb.AddSphere(-0.2, -0.2, 1.4, 0.025);
			asb.AddSurface(WavesSurfaceFunc, 2.3, 1.7);

			Bitmap bitmap = (Bitmap)Image.FromFile(GetPathAtAssemblyLocation("input\\hi.png"));
			asb.AddModelByImage(bitmap, new Point3D { X = 0.1, Y = -0.3, Z = 1.8 }, new Vector3D { X = 1, Y = 0.5, Z = -1 },
			 new Vector3D { X = -0.5, Y = 1, Z = -1 }, new Vector3D { X = 1, Y = 1, Z = 1 }, 0.8, 0.8, 0.05);

			Bitmap backgroundImage = (Bitmap)Image.FromFile(GetPathAtAssemblyLocation("input\\background.jpg"));
			asb.GenerateBitmap(backgroundImage).Save(GetPathAtAssemblyLocation("output\\hello.jpg"), ImageFormat.Jpeg);
		}

		private static void Main(string[] args)
		{
			string outputFolder = GetPathAtAssemblyLocation("output");
			if (!Directory.Exists(outputFolder))
				Directory.CreateDirectory(outputFolder);

			Console.WriteLine("CreateSpheres(1)...");
			CreateSpheres(1);

			Console.WriteLine("CreateSpheres(3)...");
			CreateSpheres(3);

			Console.WriteLine("CreateFigures...");
			CreateFigures();

			Console.WriteLine("CreatePipe...");
			CreatePipe();

			Console.WriteLine("CreateObjectsAtTwoSides...");
			CreateObjectsAtTwoSides();

			Console.WriteLine("CreateSurface...");
			CreateSurface();

			Console.WriteLine("CreateWavesWithImage...");
			CreateWavesWithImage();
		}
	}
}

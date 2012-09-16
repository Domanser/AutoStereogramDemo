using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace AutoStereogramDemo
{
	public class Rectangle2D<T>
	{
		public T X1 { get; set; }
		public T Y1 { get; set; }
		public T X2 { get; set; }
		public T Y2 { get; set; }
	}

	public class Point3D
	{
		public double X { get; set; }
		public double Y { get; set; }
		public double Z { get; set; }

		public static Point3D operator +(Point3D p, Vector3D v)
		{
			return new Point3D { X = p.X + v.X, Y = p.Y + v.Y, Z = p.Z + v.Z };
		}

		public static Point3D operator -(Point3D p, Vector3D v)
		{
			return new Point3D { X = p.X - v.X, Y = p.Y - v.Y, Z = p.Z - v.Z };
		}

		public static Vector3D operator -(Point3D p1, Point3D p2)
		{
			return new Vector3D { X = p1.X - p2.X, Y = p1.Y - p2.Y, Z = p1.Z - p2.Z };
		}
	}

	public class Vector3D
	{
		public double X { get; set; }
		public double Y { get; set; }
		public double Z { get; set; }

		public double Abs()
		{
			return Math.Sqrt(Abs2());
		}

		public double Abs2()
		{
			return X * X + Y * Y + Z * Z;
		}

		public Vector3D Normalize()
		{
			double abs = Abs();
			if (abs < 1e-8)
				return this;

			X /= abs;
			Y /= abs;
			Z /= abs;

			return this;
		}

		public Vector3D GetSomeOrtVector()
		{
			if (Math.Abs(X) >= Math.Abs(Y) && Math.Abs(X) >= Math.Abs(Z))
				return new Vector3D { X = -(Y + Z) / X, Y = 1, Z = 1 };
			else if (Math.Abs(Y) >= Math.Abs(X) && Math.Abs(Y) >= Math.Abs(Z))
				return new Vector3D { X = 1, Y = -(X + Z) / Y, Z = 1 };
			else
				return new Vector3D { X = 1, Y = 1, Z = -(X + Y) / Z };
		}

		public static Vector3D operator +(Vector3D v1, Vector3D v2)
		{
			return new Vector3D { X = v1.X + v2.X, Y = v1.Y + v2.Y, Z = v1.Z + v2.Z };
		}

		public static Vector3D operator -(Vector3D v1, Vector3D v2)
		{
			return new Vector3D { X = v1.X - v2.X, Y = v1.Y - v2.Y, Z = v1.Z - v2.Z };
		}

		public static Vector3D operator *(Vector3D v, double value)
		{
			return new Vector3D { X = v.X * value, Y = v.Y * value, Z = v.Z * value };
		}

		public static Vector3D operator /(Vector3D v, double value)
		{
			return new Vector3D { X = v.X / value, Y = v.Y / value, Z = v.Z / value };
		}

		public static Vector3D CrossProduct(Vector3D v1, Vector3D v2)
		{
			return new Vector3D { X = v1.Y * v2.Z - v1.Z * v2.Y, Y = v1.Z * v2.X - v1.X * v2.Z, Z = v1.X * v2.Y - v1.Y * v2.X };
		}

		public static double ScalarProduct(Vector3D v1, Vector3D v2)
		{
			return v1.X * v2.X + v1.Y * v2.Y + v1.Z * v2.Z;
		}
	}

	public class Plane
	{
		public double A { get; set; }
		public double B { get; set; }
		public double C { get; set; }
		public double D { get; set; }

		public Plane()
		{
		}

		public Plane(Vector3D ortVec, Point3D point)
		{
			A = ortVec.X;
			B = ortVec.Y;
			C = ortVec.Z;
			D = -A * point.X - B * point.Y - C * point.Z;
		}

		public Plane Normalize()
		{
			double abs = Math.Sqrt(A * A + B * B + C * C);
			if (abs < 1e-8)
				return this;

			A /= abs;
			B /= abs;
			C /= abs;
			D /= abs;
			
			return this;
		}
	}
}

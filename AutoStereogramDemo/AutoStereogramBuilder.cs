using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.Text;

namespace AutoStereogramDemo
{
	public class AutoStereogramBuilder
	{
		public delegate double? SurfaceFunc(double x, double y);
		public delegate bool GetRayIntersection(int xProj, int yProj, double eyeXPos, out double x, out double z);
		public delegate bool CheckPoint(double x, double y, double z);
		public delegate void ColorGenerator(Bitmap bitmap, int y, int uniqueColorsCount, int[] colors);

		private struct SurfacePoint
		{
			public double X, Y;
			public double? Z;
			public int XLeftProj, XRightProj, YProj;
		}

		private class SurfaceRenderStackFrame
		{
			public int Step;
			public readonly SurfacePoint[,] Points = new SurfacePoint[3, 3];
		}

		private readonly double eyeLeftXPos, eyeRightXPos, eyeYPos;
		private readonly int?[,] nearestLeft, nearestRight;
		private readonly Random random = new Random();
		
		/// <param name="xResolution">Horizontal resolution of stereogram</param>
		/// <param name="yResolution">Vertical resolution of stereogram</param>
		/// <param name="subpixelsCount">Number of subpixels</param>
		/// <param name="pixelWidth">Width of pixel, in meters</param>
		/// <param name="pixelHeight">Height of pixel, in meters</param>
		/// <param name="distanceToEyes">Distance between stereogram and eyes, in meters</param>
		/// <param name="distanceBetweenEyes">Distance between eyes, in meters</param>
		public AutoStereogramBuilder(int xResolution, int yResolution, int subpixelsCount = 1, double pixelWidth = 0.0003, double pixelHeight = 0.0003, 
		 double distanceToEyes = 0.5, double distanceBetweenEyes = 0.065)
		{
			XResolution = xResolution;
			YResolution = yResolution;
			SubpixelsCount = subpixelsCount;
			PixelWidth = pixelWidth;
			PixelHeight = pixelHeight;
			DistanceToEyes = distanceToEyes;
			DistanceBetweenEyes = distanceBetweenEyes;

			eyeLeftXPos = XResolution * PixelWidth / 2 - distanceBetweenEyes / 2;
			eyeRightXPos = XResolution * PixelWidth / 2 + distanceBetweenEyes / 2;
			eyeYPos = YResolution * PixelHeight / 2;

			nearestLeft = new int?[XResolutionInternal, YResolution];
			nearestRight = new int?[XResolutionInternal, YResolution];
			Clear();
		}
		
		public int XResolution { get; private set; }
		public int YResolution { get; private set; }
		public int SubpixelsCount { get; private set; }
		public double PixelWidth { get; private set; }
		public double PixelHeight { get; private set; }
		public double DistanceToEyes { get; private set; }
		public double DistanceBetweenEyes { get; private set; }

		public double Width
		{
			get
			{
				return XResolution * PixelWidth;
			}
		}

		public double Height
		{
			get
			{
				return YResolution * PixelHeight;
			}
		}

		public int XResolutionInternal
		{
			get
			{
				return XResolution * SubpixelsCount;
			}
		}
		
		public double PixelWidthInternal
		{
			get
			{
				return PixelWidth / SubpixelsCount;
			}
		}

		/// <summary>
		/// Width + GetAdditionalWidth(z) * 2 = width of visible area on range Z from the screen
		/// </summary>
		public double GetAdditionalWidth(double z)
		{
			return z / DistanceToEyes * eyeRightXPos;
		}

		/// <summary>
		/// Height + GetAdditionalHeight(z) * 2 = height of visible area on range Z from the screen
		/// </summary>
		public double GetAdditionalHeight(double z)
		{
			return z / DistanceToEyes * eyeYPos;
		}
		
		private double GetEyeXPos(bool isLeft)
		{
			return isLeft ? eyeLeftXPos : eyeRightXPos;
		}
		
		/// <summary>
		/// x coordinate of point projection on the screen
		/// </summary>
		private int GetXProj(double x, double z, bool isLeft)
		{
			double eyeXPos = GetEyeXPos(isLeft);
			return (int)Math.Floor(((x - eyeXPos) / (z + DistanceToEyes) * DistanceToEyes + eyeXPos) / PixelWidthInternal);
		}

		/// <summary>
		/// y coordinate of point projection on the screen
		/// </summary>
		private int GetYProj(double y, double z)
		{
			return (int)Math.Floor(((y - eyeYPos) / (z + DistanceToEyes) * DistanceToEyes + eyeYPos) / PixelHeight);
		}

		/// <summary>
		/// update coordinate of nearest point, which visible through given projection on the screen
		/// </summary>
		private void AdjustNearest(int xLeftProj, int xRightProj, int yProj)
		{
			if (yProj < 0 || yProj >= YResolution)
				return;
			
			if (xLeftProj >= 0 && xLeftProj < XResolutionInternal)
			{
				if (nearestLeft[xLeftProj, yProj] == null)
					nearestLeft[xLeftProj, yProj] = xRightProj;
				else
					nearestLeft[xLeftProj, yProj] = Math.Min(nearestLeft[xLeftProj, yProj].Value, xRightProj);
			}

			if (xRightProj >= 0 && xRightProj < XResolutionInternal)
			{
				if (nearestRight[xRightProj, yProj] == null)
					nearestRight[xRightProj, yProj] = xLeftProj;
				else
					nearestRight[xRightProj, yProj] = Math.Max(nearestRight[xRightProj, yProj].Value, xLeftProj);
			}
		}
		
		public void Clear()
		{
			for (int x = 0; x < XResolutionInternal; x++)
				for (int y = 0; y < YResolution; y++)
					nearestLeft[x, y] = nearestRight[x, y] = null;
		}
		
		public Rectangle2D<int> GetConvexHullProjBoundaries(Point3D[] points, bool isLeft)
		{
			int xProjMin = XResolutionInternal - 1, xProjMax = 0, yProjMin = YResolution - 1, yProjMax = 0;

			foreach (Point3D point in points)
			{
				int xProj = GetXProj(point.X, point.Z, isLeft);
				int yProj = GetYProj(point.Y, point.Z);

				xProjMin = Math.Min(xProj, xProjMin);
				xProjMax = Math.Max(xProj, xProjMax);
				yProjMin = Math.Min(yProj, yProjMin);
				yProjMax = Math.Max(yProj, yProjMax);
			}

			return new Rectangle2D<int> { X1 = Math.Max(xProjMin, 0), X2 = Math.Min(xProjMax, XResolutionInternal - 1), 
			 Y1 = Math.Max(yProjMin, 0), Y2 = Math.Min(yProjMax, YResolution - 1) };
		}

		/// <summary>
		/// calculate intersections of some figure with rays from eyes through each pixel
		/// </summary>
		public void AddRayTracedObject(Rectangle2D<int> leftProjBoundaries, Rectangle2D<int> rightProjBoundaries, GetRayIntersection getRayIntersection)
		{
			for (int eye = 0; eye <= 1; eye++)
			{
				Rectangle2D<int> projBoundaries = (eye == 0 ? leftProjBoundaries : rightProjBoundaries);
				double eyeXPos = GetEyeXPos(eye == 0);

				for (int xProj = projBoundaries.X1; xProj <= projBoundaries.X2; xProj++)
					for (int yProj = projBoundaries.Y1; yProj <= projBoundaries.Y2; yProj++)
					{
						double x, z;
						if (getRayIntersection(xProj, yProj, eyeXPos, out x, out z))
						{
							if (eye == 0)
								AdjustNearest(xProj, GetXProj(x, z, false), yProj);
							else
								AdjustNearest(GetXProj(x, z, true), xProj, yProj);
						}
					}
			}
		}

		public void AddRayTracedObject(double x1, double x2, double y1, double y2, double z1, double z2, GetRayIntersection getRayIntersection)
		{
			List<Point3D> points = new List<Point3D>(8);
			for (int x = 0; x < 2; x++)
				for (int y = 0; y < 2; y++)
					for (int z = 0; z < 2; z++)
						points.Add(new Point3D { X = (x == 0 ? x1 : x2), Y = (y == 0 ? y1 : y2), Z = (z == 0 ? z1 : z2) });

			AddRayTracedObject(points.ToArray(), getRayIntersection);
		}

		public void AddRayTracedObject(Point3D[] points, GetRayIntersection getRayIntersection)
		{
			AddRayTracedObject(GetConvexHullProjBoundaries(points, true), GetConvexHullProjBoundaries(points, false), getRayIntersection);
		}

		private bool GetSphereRayIntersection(double xCenter, double yCenter, double zCenter, double radius, int xProj, int yProj, double eyeXPos,
		 out double x, out double z)
		{
			x = z = 0;
			
			double kx = xProj * PixelWidthInternal - eyeXPos;
			double ky = yProj * PixelHeight - eyeYPos;
			double kz = DistanceToEyes;

			double x0 = xCenter - eyeXPos;
			double y0 = yCenter - eyeYPos;
			double z0 = zCenter + DistanceToEyes;

			double a = kx * kx + ky * ky + kz * kz;
			double b = -2 * (kx * x0 + ky * y0 + kz * z0);
			double c = x0 * x0 + y0 * y0 + z0 * z0 - radius * radius;

			double d = b * b - 4 * a * c;
			if (d < 0)
				return false;

			double t = (-b - Math.Sqrt(d)) / 2 / a;
			x = eyeXPos + kx * t;
			z = kz * t - DistanceToEyes;
			
			return true;
		}
		
		public void AddSphere(double xCenter, double yCenter, double zCenter, double radius)
		{
			AddRayTracedObject(xCenter - radius, xCenter + radius, yCenter - radius, yCenter + radius, zCenter - radius, zCenter + radius,
			 (int xProj, int yProj, double eyeXPos, out double x, out double z) => 
			 GetSphereRayIntersection(xCenter, yCenter, zCenter, radius, xProj, yProj, eyeXPos, out x, out z));
		}
		
		private bool GetCylinderRayIntersection(Point3D p1, Point3D p2, double radius, int xProj, int yProj, double eyeXPos, out double x, out double z)
		{
			x = z = 0;

			Point3D proj = new Point3D { X = xProj * PixelWidthInternal, Y = yProj * PixelHeight, Z = 0 };
			Point3D eye = new Point3D { X = eyeXPos, Y = eyeYPos, Z = -DistanceToEyes };

			Vector3D v1 = (proj - eye).Normalize(), v2 = (p2 - p1).Normalize();
			double i = Vector3D.ScalarProduct(v1, v2), j = Vector3D.ScalarProduct(eye - p1, v2);
			Vector3D vAux1 = v1 - v2 * i, vAux2 = eye - p1 - v2 * j;
			double a = vAux1.Abs2(), b = 2 * Vector3D.ScalarProduct(vAux1, vAux2), c = vAux2.Abs2() - radius * radius;
			double discr = b * b - 4 * a * c;

			if (discr < 0 || a < 1e-8)
				return false;

			double x1 = Math.Min(p1.X, p2.X) - 1e-8, x2 = Math.Max(p1.X, p2.X) + 1e-8;
			double y1 = Math.Min(p1.Y, p2.Y) - 1e-8, y2 = Math.Max(p1.Y, p2.Y) + 1e-8;
			double z1 = Math.Min(p1.Z, p2.Z) - 1e-8, z2 = Math.Max(p1.Z, p2.Z) + 1e-8;

			for (int r = 0; r <= 1; r++)
			{
				double t = (-b + (r == 0 ? -Math.Sqrt(discr) : Math.Sqrt(discr))) / 2 / a;
				Point3D pAxis = p1 + v2 * (i * t + j);

				if (pAxis.X >= x1 && pAxis.X <= x2 && pAxis.Y >= y1 && pAxis.Y <= y2 && pAxis.Z >= z1 && pAxis.Z <= z2)
				{
					x = eyeXPos + v1.X * t;
					z = -DistanceToEyes + v1.Z * t;
					return true;
				}
			}

			return false;
		}

		public void AddCylinder(Point3D p1, Point3D p2, double radius)
		{
			Vector3D v1 = (p1 - p2).GetSomeOrtVector().Normalize() * Math.Sqrt(2);
			Vector3D v2 = Vector3D.CrossProduct(v1, p1 - p2).Normalize() * Math.Sqrt(2);
			Point3D[] points = new Point3D[] { p1 - v1 * radius, p1 + v1 * radius, p1 - v2 * radius, p1 + v2 * radius, 
			 p2 - v1 * radius, p2 + v1 * radius, p2 - v2 * radius, p2 + v2 * radius };

			AddRayTracedObject(points, 
			 (int xProj, int yProj, double eyeXPos, out double x, out double z) => GetCylinderRayIntersection(p1, p2, radius, xProj, yProj, eyeXPos, out x, out z));
		}

		private bool GetPlaneRayIntersection(int xProj, int yProj, double eyeXPos, double a, double b, double c, double d, out double x, out double z, CheckPoint checkPoint)
		{
			x = z = 0;
			
			double dx = xProj * PixelWidthInternal - eyeXPos;
			double dy = yProj * PixelHeight - eyeYPos;
			double dz = DistanceToEyes;
			
			double t;
			if (!MathHelper.Divide(-(a * eyeXPos + b * eyeYPos - c * DistanceToEyes + d), a * dx + b * dy + c * dz, out t))
				return false;

			x = eyeXPos + t * dx;
			z = -DistanceToEyes + t * dz;
			
			return checkPoint != null ? checkPoint(x, eyeYPos + t * dy, z) : true;
		}

		public void AddPlane(double x1, double x2, double y1, double y2, double z1, double z2, double a, double b, double c, double d, CheckPoint checkPoint = null)
		{
			AddRayTracedObject(x1, x2, y1, y2, z1, z2,
			 (int xProj, int yProj, double eyeXPos, out double x, out double z) => 
			 GetPlaneRayIntersection(xProj, yProj, eyeXPos, a, b, c, d, out x, out z, checkPoint));
		}

		public void AddPlane(double x1, double x2, double y1, double y2, double z1, double z2, double a, double b, double c, double x0, double y0, double z0, 
		 CheckPoint checkPoint = null)
		{
			AddPlane(x1, x2, y1, y2, z1, z2, a, b, c, -a * x0 - b * y0 - c * z0, checkPoint);
		}

		public void AddHorizontalPlane(double x1, double x2, double y1, double y2, double z, CheckPoint checkPoint = null)
		{
			AddPlane(x1, x2, y1, y2, z, z, 0, 0, 1, 0, 0, z, checkPoint);
		}

		public void AddHorizontalPlane(double z, CheckPoint checkPoint = null)
		{
			double addWidth = GetAdditionalWidth(z), addHeight = GetAdditionalHeight(z);
			AddHorizontalPlane(-addWidth, Width + addWidth, -addHeight, Height + addHeight, z, checkPoint);
		}
		
		public void AddPolygon(Point3D[] points)
		{
			if (points.Length < 3)
				throw new ArgumentException("points.Length < 3");

			Point3D innerPoint = new Point3D { X = points.Average(p => p.X), Y = points.Average(p => p.Y), Z = points.Average(p => p.Z) };
			Vector3D ortVec = Vector3D.CrossProduct(innerPoint - points[0], points[1] - points[0]);
			Plane plane = new Plane(ortVec, points[0]);
			
			Plane[] ortPlanes = new Plane[points.Length];

			for (int n = 0; n < points.Length; n++)
			{
				ortPlanes[n] = new Plane(Vector3D.CrossProduct(points[(n + 1) % points.Length] - points[n], ortVec), points[n]);
				ortPlanes[n].Normalize();
			}

			AddRayTracedObject(points, (int xProj, int yProj, double eyeXPos, out double x, out double z) =>
			 GetPlaneRayIntersection(xProj, yProj, eyeXPos, plane.A, plane.B, plane.C, plane.D, out x, out z, (x1, y1, z1) => ortPlanes.All(p => p.A * x1 + p.B * y1 + p.C * z1 + p.D > -1e-8)));
		}

		public void AddTriangle(double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3)
		{
			AddPolygon(new Point3D[] { new Point3D { X = x1, Y = y1, Z = z1 }, new Point3D { X = x2, Y = y2, Z = z2 }, new Point3D { X = x3, Y = y3, Z = z3 } });
		}

		private static bool IsDivisionRequired(int x1, int y1, int x2, int y2, int x3, int y3, int x4, int y4)
		{
			if (Math.Abs(x1 - x2) > 1 || Math.Abs(y1 - y2) > 1 || Math.Abs(x2 - x3) > 1 || Math.Abs(y2 - y3) > 1 ||
			 Math.Abs(x3 - x4) > 1 || Math.Abs(y3 - y4) > 1 || Math.Abs(x4 - x1) > 1 || Math.Abs(y4 - y1) > 1)
				return true;

			if ((Math.Abs(x1 - x3) == 0 && Math.Abs(y1 - y3) == 2 && Math.Abs(x2 - x4) == 2 && Math.Abs(y2 - y4) == 0) ||
			 (Math.Abs(x1 - x3) == 2 && Math.Abs(y1 - y3) == 0 && Math.Abs(x2 - x4) == 0 && Math.Abs(y2 - y4) == 2))
				return true;

			return false;
		}

		private bool IsDivisionRequired(SurfaceRenderStackFrame stackFrame, SurfaceFunc func)
		{
			// if some points are out of func domain (when grid cell located on domain edge), then use interpolated values

			List<double> zValues = new List<double>(3);

			for (int x1 = 0; x1 <= 2; x1 += 2)
				for (int y1 = 0; y1 <= 2; y1 += 2)
				{
					if (stackFrame.Points[x1, y1].Z != null)
						continue;

					zValues.Clear();

					for (int x2 = 0; x2 <= 2; x2 += 2)
						for (int y2 = 0; y2 <= 2; y2 += 2)
						{
							if (stackFrame.Points[x2, y2].Z == null)
								continue;
							
							double? z = func(stackFrame.Points[x2, y2].X + (stackFrame.Points[x2, y2].X - stackFrame.Points[x1, y1].X) / 100, 
							 stackFrame.Points[x2, y2].Y + (stackFrame.Points[x2, y2].Y - stackFrame.Points[x1, y1].Y) / 100);

							if (z != null)
								zValues.Add(stackFrame.Points[x2, y2].Z.Value + (stackFrame.Points[x2, y2].Z.Value - z.Value) * 100);
						}

					if (zValues.Count == 0)
						return false;

					double zAvg = zValues.Average();

					stackFrame.Points[x1, y1].XLeftProj = GetXProj(stackFrame.Points[x1, y1].X, zAvg, true);
					stackFrame.Points[x1, y1].XRightProj = GetXProj(stackFrame.Points[x1, y1].X, zAvg, false);
					stackFrame.Points[x1, y1].YProj = GetYProj(stackFrame.Points[x1, y1].Y, zAvg);
				}

			// if adjacent points of grid with same x or y coordinate are projected to non-adjacent points on the screen, then divide grid cell on subcells

			return IsDivisionRequired(stackFrame.Points[0, 0].XLeftProj, stackFrame.Points[0, 0].YProj, stackFrame.Points[2, 0].XLeftProj, stackFrame.Points[2, 0].YProj, 
			 stackFrame.Points[2, 2].XLeftProj, stackFrame.Points[2, 2].YProj, stackFrame.Points[0, 2].XLeftProj, stackFrame.Points[0, 2].YProj) ||
			 IsDivisionRequired(stackFrame.Points[0, 0].XRightProj, stackFrame.Points[0, 0].YProj, stackFrame.Points[2, 0].XRightProj, stackFrame.Points[2, 0].YProj, 
			 stackFrame.Points[2, 2].XRightProj, stackFrame.Points[2, 2].YProj, stackFrame.Points[0, 2].XRightProj, stackFrame.Points[0, 2].YProj);
		}

		/// <summary>
		/// calculate point projection for some cell and, if required, divide that cell on subcells and repeat procedure for them
		/// </summary>
		private void DivideCell(SurfaceRenderStackFrame[] stackFrames, SurfaceFunc func)
		{
			int depth = 0;
			stackFrames[0].Step = 0;

			while (depth >= 0)
			{
				var currentFrame = stackFrames[depth];

				if (currentFrame.Step == 0)
				{
					if (depth == stackFrames.Length - 1 || !IsDivisionRequired(currentFrame, func))  // no division required or division limit reached?
					{
						// save calculated point projections

						for (int x = 0; x <= 2; x++)
							for (int y = 0; y <= 2; y++)
								if (currentFrame.Points[x, y].Z != null)
									AdjustNearest(currentFrame.Points[x, y].XLeftProj, currentFrame.Points[x, y].XRightProj, currentFrame.Points[x, y].YProj);

						depth--;
					}
					else
					{
						// calculate additional point coordinates and projections (these values will be used for grid subcells)

						double[] xCoords = { currentFrame.Points[0, 0].X, (currentFrame.Points[0, 0].X + currentFrame.Points[2, 0].X) / 2, currentFrame.Points[2, 0].X };
						double[] yCoords = { currentFrame.Points[0, 0].Y, (currentFrame.Points[0, 0].Y + currentFrame.Points[0, 2].Y) / 2, currentFrame.Points[0, 2].Y };

						for (int x = 0; x <= 2; x++)
							for (int y = 0; y <= 2; y++)
								if (x == 1 || y == 1)
									InitPoint(ref currentFrame.Points[x, y], xCoords[x], yCoords[y], func(xCoords[x], yCoords[y]));

						currentFrame.Step++;
					}
				}
				else if (currentFrame.Step >= 1 && currentFrame.Step <= 4)
				{
					// initialize next grid subcell and pass to processing of it

					var nextFrame = stackFrames[depth + 1];
					nextFrame.Step = 0;

					for (int x = 0; x <= 1; x++)
						for (int y = 0; y <= 1; y++)
							nextFrame.Points[x * 2, y * 2] = currentFrame.Points[(currentFrame.Step - 1) / 2 + x, (currentFrame.Step - 1) % 2 + y];

					currentFrame.Step++;
					depth++;
				}
				else
					depth--;
			}
		}

		private void InitPoint(ref SurfacePoint point, double x, double y, double? z)
		{
			point.X = x;
			point.Y = y;
			point.Z = z;

			if (z != null)
			{
				point.XLeftProj = GetXProj(x, z.Value, true);
				point.XRightProj = GetXProj(x, z.Value, false);
				point.YProj = GetYProj(y, z.Value);
			}
		}

		/// <summary>
		/// create figure using surface function
		/// </summary>
		public void AddSurface(SurfaceFunc func, double x1, double x2, double y1, double y2, double zClosest = 0)
		{
			// calculate surface function values for each point of coordinate grid, and divide grid cell on subcells, if required (for more accurate calculations)

			int xCellCount = (int)Math.Ceiling((x2 - x1) / (PixelWidthInternal / DistanceToEyes * (DistanceToEyes + zClosest)));
			int yCellCount = (int)Math.Ceiling((y2 - y1) / (PixelHeight / DistanceToEyes * (DistanceToEyes + zClosest)));
			
			var stackFrames = new SurfaceRenderStackFrame[10];
			for (int n = 0; n < stackFrames.Length; n++)
				stackFrames[n] = new SurfaceRenderStackFrame();

			double?[] lineZValues = new double?[xCellCount + 1];
			double yPrev = 0;

			for (int yCell = 0; yCell <= yCellCount; yCell++)
			{
				double y = (yCell != yCellCount ? y1 + yCell * (y2 - y1) / yCellCount : y2);

				for (int xCell = 0; xCell <= xCellCount; xCell++)
				{
					double x = (xCell != xCellCount ? x1 + xCell * (x2 - x1) / xCellCount : x2);

					double? z = func(x, y);

					if (yCell != 0)
					{
						var firstFrame = stackFrames[0];

						if (xCell <= 1)
						{
							InitPoint(ref firstFrame.Points[xCell == 0 ? 0 : 2, 0], x, yPrev, lineZValues[xCell]);
							InitPoint(ref firstFrame.Points[xCell == 0 ? 0 : 2, 2], x, y, z);
						}
						else
						{
							firstFrame.Points[0, 0] = firstFrame.Points[2, 0];
							firstFrame.Points[0, 2] = firstFrame.Points[2, 2];

							InitPoint(ref firstFrame.Points[2, 0], x, yPrev, lineZValues[xCell]);
							InitPoint(ref firstFrame.Points[2, 2], x, y, z);
						}

						if (xCell >= 1)
							DivideCell(stackFrames, func);
					}

					lineZValues[xCell] = z;
				}

				yPrev = y;
			}
		}

		/// <param name="zFarthest">determines bounds of the grid (grid must coincide with visible region)</param>
		/// <param name="zClosest">determines base grid cell size</param>
		public void AddSurface(SurfaceFunc func, double zFarthest, double zClosest = 0)
		{
			double addWidth = GetAdditionalWidth(zFarthest), addHeight = GetAdditionalHeight(zFarthest);
			AddSurface(func, -addWidth, Width + addWidth, -addHeight, Height + addHeight, zClosest);
		}

		public void AddModelByImage(Bitmap image, Point3D origin, Vector3D xVec, Vector3D yVec, Vector3D zVec, double xSize, double ySize, double zSize, 
		 bool whiteIsBackground = true)
		{
			xVec.Normalize();
			yVec.Normalize();
			zVec.Normalize();

			Func<int, int, Point3D> getPoint = ((x, y) =>
				{
					Color color = image.GetPixel(x, y);
					int depth = (whiteIsBackground ? color.R + color.G + color.B : 765 - (color.R + color.G + color.B));

					if (depth != 765)
						return origin + xVec * xSize * x / (image.Width - 1) + yVec * ySize * y / (image.Height - 1) + zVec * zSize * depth / 765;
					else
						return null;
				});

			for (int x = 0; x < image.Width - 1; x++)
				for (int y = 0; y < image.Height - 1; y++)
				{
					Point3D p1 = getPoint(x, y), p2 = getPoint(x + 1, y), p3 = getPoint(x, y + 1), p4 = getPoint(x + 1, y + 1);

					if (p1 != null && p2 != null && p3 != null)
						AddPolygon(new Point3D[] { p1, p2, p3 });
					if (p2 != null && p3 != null && p4 != null)
						AddPolygon(new Point3D[] { p2, p3, p4 });
					if (p2 == null && p1 != null && p3 != null && p4 != null)
						AddPolygon(new Point3D[] { p1, p3, p4 });
					if (p3 == null && p1 != null && p2 != null && p4 != null)
						AddPolygon(new Point3D[] { p1, p2, p4 });
				}
		}

		private void DefaultColorGenerator(Bitmap bitmap, int y, int uniqueColorsCount, int[] colors)
		{
			int?[] bitmapColors = new int?[uniqueColorsCount];

			for (int x = 0; x < XResolutionInternal; x++)
			{
				if (bitmapColors[colors[x]] == null)
				{
					bool colorSet = false;

					if (x > 0 && x < XResolutionInternal - 1)
					{
						int n;
						for (n = 1; n <= SubpixelsCount && x + n < XResolutionInternal; n++)
							if (bitmapColors[colors[x + n]] != null)
							{
								colorSet = true;
								break;
							}

						if (colorSet)
						{
							int color1 = bitmapColors[colors[x - 1]].Value, color2 = bitmapColors[colors[x + n]].Value;

							for (int m = 0; m < n; m++)
								bitmapColors[colors[x + m]] = (int)Math.Round((double)(color1 * (n - m) + color2 * (m + 1)) / (n + 1));
						}
					}

					if (!colorSet)
					{
						int randMax = 115 + SubpixelsCount * 35;
						int c = ((x / SubpixelsCount > 0 ? bitmap.GetPixel(x / SubpixelsCount - 1, y).R : 128) +
						 (y > 0 ? bitmap.GetPixel(x / SubpixelsCount, y - 1).R : 128)) / 2 + random.Next(randMax) - randMax / 2 + 1;
						bitmapColors[colors[x]] = Math.Min(Math.Max(c, 0), 255);
					}
				}

				if (x % SubpixelsCount == SubpixelsCount - 1)
				{
					int cAvg = 0;
					for (int n = 0; n < SubpixelsCount; n++)
						cAvg += bitmapColors[colors[x - n]].Value;
					cAvg = (int)Math.Round((double)cAvg / SubpixelsCount);

					cAvg = Math.Max(Math.Min((int)((cAvg - 128) * (0.95 + SubpixelsCount * 0.05) + 128), 255), 0);

					bitmap.SetPixel(x / SubpixelsCount, y, Color.FromArgb(cAvg, cAvg, cAvg));
				}
			}
		}

		public Bitmap GenerateBitmap(ColorGenerator colorGenerator = null)
		{
			colorGenerator = colorGenerator ?? DefaultColorGenerator;

			Bitmap bitmap = new Bitmap(XResolution, YResolution);
			int[] colors = new int[XResolutionInternal];
			
			for (int y = 0; y < YResolution; y++)
			{
				int currentColor = 0, xBeg = -1;
				
				for (int n = 0; n < colors.Length; n++)
					colors[n] = -1;
					
				while (true)
				{
					do
						xBeg++;
					while (xBeg < XResolutionInternal && colors[xBeg] != -1);
						
					if (xBeg == XResolutionInternal)
						break;

					// enumerate projection points from left to right

					bool hasCycle = false;
					int x = xBeg;

					while (true)
					{
						colors[x] = currentColor;
						int? xNext = nearestLeft[x, y];
						if (xNext == null || xNext < 0 || xNext >= XResolutionInternal || nearestRight[xNext.Value, y] != x)
							break;

						x = xNext.Value;

						if (x == xBeg)
						{
							hasCycle = true;
							break;
						}
					}
					
					if (!hasCycle)
					{
						// enumerate projection points from right to left

						x = xBeg;

						while (true)
						{
							colors[x] = currentColor;
							int? xNext = nearestRight[x, y];
							if (xNext == null || xNext < 0 || xNext >= XResolutionInternal || nearestLeft[xNext.Value, y] != x)
								break;

							x = xNext.Value;
						}
					}

					currentColor++;
				}
				
				colorGenerator(bitmap, y, currentColor, colors);
			}
			
			return bitmap;
		}

		private void SampleImageColorGenerator(Bitmap bitmap, Bitmap backgroundImage, int y, int uniqueColorsCount, int[] colors)
		{
			Color?[] bitmapColors = new Color?[uniqueColorsCount];
			int nextXWithColor = -1;

			for (int x = 0; x < colors.Length; x++)
			{
				if (bitmapColors[colors[x]] == null)
				{
					if (nextXWithColor < x)
					{
						for (int n = 1; n <= SubpixelsCount * 5; n++)
							if (x + n < XResolutionInternal && bitmapColors[colors[x + n]] != null)
							{
								nextXWithColor = x + n;
								break;
							}
					}

					Color currentColor;
					if (y % backgroundImage.Height < 5 && y / backgroundImage.Height != 0)
					{
						Color c1 = backgroundImage.GetPixel((x / SubpixelsCount) % backgroundImage.Width, backgroundImage.Height - 1);
						Color c2 = backgroundImage.GetPixel((x / SubpixelsCount) % backgroundImage.Width, y % backgroundImage.Height);
						double q = (y % backgroundImage.Height + 1) / 6.0;
						currentColor = Color.FromArgb((int)Math.Round(c1.R * (1 - q) + c2.R * q), (int)Math.Round(c1.G * (1 - q) + c2.G * q), (int)Math.Round(c1.B * (1 - q) + c2.B * q));
					}
					else
						currentColor = backgroundImage.GetPixel((x / SubpixelsCount) % backgroundImage.Width, y % backgroundImage.Height);

					if (nextXWithColor > x)
					{
						Color nextColor = bitmapColors[colors[nextXWithColor]].Value;
						double q = (double)(nextXWithColor - x) / (SubpixelsCount * 5 + 1);
						bitmapColors[colors[x]] = Color.FromArgb((int)Math.Round(currentColor.R * q + nextColor.R * (1 - q)), (int)Math.Round(currentColor.G * q + nextColor.G * (1 - q)),
						 (int)Math.Round(currentColor.B * q + nextColor.B * (1 - q)));
					}
					else
						bitmapColors[colors[x]] = currentColor;
				}

				if (x % SubpixelsCount == SubpixelsCount - 1)
				{
					int r = 0, g = 0, b = 0;
					for (int n = 0; n < SubpixelsCount; n++)
					{
						r += bitmapColors[colors[x - x % SubpixelsCount + n]].Value.R;
						g += bitmapColors[colors[x - x % SubpixelsCount + n]].Value.G;
						b += bitmapColors[colors[x - x % SubpixelsCount + n]].Value.B;
					}

					Color color = Color.FromArgb((int)Math.Round((double)r / SubpixelsCount), (int)Math.Round((double)g / SubpixelsCount), (int)Math.Round((double)b / SubpixelsCount));
					bitmap.SetPixel(x / SubpixelsCount, y, color);
				}
			}
		}

		public Bitmap GenerateBitmap(Bitmap backgroundImage)
		{
			return GenerateBitmap((bitmap, y, uniqueColorsCount, colors) => SampleImageColorGenerator(bitmap, backgroundImage, y, uniqueColorsCount, colors));
		}
	}
}

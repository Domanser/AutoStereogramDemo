using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace AutoStereogramDemo
{
	public class MathHelper
	{
		public static bool Divide(double num, double denom, out double x)
		{
			x = 0;

			if (Math.Abs(denom) < Math.Abs(num))
			{
				if (Math.Abs(denom / num) < 1e-7)
					return false;
			}
			else if (denom == 0 && num == 0)
				return false;
				
			x = num / denom;
			return true;
		}
	}
}

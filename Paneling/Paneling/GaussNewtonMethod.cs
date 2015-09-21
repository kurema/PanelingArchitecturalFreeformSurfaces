using Rhino;
using Rhino.Geometry;
using Rhino.DocObjects;
using Rhino.Collections;

using System;
using System.IO;
using System.Xml;
using System.Xml.Linq;
using System.Linq;
using System.Data;
using System.Drawing;
using System.Reflection;
using System.Collections;
using System.Windows.Forms;
using System.Collections.Generic;
using System.Runtime.InteropServices;

namespace Paneling
{
    public static class GaussNewtonMethod
    {
        static public Matrix SimpleSingle(Matrix Jacobian, Matrix EnergyVector, Matrix StatusVector,double ZeroTolerance=1e-10)
        {
            Matrix M1 = Jacobian.Duplicate();
            M1.Transpose();
            Matrix M2 = M1 * Jacobian;
            M2.Invert(ZeroTolerance);
            M1 = M1 * EnergyVector;
            M1 = M2 * M1;
            M1.Scale(-1);
            return M1;
        }
    }
}

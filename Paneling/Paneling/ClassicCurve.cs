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
    class ClassicCurve
    {
        public class CurveOptimization
        {
            public List<Mold> MoldList = new List<Mold>();
            public Rhino.Geometry.Curve TargetCurve;
            public int EvaluationCount = 300;
            public int ContinuousStepCount = 100;
            public int Strokes = 5;

            public MessageManager Messager = new MessageManager();

            public List<Rhino.Geometry.Curve> CurveHistory = new List<Rhino.Geometry.Curve>();

            public void Solve()
            {
                Messager.Show();
                int count = Strokes;
                for (int i = 0; i < count; i++)
                {
                    Vector3d vector = TargetCurve.PointAtNormalizedLength((double)(i + 1) / (double)count) - TargetCurve.PointAtNormalizedLength((double)i / (double)count);
                    var m = new Molds.Arc();
                    m.Length = vector.Length;
                    m.Angle = Math.Acos(vector.X / vector.Length) * (vector.Y > 0 ? 1 : -1);
                    MoldList.Add(m);
                }
                ContinuousStep(ContinuousStepCount);
            }

            public Rhino.Geometry.Matrix GaussNewton(Rhino.Geometry.Matrix Jacobian, Rhino.Geometry.Matrix EnergyVector, Rhino.Geometry.Matrix StatusVector)
            {
                //Jacobian.Transpose();
                Rhino.Geometry.Matrix M1 = Jacobian.Duplicate();
                M1.Transpose();
                Rhino.Geometry.Matrix M2 = M1 * Jacobian;
                M2.Invert(1e-10);//(ZeroTolerance)
                M1 = M1 * EnergyVector;
                M1 = M2 * M1;
                M1.Scale(-1);
                return M1;
            }

            int step = 0;
            public void ContinuousStep(int count)
            {
                var maxdiff = GetMaxStatusDiff();
                var minimumStatus = GetStatus();
                double minimumValue = double.MaxValue;
                int failureCount = 0;

                string outcsv = "i,minimumValue,sumenergy,failureCount,result\r\n";

                for (int i = 0; i < count; i++)
                {
                    CurveHistory.Add(GetCurve());

                    step = i;
                    Messager.SetProgress(i * 100 / count);
                    var status = GetStatus();
                    var energy = GetEnergy();
                    double sumenergy = 0;
                    for (int j = 0; j < energy.RowCount; j++)
                    {
                        sumenergy += Math.Pow(energy[j, 0], 2);
                    }
                    var origStatus = status.Duplicate();
                    var statusDif = GaussNewton(GetJacobian(), energy, status);

                    double baseValue = 2.0;
                    double MinEnergySum = sumenergy;
                    double MinStatusDifScale = 1.0;

                    for (int j = 1; j < 10; j++)
                    {
                        double tmpES = 0;
                        var tmpStatusDif = statusDif.Duplicate();
                        tmpStatusDif.Scale(Math.Pow(baseValue, j));
                        SetStatus(status + tmpStatusDif);
                        if ((tmpES = this.GetEnergySum()) < MinEnergySum)
                        {
                            MinEnergySum = tmpES;
                            MinStatusDifScale = Math.Pow(baseValue, j);
                        }

                        tmpStatusDif = statusDif.Duplicate();
                        tmpStatusDif.Scale(Math.Pow(baseValue, -j));
                        SetStatus(status + tmpStatusDif);
                        if ((tmpES = this.GetEnergySum()) < MinEnergySum)
                        {
                            MinEnergySum = tmpES;
                            MinStatusDifScale = Math.Pow(baseValue, -j);
                        }
                    }
                    for (int j = 0; j < 10; j++)
                    {
                        baseValue = Math.Sqrt(baseValue);

                        var tmpStatusDif = statusDif.Duplicate();
                        tmpStatusDif.Scale(MinStatusDifScale * baseValue);
                        SetStatus(status + tmpStatusDif);
                        double tmpES = 0;
                        if ((tmpES = this.GetEnergySum()) < MinEnergySum)
                        {
                            MinEnergySum = tmpES;
                            MinStatusDifScale = MinStatusDifScale * baseValue;
                        }

                        tmpStatusDif = statusDif.Duplicate();
                        tmpStatusDif.Scale(MinStatusDifScale / baseValue);
                        SetStatus(status + tmpStatusDif);
                        tmpES = 0;
                        if ((tmpES = this.GetEnergySum()) < MinEnergySum)
                        {
                            MinEnergySum = tmpES;
                            MinStatusDifScale = MinStatusDifScale / baseValue;
                        }
                    }

                    if (MinStatusDifScale == 1.0) { SetStatus(status); break; }
                    statusDif.Scale(MinStatusDifScale);
                    Messager.WriteLine("DifScale: " + MinStatusDifScale);
                    SetStatus(status + statusDif);
                }
                Messager.WriteLine(outcsv);
                Messager.SetProgress(100);
            }

            public double GetEnergySum()
            {
                var energy = GetEnergy();
                double sumenergy = 0;
                for (int j = 0; j < energy.RowCount; j++)
                {
                    sumenergy += Math.Pow(energy[j, 0], 2);
                }
                return sumenergy;
            }

            public Rhino.Geometry.Matrix GetJacobian()
            {
                Rhino.Geometry.Matrix status = GetStatus();
                Rhino.Geometry.Matrix originEnergy = GetEnergy();
                Rhino.Geometry.Matrix ret = new Rhino.Geometry.Matrix(originEnergy.RowCount, status.RowCount);

                //string delme = "";
                for (int i = 0; i < status.RowCount; i++)
                {
                    double diff = Math.Max(status[i, 0] * 1e-6, 1e-6);
                    status[i, 0] += diff;
                    SetStatus(status);
                    Rhino.Geometry.Matrix newEnergy = GetEnergy();
                    for (int j = 0; j < originEnergy.RowCount; j++)
                    {
                        ret[j, i] = (newEnergy[j, 0] - originEnergy[j, 0]) / diff;
                        //delme += ret[j, i] + "\t";
                    }
                    //delme += "\r\n";
                    status[i, 0] -= diff;
                }
                //Messager.WriteLine(delme);
                SetStatus(status);
                return ret;
            }

            public Rhino.Geometry.Matrix GetEnergy()
            {
                Rhino.Geometry.Curve curve = GetCurve();
                //      Rhino.Geometry.Matrix ret = new Rhino.Geometry.Matrix(EvaluationCount + MoldList.Count(), 1);
                Rhino.Geometry.Matrix ret = new Rhino.Geometry.Matrix(EvaluationCount * 2, 1);
                //Rhino.Geometry.Matrix ret = new Rhino.Geometry.Matrix(EvaluationCount * 2 + MoldList.Count() - 1, 1);

                for (int i = 1; i <= EvaluationCount; i++)
                {
                    double t;
                    Point3d point = TargetCurve.PointAtNormalizedLength(((double)i) / ((double)EvaluationCount));
                    curve.ClosestPoint(point, out t);
                    ret[i - 1, 0] = point.DistanceTo(curve.PointAt(t)) * i;
                }
                for (int i = 1; i <= EvaluationCount; i++)
                {
                    double t;
                    Point3d point = curve.PointAtNormalizedLength(((double)i) / ((double)EvaluationCount));
                    TargetCurve.ClosestPoint(point, out t);
                    ret[EvaluationCount + i - 1, 0] = point.DistanceTo(TargetCurve.PointAt(t)) * i;
                }/*
      var curves = GetCurves();
      for(int i = 0;i < MoldList.Count() - 1;i++){
      double t;
      ret[EvaluationCount * 2 + i, 0] = (curves[i].TangentAt(curves[i].Domain.Max) - curves[i + 1].TangentAt(curves[i].Domain.Min)).Length * ((double) EvaluationCount / (double) MoldList.Count()) * 35;
      }*/

                /*
                int cnt = EvaluationCount;
                foreach(Mold m in MoldList){
                ret[cnt, 0] = m.GetCost() * 2;
                cnt++;
                }
                */
                return ret;
            }

            public Rhino.Geometry.Matrix GetMaxStatusDiff()
            {
                List<double> ret = new List<double>();
                foreach (Mold m in MoldList)
                {
                    ret.AddRange(m.GetMaxParamaters());
                }
                Rhino.Geometry.Matrix retmatrix = new Rhino.Geometry.Matrix(ret.Count(), 1);
                for (int i = 0; i < ret.Count(); i++)
                {
                    retmatrix[i, 0] = ret[i];
                }
                return retmatrix;
            }

            public Rhino.Geometry.Matrix GetStatus()
            {
                List<double> ret = new List<double>();
                foreach (Mold m in MoldList)
                {
                    ret.AddRange(m.GetParamaters());
                }
                Rhino.Geometry.Matrix retmatrix = new Rhino.Geometry.Matrix(ret.Count(), 1);
                for (int i = 0; i < ret.Count(); i++)
                {
                    retmatrix[i, 0] = ret[i];
                }
                return retmatrix;
            }

            public void SetStatus(Rhino.Geometry.Matrix status)
            {
                int cnt = 0;
                foreach (Mold m in MoldList)
                {
                    List<double> setmatrix = new List<double>();
                    for (int i = 0; i < m.GetParamaterCount(); i++)
                    {
                        setmatrix.Add(status[cnt, 0]);
                        cnt++;
                    }
                    m.SetParamaters(setmatrix);
                }
            }

            public Rhino.Geometry.Curve GetCurve()
            {
                Point3d CurrentPoint = TargetCurve.PointAtNormalizedLength(0);
                List<Rhino.Geometry.Curve> curves = new List<Rhino.Geometry.Curve>();
                foreach (Mold m in MoldList)
                {
                    Rhino.Geometry.Curve curve = m.GetCurve(CurrentPoint);
                    CurrentPoint = curve.PointAtEnd;
                    curves.Add(curve);
                }
                //if((Rhino.Geometry.Curve.JoinCurves(curves, 0.1)).Count() > 1){Messager.WriteLine("wow " + (Rhino.Geometry.Curve.JoinCurves(curves)).Count() + "\r\nstep" + step);throw new Exception();}
                return (Rhino.Geometry.Curve.JoinCurves(curves, 0.1))[0];
            }

            public Rhino.Geometry.Curve[] GetCurves()
            {
                Point3d CurrentPoint = TargetCurve.PointAtNormalizedLength(0);
                List<Rhino.Geometry.Curve> curves = new List<Rhino.Geometry.Curve>();
                foreach (Mold m in MoldList)
                {
                    Rhino.Geometry.Curve curve = m.GetCurve(CurrentPoint);
                    CurrentPoint = curve.PointAtEnd;
                    curves.Add(curve);
                }
                //      return (Rhino.Geometry.Curve.JoinCurves(curves));
                return curves.ToArray();
            }

            public interface Mold
            {
                List<double> GetParamaters();
                List<double> GetMaxParamaters();
                void SetParamaters(List<double> matrix);
                int GetParamaterCount();
                Rhino.Geometry.Curve GetCurve(Point3d Start);
                double GetCost();
            }

            public class Molds
            {
                public class Line : Mold
                {
                    public double Length = 1.0;
                    public double Angle = 0.0;

                    public List<double> GetParamaters()
                    {
                        List<double> ret = new List<double>();
                        ret.Add(this.Length);
                        ret.Add(this.Angle);
                        return ret;
                    }

                    public List<double> GetMaxParamaters()
                    {
                        List<double> ret = new List<double>();
                        ret.Add(0.1);
                        ret.Add(0.02);
                        return ret;
                    }

                    public void SetParamaters(List<double> matrix)
                    {
                        this.Length = matrix[0];
                        this.Angle = matrix[1];
                    }

                    public int GetParamaterCount()
                    {
                        return 2;
                    }

                    public Rhino.Geometry.Curve GetCurve(Point3d Start)
                    {
                        Rhino.Geometry.Curve ret = new Rhino.Geometry.LineCurve(Start, Start + new Vector3d(Length * Math.Cos(Angle), Length * Math.Sin(Angle), 0));
                        return ret;
                    }

                    public double GetCost() { return 1.0 * (GetCurve(new Point3d(0, 0, 0)).GetLength()); }
                }
                public class Arc : Mold
                {
                    public double Length = 1.0;
                    public double Angle = 0.0;
                    public double InteriorLength = 0.0;

                    public List<double> GetParamaters()
                    {
                        List<double> ret = new List<double>();
                        ret.Add(this.Length);
                        ret.Add(this.Angle);
                        ret.Add(this.InteriorLength);
                        return ret;
                    }

                    public List<double> GetMaxParamaters()
                    {
                        List<double> ret = new List<double>();

                        /*
                        ret.Add(0.3);
                        ret.Add(0.1);
                        ret.Add(0.3);
                        */

                        /*
                        ret.Add(0.1);
                        ret.Add(0.03);
                        ret.Add(0.7);
                        */
                        ret.Add(1);
                        ret.Add(0.05);
                        ret.Add(2);
                        /*
                        ret.Add(100);
                        ret.Add(100);
                        ret.Add(100);
                      */
                        return ret;
                    }

                    public void SetParamaters(List<double> matrix)
                    {
                        this.Length = matrix[0];
                        this.Angle = matrix[1];
                        this.InteriorLength = matrix[2];
                    }

                    public int GetParamaterCount()
                    {
                        return 3;
                    }

                    public Rhino.Geometry.Curve GetCurve(Point3d Start)
                    {
                        Point3d EndPoint = Start + new Vector3d(Length * Math.Cos(Angle), Length * Math.Sin(Angle), 0);
                        Point3d InteriorPoint = Start + new Vector3d(Length / 2.0 * Math.Cos(Angle), Length / 2.0 * Math.Sin(Angle), 0) + new Vector3d(InteriorLength * Math.Sin(Angle), InteriorLength * Math.Cos(Angle), 0);
                        Rhino.Geometry.Arc ret = new Rhino.Geometry.Arc(Start, InteriorPoint, EndPoint);
                        if (ret.ToNurbsCurve() == null)
                        {
                            return new Rhino.Geometry.LineCurve(Start, Start + new Vector3d(Length * Math.Cos(Angle), Length * Math.Sin(Angle), 0));
                        }
                        else
                        {
                            return ret.ToNurbsCurve();
                        }
                    }

                    public double GetCost() { return 2.0 * (GetCurve(new Point3d(0, 0, 0)).GetLength()); }
                }
            }

        }

        /// <summary>
        /// 現在の進行状況を報告するクラスです。
        /// </summary>
        public class MessageManager
        {
            private System.Windows.Forms.Form form;
            private System.Windows.Forms.TextBox textBox;
            private System.Windows.Forms.ProgressBar progressBar;
            private System.Windows.Forms.Button abortButton;
            public string Text { get { return this.textBox.Text; } set { this.textBox.Text = value; } }
            private System.DateTime startTime;

            public MessageManager()
            {
                form = new System.Windows.Forms.Form();
                form.Width = 500;
                form.Height = 500;

                textBox = new System.Windows.Forms.TextBox();
                textBox.Anchor = ((System.Windows.Forms.AnchorStyles)((((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Bottom)
                  | System.Windows.Forms.AnchorStyles.Left)
                  | System.Windows.Forms.AnchorStyles.Right)));
                textBox.Location = new System.Drawing.Point(0, 0);
                //      textBox.Height = 400;
                textBox.Dock = DockStyle.Fill;
                textBox.Multiline = true;
                //textBox.ReadOnly = true;
                textBox.ScrollBars = System.Windows.Forms.ScrollBars.Both;

                progressBar = new System.Windows.Forms.ProgressBar();
                progressBar.Anchor = ((System.Windows.Forms.AnchorStyles)((((System.Windows.Forms.AnchorStyles.Bottom)
                  | System.Windows.Forms.AnchorStyles.Left)
                  | System.Windows.Forms.AnchorStyles.Right)));
                progressBar.Dock = DockStyle.Bottom;
                progressBar.Height = 50;

                abortButton = new System.Windows.Forms.Button();
                abortButton.Anchor = ((System.Windows.Forms.AnchorStyles)((((System.Windows.Forms.AnchorStyles.Bottom)
                  | System.Windows.Forms.AnchorStyles.Left)
                  | System.Windows.Forms.AnchorStyles.Right)));
                abortButton.Text = "Abort";
                abortButton.Dock = DockStyle.Bottom;
                abortButton.Height = 50;
                abortButton.Click += (a, b) => { abortButton.Enabled = false; throw new Exception(); };

                form.Controls.Add(textBox);
                form.Controls.Add(progressBar);
                form.Controls.Add(abortButton);

                startTime = System.DateTime.Now;
            }

            public void Show()
            {
                form.Show();
            }

            public void WriteLine(string text)
            {
                this.Text += (System.DateTime.Now).ToString() + "\r\n" + text + "\r\n\r\n";

                textBox.Select(textBox.Text.Length, 0);
                textBox.Focus();
                textBox.ScrollToCaret();

                System.Windows.Forms.Application.DoEvents();
            }

            public void SetProgress(int parcentage)
            {
                progressBar.Value = parcentage;
                System.Windows.Forms.Application.DoEvents();
            }
        }
    }
}

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
    public class ClassicPanelling
    {


        /// <summary>
        /// 連続最適化ステップを処理するクラスです。
        /// </summary>
        public class ContinuousOptimizationStep
        {
            /// <summary>
            /// 現在のパネリングにかかわる状態です。
            /// </summary>
            public CurrentStatus Status;
            /// <summary>
            /// 与えられた条件です。
            /// </summary>
            public InitialStatus Argument;
            /// <summary>
            /// 最小二乗法における各要素の重みです。
            /// </summary>
            public List<double> Alpha { get { return _Alpha; } }
            private List<double> _Alpha = new List<double> { 1, 1000, 1, 1, 10 };

            /// <summary>
            /// 閾値に応じたAlpha_Kinkに再設定します。
            /// </summary>
            /// <returns>なし</returns>
            public void UpdateAlphaKink()
            {
                _Alpha[2] = Math.Pow(Argument.DivergenceThreshold / (Argument.KinkAngleThreshold / Math.PI * 180), 2) * Alpha[1];
                //_Alpha[2] = Math.Pow(Status.DivergenceThreshold / Status.KinkAngleThreshold, 2) * Alpha[1];
            }

            /// <summary>
            /// 与えられた条件に従って連続最適化ステップを実行します。
            /// </summary>
            /// <param name="count">最大の実行回数</param>
            public void Evaluate(int count)
            {
                //double BestValue = double.MaxValue;
                double[] BestParamaters = Status.Paramaters;

                for (int k = 0; k < count; k++)
                {
                    Status.StatusMessageManager.SetProgress((k) * 100 / (count));

                    List<List<double>> JacobianTmp = GetJacobianMatrix();
                    //Jacobian's Row and Column is swaped due to bag in function.
                    Rhino.Geometry.Matrix Jacobian = new Rhino.Geometry.Matrix(JacobianTmp[0].Count(), JacobianTmp.Count());
                    for (int i = 0; i < JacobianTmp.Count(); i++)
                    {
                        for (int j = 0; j < JacobianTmp[i].Count(); j++)
                        {
                            Jacobian[j, i] = JacobianTmp[i][j];
                        }
                    }


                    //following is to check the error. fix me.
                    /*
                    List<double> EnergyVectorTmp1 = GetEnergyVector();
                    Status.Paramaters = Status.4;
                    List<double> EnergyVectorTmp2 = GetEnergyVector();
                    double dfa = 0;
                    for(int i = 0;i < EnergyVectorTmp1.Count();i++){
                    dfa += EnergyVectorTmp2[i] - EnergyVectorTmp1[i];
                    }
                    Status.StatusMessageManager.WriteLine("Errorcheck " + dfa);
                    */


                    double EnergyValue = 0.0;
                    List<double> EnergyVectorTmp = GetEnergyVector();
                    Rhino.Geometry.Matrix EnergyVector = new Rhino.Geometry.Matrix(EnergyVectorTmp.Count(), 1);
                    for (int i = 0; i < EnergyVectorTmp.Count(); i++)
                    {
                        EnergyVector[i, 0] = EnergyVectorTmp[i];
                        EnergyValue += EnergyVectorTmp[i] * EnergyVectorTmp[i];
                    }

                    double[] StatusVectorTmp = Status.Paramaters;
                    Rhino.Geometry.Matrix StatusVector = new Rhino.Geometry.Matrix(StatusVectorTmp.GetLength(0), 1);
                    for (int i = 0; i < StatusVectorTmp.GetLength(0); i++)
                    {
                        StatusVector[i, 0] = StatusVectorTmp[i];
                    }


                    Rhino.Geometry.Matrix M1 = Jacobian.Duplicate();
                    M1.Transpose();
                    Rhino.Geometry.Matrix M2 = M1 * Jacobian;
                    M2.Invert(1e-5);//(ZeroTolerance)
                    M1 = M1 * EnergyVector;
                    M1 = M2 * M1;
                    M1.Scale(-1);


                    //New code of length estimation
                    var origStatus = StatusVector.Duplicate();
                    var statusDif = M1.Duplicate();
                    var status = StatusVector.Duplicate();

                    double baseValue = 10.0;
                    double MinEnergySum = GetEnergySum();
                    double MinStatusDifScale = 1.0;

                    for (int j = 1; j < 20; j++)
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
                    for (int j = 0; j < 20; j++)
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
                    SetStatus(status + statusDif);
                    //End of the new code.

                    Status.StatusMessageManager.WriteLine(Status.GetCurrentInformation() + "\r\n" + this.GetEnergySum() + "\r\n" + MinStatusDifScale);
                }

                //Status.StatusMessageManager.WriteLine(Status.GetCurrentInformation());
            }

            private void SetStatus(Rhino.Geometry.Matrix status)
            {
                List<double> ret = new List<double>();
                for (int i = 0; i < status.RowCount; i++)
                {
                    ret.Add(status[i, 0]);
                }
                Status.Paramaters = ret.ToArray();
            }

            /// <summary>
            /// ヤコビアンを求めます。
            /// </summary>
            /// <returns>Rhino.Geometry.Matrixではなく配列の配列の形で返されるヤコビアン。</returns>
            public List<List<double>> GetJacobianMatrix()
            {
                List<List<double>> ret = new List<List<double>>();

                List<double> originalEV = GetEnergyVector();

                double diffpntBase = 1e-6;// 0.0001;
                double diffpnt = 0;// 0.0001;
                //CurveNetwork
                if (Status.FixCurveNetwork == false)
                {
                    for (int i = 0; i < Status.CurveNetworks.Count(); i++)
                    {
                        for (int j = 0; j < Status.CurveNetworks[i].ControlPoints.Count(); j++)
                        {
                            double tmpsum = Math.Abs(Status.CurveNetworks[i].ControlPoints[j].X) + Math.Abs(Status.CurveNetworks[i].ControlPoints[j].Y) + Math.Abs(Status.CurveNetworks[i].ControlPoints[j].Z);
                            diffpnt = tmpsum * diffpntBase;

                            Status.CurveNetworks[i].ControlPoints[j] = new Rhino.Geometry.Point3d(
                              Status.CurveNetworks[i].ControlPoints[j].X + diffpnt, Status.CurveNetworks[i].ControlPoints[j].Y, Status.CurveNetworks[i].ControlPoints[j].Z);
                            ret.Add(GetDifference(GetEnergyVector(), originalEV, 1.0 / diffpnt));
                            Status.CurveNetworks[i].ControlPoints[j] = new Rhino.Geometry.Point3d(
                              Status.CurveNetworks[i].ControlPoints[j].X - diffpnt, Status.CurveNetworks[i].ControlPoints[j].Y, Status.CurveNetworks[i].ControlPoints[j].Z);

                            Status.CurveNetworks[i].ControlPoints[j] = new Rhino.Geometry.Point3d(
                              Status.CurveNetworks[i].ControlPoints[j].X, Status.CurveNetworks[i].ControlPoints[j].Y + diffpnt, Status.CurveNetworks[i].ControlPoints[j].Z);
                            ret.Add(GetDifference(GetEnergyVector(), originalEV, 1.0 / diffpnt));
                            Status.CurveNetworks[i].ControlPoints[j] = new Rhino.Geometry.Point3d(
                              Status.CurveNetworks[i].ControlPoints[j].X, Status.CurveNetworks[i].ControlPoints[j].Y - diffpnt, Status.CurveNetworks[i].ControlPoints[j].Z);

                            Status.CurveNetworks[i].ControlPoints[j] = new Rhino.Geometry.Point3d(
                              Status.CurveNetworks[i].ControlPoints[j].X, Status.CurveNetworks[i].ControlPoints[j].Y, Status.CurveNetworks[i].ControlPoints[j].Z + diffpnt);
                            ret.Add(GetDifference(GetEnergyVector(), originalEV, 1.0 / diffpnt));
                            Status.CurveNetworks[i].ControlPoints[j] = new Rhino.Geometry.Point3d(
                              Status.CurveNetworks[i].ControlPoints[j].X, Status.CurveNetworks[i].ControlPoints[j].Y, Status.CurveNetworks[i].ControlPoints[j].Z - diffpnt);
                        }
                    }
                }

                //Rigid Transformation
                //double diffrtpointBase = 1e-3;// 0.0001;
                double diffrtangle = 1e-3;// 0.00001;
                double diffrtpoint = 0;

                for (int i = 0; i < Status.Panels.Count(); i++)
                {
                    double tmpsum = Math.Abs(Status.RigidTransformations[Status.Panels[i]].TargetPosition.X) + Math.Abs(Status.RigidTransformations[Status.Panels[i]].TargetPosition.Y) + Math.Abs(Status.RigidTransformations[Status.Panels[i]].TargetPosition.Z);
                    //diffrtpoint = tmpsum > 1e3 ? tmpsum * 1e-6 : diffrtpointBase;
                    diffrtpoint = Math.Max(tmpsum * 1e-6, 1e-3);
                    Status.RigidTransformations[Status.Panels[i]].TargetPosition = new Rhino.Geometry.Point3d(
                      Status.RigidTransformations[Status.Panels[i]].TargetPosition.X + diffrtpoint, Status.RigidTransformations[Status.Panels[i]].TargetPosition.Y, Status.RigidTransformations[Status.Panels[i]].TargetPosition.Z);
                    ret.Add(GetDifference(GetEnergyVector(), originalEV, 1.0 / diffrtpoint));
                    Status.RigidTransformations[Status.Panels[i]].TargetPosition = new Rhino.Geometry.Point3d(
                      Status.RigidTransformations[Status.Panels[i]].TargetPosition.X - diffrtpoint, Status.RigidTransformations[Status.Panels[i]].TargetPosition.Y, Status.RigidTransformations[Status.Panels[i]].TargetPosition.Z);

                    Status.RigidTransformations[Status.Panels[i]].TargetPosition = new Rhino.Geometry.Point3d(
                      Status.RigidTransformations[Status.Panels[i]].TargetPosition.X, Status.RigidTransformations[Status.Panels[i]].TargetPosition.Y + diffrtpoint, Status.RigidTransformations[Status.Panels[i]].TargetPosition.Z);
                    ret.Add(GetDifference(GetEnergyVector(), originalEV, 1.0 / diffrtpoint));
                    Status.RigidTransformations[Status.Panels[i]].TargetPosition = new Rhino.Geometry.Point3d(
                      Status.RigidTransformations[Status.Panels[i]].TargetPosition.X, Status.RigidTransformations[Status.Panels[i]].TargetPosition.Y - diffrtpoint, Status.RigidTransformations[Status.Panels[i]].TargetPosition.Z);

                    Status.RigidTransformations[Status.Panels[i]].TargetPosition = new Rhino.Geometry.Point3d(
                      Status.RigidTransformations[Status.Panels[i]].TargetPosition.X, Status.RigidTransformations[Status.Panels[i]].TargetPosition.Y, Status.RigidTransformations[Status.Panels[i]].TargetPosition.Z + diffrtpoint);
                    ret.Add(GetDifference(GetEnergyVector(), originalEV, 1.0 / diffrtpoint));
                    Status.RigidTransformations[Status.Panels[i]].TargetPosition = new Rhino.Geometry.Point3d(
                      Status.RigidTransformations[Status.Panels[i]].TargetPosition.X, Status.RigidTransformations[Status.Panels[i]].TargetPosition.Y, Status.RigidTransformations[Status.Panels[i]].TargetPosition.Z - diffrtpoint);

                    Status.RigidTransformations[Status.Panels[i]].AngleAroundX += diffrtangle;
                    ret.Add(GetDifference(GetEnergyVector(), originalEV, 1.0 / diffrtangle));
                    Status.RigidTransformations[Status.Panels[i]].AngleAroundX -= diffrtangle;

                    Status.RigidTransformations[Status.Panels[i]].AngleAroundY += diffrtangle;
                    ret.Add(GetDifference(GetEnergyVector(), originalEV, 1.0 / diffrtangle));
                    Status.RigidTransformations[Status.Panels[i]].AngleAroundY -= diffrtangle;

                    Status.RigidTransformations[Status.Panels[i]].AngleAroundZ += diffrtangle;
                    ret.Add(GetDifference(GetEnergyVector(), originalEV, 1.0 / diffrtangle));
                    Status.RigidTransformations[Status.Panels[i]].AngleAroundZ -= diffrtangle;
                }

                //Mold Paramaters
                double diffmoldBase = 1e-6;// 0.000000001;
                double diffmold = 0;
                for (int i = 0; i < Status.Molds.Count(); i++)
                {
                    PanellingInfo.Mold tmpm = Status.Molds[i];
                    for (int j = 1; j < tmpm.Paramaters.Count(); j++)
                    {
                        diffmold = Math.Max(Math.Abs(tmpm.Paramaters[j]) * diffmoldBase, 1e-12);
                        tmpm.SetParamater(j, tmpm.Paramaters[j] + diffmold);
                        ret.Add(GetDifference(GetEnergyVector(), originalEV, 1.0 / diffmold));
                        /*
                                  var list = GetDifference(GetEnergyVector(), originalEV, 1.0 / diffmold);
                                  string ttt = "";
                                  foreach(double d in list){
                                    ttt += d + ",";
                                  }
                                  Status.StatusMessageManager.WriteLine(ttt);
                        */
                        tmpm.SetParamater(j, tmpm.Paramaters[j] - diffmold);
                    }
                }
                return ret;
            }

            public double GetEnergySum()
            {
                List<double> ene = GetEnergyVector();
                double ret = 0;
                for (int i = 0; i < ene.Count(); i++)
                {
                    ret += Math.Pow(ene[i], 2);
                }
                return ret;
            }

            /// <summary>
            /// 最小二乗法の対象となる各要素を一覧を求めます。
            /// </summary>
            /// <param name="mul">全体にかけるスカラー値。</param>
            /// <returns>要素の行列。</returns>
            public List<double> GetEnergyVector(double mul)
            {
                List<double> ret = GetEnergyVector();
                for (int i = 0; i < ret.Count(); i++)
                {
                    ret[i] *= mul;
                }
                return ret;
            }

            public List<double> GetDifference(List<double> a, List<double> b, double times)
            {
                for (int i = 0; i < a.Count(); i++)
                {
                    a[i] -= b[i];
                    a[i] *= times;
                }
                return a;
            }

            /// <summary>
            /// 最小二乗法の対象となる各要素を一覧を求めます。
            /// </summary>
            /// <returns>要素の行列。</returns>
            public List<double> GetEnergyVector()
            {
                Rhino.Geometry.Surface sf = Argument.TargetSurface;
                List<double> ret = new List<double>();

                var actualPanels = Status.GetActualPanels();

                for (int i = 0; i < Status.CurveNetworks.Count(); i++)
                {
                    PanellingInfo.CurveNetwork cn = Status.CurveNetworks[i];
                    PanellingInfo.CurveNetwork cnorig = Argument.CurveNetworks[i];

                    Rhino.Geometry.Surface panel1 = actualPanels[Status.CurveNetworkPanelAssignment1[cn]];
                    /*
                    Rhino.Geometry.Surface panel1 = Status.PanelMoldAssignment[Status.CurveNetworkPanelAssignment1[cn]].GetPanel();
                    PanellingInfo.RigidTransformation rg1 = Status.RigidTransformations[Status.CurveNetworkPanelAssignment1[cn]];
                    panel1.Rotate(rg1.AngleAroundX, new Vector3d(1, 0, 0), new Point3d(0, 0, 0));
                    panel1.Rotate(rg1.AngleAroundY, new Vector3d(0, 1, 0), new Point3d(0, 0, 0));
                    panel1.Rotate(rg1.AngleAroundZ, new Vector3d(0, 0, 1), new Point3d(0, 0, 0));
                    panel1.Translate(rg1.TargetPosition.X, rg1.TargetPosition.Y, rg1.TargetPosition.Z);
                    */


                    Rhino.Geometry.Surface panel2 = actualPanels[Status.CurveNetworkPanelAssignment2[cn]];
                    /*
                    Rhino.Geometry.Surface panel2 = Status.PanelMoldAssignment[Status.CurveNetworkPanelAssignment2[cn]].GetPanel();
                    PanellingInfo.RigidTransformation rg2 = Status.RigidTransformations[Status.CurveNetworkPanelAssignment2[cn]];
                    panel2.Rotate(rg2.AngleAroundX, new Vector3d(1, 0, 0), new Point3d(0, 0, 0));
                    panel2.Rotate(rg2.AngleAroundY, new Vector3d(0, 1, 0), new Point3d(0, 0, 0));
                    panel2.Rotate(rg2.AngleAroundZ, new Vector3d(0, 0, 1), new Point3d(0, 0, 0));
                    panel2.Translate(rg2.TargetPosition.X, rg2.TargetPosition.Y, rg2.TargetPosition.Z);
                    */

                    for (int j = 0; j < cn.ControlPoints.Count(); j++)
                    {
                        double u, v;
                        Point3d point = cn.ControlPoints[j];
                        Point3d pointorig = cnorig.ControlPoints[j];

                        Vector3d[] dummyv3d;
                        Point3d cp;

                        if (Status.FixCurveNetwork == false)
                        {
                            //Surface fitting
                            sf.ClosestPoint(point, out u, out v);
                            sf.Evaluate(u, v, 0, out cp, out dummyv3d);
                            ret.Add(cp.DistanceTo(point) * Math.Sqrt(Alpha[0]));
                            //Status.StatusMessageManager.WriteLine("SF " + cp.X + "," + cp.Y + "," + cp.Z + "  " + point.X + "," + point.Y + "," + point.Z + " " + ( cp.DistanceTo(point) * Math.Sqrt(Alpha[0])));
                        }

                        //Divergence
                        double u1, v1, u2, v2;
                        Point3d cp1, cp2;
                        panel1.ClosestPoint(point, out u1, out v1);
                        panel1.Evaluate(u1, v1, 0, out cp1, out dummyv3d);
                        ret.Add(cp1.DistanceTo(point) * Math.Sqrt(Alpha[1]));

                        panel2.ClosestPoint(point, out u2, out v2);
                        panel2.Evaluate(u2, v2, 0, out cp2, out dummyv3d);
                        ret.Add(cp2.DistanceTo(point) * Math.Sqrt(Alpha[1]));

                        //Kink angles
                        /*
                        double tmp = Math.Acos(panel1.NormalAt(u1, v1) * panel2.NormalAt(u2, v2)) * Math.Sqrt(Alpha[2]);
                        if(Double.IsNaN(tmp)){
                        ret.Add(0);
                        }else{
                        ret.Add(tmp);
                        }
                        */
                        ret.Add((panel1.NormalAt(u1, v1) - panel2.NormalAt(u2, v2)).Length * Math.Sqrt(Alpha[2]));

                        if (Status.FixCurveNetwork == false)
                        {
                            //Curve fairness
                            ret.Add(point.DistanceTo(pointorig) * Math.Sqrt(Alpha[3]));
                        }

                        //ToDo: impliment Panel centering.
                    }
                }
                return ret;
            }

        }


        /// <summary>
        /// 離散最適化ステップです。
        /// </summary>
        public class DiscreateOptimizationStep
        {
            /// <summary>
            /// 現在のパネリングにかかわる状態です。
            /// </summary>
            public CurrentStatus Status;
            /// <summary>
            /// 与えられた条件です。
            /// </summary>
            public InitialStatus Argument;

            /// <summary>
            /// Sets initializationにおいて各パネルにCubicMoldを当てはめて計算するためのStatusです。外部で初期化されます。
            /// </summary>
            public CurrentStatus SetInitStatus;

            public void Evaluate()
            {
                List<PanellingInfo.Panel> Sdash = new List<PanellingInfo.Panel>(Status.Panels);
                Dictionary<PanellingInfo.Mold, List<PanellingInfo.Panel>> SigmaDash = new Dictionary<PanellingInfo.Mold, List<PanellingInfo.Panel>>();
                Dictionary<PanellingInfo.Mold, List<PanellingInfo.Panel>> Sigma = new Dictionary<PanellingInfo.Mold, List<PanellingInfo.Panel>>();

                //下のtoleranceは要調節です。
                this.Initialize(out SigmaDash, 3.0);

                /*
                List<PanellingInfo.Panel> SumPanel = new List<PanellingInfo.Panel>();

                foreach(KeyValuePair<PanellingInfo.Mold, List<PanellingInfo.Panel>> kvp in SigmaDash){
                  SumPanel.AddRange(kvp.Value);
                  this.Status.StatusMessageManager.WriteLine("Sigma Count" + kvp.Value.Count());
                }
                foreach(PanellingInfo.Panel panel in this.Status.Panels){
                  if(!SumPanel.Contains(panel)){
                    PanellingInfo.Molds.PolynomialMold pmold = new PanellingInfo.Molds.PolynomialMold();
                    this.Status.Molds.Add(pmold);
                    SigmaDash.Add(pmold, new List<PanellingInfo.Panel>(){panel});
                  }
                }
                */

                var tmpMolds = new List<PanellingInfo.Mold>();
                Status.PanelMoldAssignment = new Dictionary<PanellingInfo.Panel, PanellingInfo.Mold>();
                while (Sdash.Count() > 0)
                {
                    PanellingInfo.Mold mold = this.GetMaximumEfficiency(SigmaDash, Sdash);
                    tmpMolds.Add(mold);
                    Sigma.Add(mold, SigmaDash[mold]);
                    foreach (PanellingInfo.Panel panel in SigmaDash[mold])
                    {
                        Sdash.Remove(panel);
                        if (!Status.PanelMoldAssignment.ContainsKey(panel))
                        {
                            Status.PanelMoldAssignment.Add(panel, mold);
                        }
                    }
                    SigmaDash.Remove(mold);
                }
                Status.Molds = tmpMolds;

            }

            /// <summary>
            /// Efficiencyを取得します。
            /// </summary>
            /// <param name="SigmaDash">SigmaDashです。</param>
            /// <param name="Sdash">Sdashです。</param>
            /// <returns>Efficiencyです。</returns>
            public Dictionary<PanellingInfo.Mold, double> GetEfficiency(Dictionary<PanellingInfo.Mold, List<PanellingInfo.Panel>> SigmaDash, List<PanellingInfo.Panel> Sdash)
            {
                Dictionary<PanellingInfo.Mold, double> Efficiency = new Dictionary<PanellingInfo.Mold, double>();
                foreach (PanellingInfo.Mold mold in SigmaDash.Keys)
                {
                    double cost = mold.GetCost();
                    double count = 0;
                    foreach (PanellingInfo.Panel p in SigmaDash[mold])
                    {
                        if (Sdash.Contains(p))
                        {
                            cost += mold.GetPanelCost();
                            count++;
                        }
                        Efficiency.Add(mold, count / cost);
                    }
                }
                return Efficiency;
            }

            /// <summary>
            /// 実際に必要なのは最大のEfficiencyを持つ型なので、その値を取得します。
            /// </summary>
            /// <param name="SigmaDash">SigmaDashです。</param>
            /// <param name="Sdash">Sdashです。</param>
            /// <returns>最大のEfficiencyを持つ型です。</returns>
            public PanellingInfo.Mold GetMaximumEfficiency(Dictionary<PanellingInfo.Mold, List<PanellingInfo.Panel>> SigmaDash, List<PanellingInfo.Panel> Sdash)
            {
                double efficiency = 0;
                PanellingInfo.Mold ret = this.Status.Molds[0];
                foreach (PanellingInfo.Mold mold in SigmaDash.Keys)
                {
                    double cost = mold.GetCost();
                    double count = 0;
                    foreach (PanellingInfo.Panel p in SigmaDash[mold])
                    {
                        if (Sdash.Contains(p))
                        {
                            cost += mold.GetPanelCost();
                            count++;
                        }
                        if (count / cost > efficiency)
                        {
                            efficiency = count / cost;
                            ret = mold;
                        }
                    }
                }
                return ret;
            }

            /// <summary>
            /// 型の六次元評価を用いて、パネルと型の対応を初期化します。この初期化において一度連続最適化ステップを実行するので実行時間に注意してください。
            /// </summary>
            /// <param name="Sk">Skです。</param>
            /// <param name="tolerance">六次元で評価する際の許容値です。</param>
            /// <returns>なし</returns>
            public void Initialize(out Dictionary<PanellingInfo.Mold, List<PanellingInfo.Panel>> Sk, double tolerance)
            {
                Sk = new Dictionary<PanellingInfo.Mold, List<PanellingInfo.Panel>>();
                SetInitStatus.CurveNetworks = this.Status.CurveNetworks;

                ContinuousOptimizationStep cos = new ContinuousOptimizationStep();
                cos.Argument = this.Argument;
                cos.Status = this.SetInitStatus;
                cos.UpdateAlphaKink();
                this.Status.StatusMessageManager.WriteLine("DOS/COS start");
                cos.Evaluate(50);
                this.Status.StatusMessageManager.WriteLine("DOS/COS finished");

                /*
                foreach(var mold in cos.Status.Molds){
                  if(!Status.Molds.Contains(mold)){
                    Status.Molds.Add(mold);
                  }
                }
                */

                Dictionary<PanellingInfo.Panel, double> MinimumTolerance = new Dictionary<PanellingInfo.Panel, double>();
                foreach (PanellingInfo.Panel p in Status.Panels)
                {
                    MinimumTolerance.Add(p, double.MaxValue);
                }
                foreach (PanellingInfo.Mold m in Status.Molds)
                {
                    foreach (PanellingInfo.Panel p in Status.Panels)
                    {
                        MinimumTolerance[p] = Math.Min(MinimumTolerance[p], cos.Status.PanelMoldAssignment[p].DifferenceFrom(m));
                    }
                }
                double AssumedTolerance = 0;
                foreach (PanellingInfo.Panel p in Status.Panels)
                {
                    AssumedTolerance = Math.Max(AssumedTolerance, MinimumTolerance[p]);
                }
                this.Status.StatusMessageManager.WriteLine("Assumed Tolerance is " + AssumedTolerance);
                foreach (PanellingInfo.Mold m in Status.Molds)
                {
                    Sk.Add(m, new List<PanellingInfo.Panel>());

                    string delme = "Tolerance ";
                    foreach (PanellingInfo.Panel p in Status.Panels)
                    {
                        delme += (cos.Status.PanelMoldAssignment[p].DifferenceFrom(m)) + " ";
                        if (cos.Status.PanelMoldAssignment[p].DifferenceFrom(m) <= AssumedTolerance)
                        {
                            Sk[m].Add(p);
                        }
                    }
                    this.Status.StatusMessageManager.WriteLine(delme);
                }
                /*
                      foreach(PanellingInfo.Panel tmppanel in cos.Status.Panels){
                        this.Status.Molds.Add(cos.Status.PanelMoldAssignment[tmppanel]);
                        Sk[cos.Status.PanelMoldAssignment[tmppanel]] = new List<PanellingInfo.Panel>();
                        Sk[cos.Status.PanelMoldAssignment[tmppanel]].Add(tmppanel);
                      }
                      */
            }
        }

        public class InitializeTransformation
        {
            /// <summary>
            /// 現在のパネリングにかかわる状態です。
            /// </summary>
            public CurrentStatus Status;

            public void Solve()
            {
                SolveAngle(2, 3);
                SolvePosition(1, 10, 1000000);
                SolveAngle(2, 3);
                SolvePosition(2, 15, 1000);
                SolveAngle(2, 10);
            }

            public void SolvePosition(int basecnt, int cnt, double OriginalOrder)
            {
                var minimumDiv = new List<Dictionary<PanellingInfo.Panel, double>>();
                var minimumCnt = new List<Dictionary<PanellingInfo.Panel, int>>();
                var minimumResult = new List<Dictionary<PanellingInfo.Panel, double>>();
                double currentOrder = OriginalOrder;

                var actualpanelnotra = Status.GetActualPanelsNoTransform();

                for (int i = 0; i < 3; i++)
                {
                    minimumDiv.Add(new Dictionary<PanellingInfo.Panel, double>());
                    minimumCnt.Add(new Dictionary<PanellingInfo.Panel, int>());
                    minimumResult.Add(new Dictionary<PanellingInfo.Panel, double>());
                    foreach (PanellingInfo.Panel tmppanel in Status.Panels)
                    {
                        minimumDiv[i].Add(tmppanel, double.MaxValue);
                        minimumCnt[i].Add(tmppanel, -1);
                        minimumResult[i].Add(tmppanel, -1);
                    }
                }
                foreach (PanellingInfo.Panel panel in Status.Panels)
                {
                    minimumResult[0][panel] = Status.RigidTransformations[panel].TargetPosition.X;
                    minimumResult[1][panel] = Status.RigidTransformations[panel].TargetPosition.Y;
                    minimumResult[2][panel] = Status.RigidTransformations[panel].TargetPosition.Z;
                }

                for (int Digit = 0; Digit < cnt; Digit++)
                {
                    currentOrder /= 2 * basecnt;
                    foreach (PanellingInfo.Panel panel in Status.Panels)
                    {
                        minimumDiv[0][panel] = double.MaxValue;
                        minimumDiv[1][panel] = double.MaxValue;
                        minimumDiv[2][panel] = double.MaxValue;
                    }
                    for (int i = -basecnt; i <= basecnt; i++)
                    {
                        foreach (PanellingInfo.Panel panel in Status.Panels)
                        {
                            Status.RigidTransformations[panel].TargetPosition.X = GetNumber(minimumResult[0][panel], currentOrder, i);
                        }
                        var result = Status.GetCurrentDivergenceSquare(actualpanelnotra);
                        foreach (PanellingInfo.Panel panel in Status.Panels)
                        {
                            if (result[panel] < minimumDiv[0][panel])
                            {
                                minimumDiv[0][panel] = result[panel];
                                minimumCnt[0][panel] = i;
                            }
                        }
                    }
                    foreach (PanellingInfo.Panel panel in Status.Panels)
                    {
                        Status.RigidTransformations[panel].TargetPosition.X = minimumResult[0][panel] = GetNumber(minimumResult[0][panel], currentOrder, minimumCnt[0][panel]);
                    }
                    for (int i = -basecnt; i <= basecnt; i++)
                    {
                        foreach (PanellingInfo.Panel panel in Status.Panels)
                        {
                            Status.RigidTransformations[panel].TargetPosition.Y = GetNumber(minimumResult[1][panel], currentOrder, i);
                        }
                        var result = Status.GetCurrentDivergenceSquare(actualpanelnotra);
                        foreach (PanellingInfo.Panel panel in Status.Panels)
                        {
                            if (result[panel] < minimumDiv[1][panel])
                            {
                                minimumDiv[1][panel] = result[panel];
                                minimumCnt[1][panel] = i;
                            }
                        }
                    }
                    foreach (PanellingInfo.Panel panel in Status.Panels)
                    {
                        Status.RigidTransformations[panel].TargetPosition.Y = minimumResult[1][panel] = GetNumber(minimumResult[1][panel], currentOrder, minimumCnt[1][panel]);
                    }
                    for (int i = -basecnt; i <= basecnt; i++)
                    {
                        foreach (PanellingInfo.Panel panel in Status.Panels)
                        {
                            Status.RigidTransformations[panel].TargetPosition.Z = GetNumber(minimumResult[2][panel], currentOrder, i);
                        }
                        var result = Status.GetCurrentDivergenceSquare(actualpanelnotra);
                        foreach (PanellingInfo.Panel panel in Status.Panels)
                        {
                            if (result[panel] < minimumDiv[2][panel])
                            {
                                minimumDiv[2][panel] = result[panel];
                                minimumCnt[2][panel] = i;
                            }
                        }
                    }
                    foreach (PanellingInfo.Panel panel in Status.Panels)
                    {
                        Status.RigidTransformations[panel].TargetPosition.Z = minimumResult[2][panel] = GetNumber(minimumResult[2][panel], currentOrder, minimumCnt[2][panel]);
                    }
                }
            }
            public void SolveAngle(int basecnt, int cnt)
            {
                var minimumDiv = new List<Dictionary<PanellingInfo.Panel, double>>();
                var minimumCnt = new List<Dictionary<PanellingInfo.Panel, int>>();
                var minimumResult = new List<Dictionary<PanellingInfo.Panel, double>>();
                double currentOrder = Math.PI * 2;

                var actualpanelnotra = Status.GetActualPanelsNoTransform();

                for (int i = 0; i < 3; i++)
                {
                    minimumDiv.Add(new Dictionary<PanellingInfo.Panel, double>());
                    minimumCnt.Add(new Dictionary<PanellingInfo.Panel, int>());
                    minimumResult.Add(new Dictionary<PanellingInfo.Panel, double>());
                    foreach (PanellingInfo.Panel tmppanel in Status.Panels)
                    {
                        minimumDiv[i].Add(tmppanel, double.MaxValue);
                        minimumCnt[i].Add(tmppanel, -1);
                        minimumResult[i].Add(tmppanel, -1);
                    }
                }
                foreach (PanellingInfo.Panel panel in Status.Panels)
                {
                    minimumResult[0][panel] = Status.RigidTransformations[panel].AngleAroundX;
                    minimumResult[1][panel] = Status.RigidTransformations[panel].AngleAroundY;
                    minimumResult[2][panel] = Status.RigidTransformations[panel].AngleAroundZ;
                }

                for (int Digit = 0; Digit < cnt; Digit++)
                {
                    currentOrder /= 2 * basecnt;
                    foreach (PanellingInfo.Panel panel in Status.Panels)
                    {
                        minimumDiv[0][panel] = double.MaxValue;
                        minimumDiv[1][panel] = double.MaxValue;
                        minimumDiv[2][panel] = double.MaxValue;
                    }
                    for (int i = -basecnt; i <= basecnt; i++)
                    {
                        for (int j = -basecnt; j <= basecnt; j++)
                        {
                            for (int k = -basecnt; k <= basecnt; k++)
                            {
                                foreach (PanellingInfo.Panel panel in Status.Panels)
                                {
                                    Status.RigidTransformations[panel].AngleAroundX = Math.Tan(GetNumber(minimumResult[0][panel], currentOrder, i));
                                    Status.RigidTransformations[panel].AngleAroundY = Math.Tan(GetNumber(minimumResult[1][panel], currentOrder, j));
                                    Status.RigidTransformations[panel].AngleAroundZ = Math.Tan(GetNumber(minimumResult[2][panel], currentOrder, k));
                                }
                                var result = Status.GetCurrentDivergenceVariance(actualpanelnotra);
                                foreach (PanellingInfo.Panel panel in Status.Panels)
                                {
                                    if (result[panel] < minimumDiv[0][panel])
                                    {
                                        minimumDiv[0][panel] = result[panel];
                                        minimumCnt[0][panel] = i;
                                    }
                                    if (result[panel] < minimumDiv[1][panel])
                                    {
                                        minimumDiv[1][panel] = result[panel];
                                        minimumCnt[1][panel] = j;
                                    }
                                    if (result[panel] < minimumDiv[2][panel])
                                    {
                                        minimumDiv[2][panel] = result[panel];
                                        minimumCnt[2][panel] = k;
                                    }
                                }
                            }
                        }
                    }
                    foreach (PanellingInfo.Panel panel in Status.Panels)
                    {
                        minimumResult[0][panel] = (GetNumber(minimumResult[0][panel], currentOrder, minimumCnt[0][panel]));
                        minimumResult[1][panel] = (GetNumber(minimumResult[1][panel], currentOrder, minimumCnt[1][panel]));
                        minimumResult[2][panel] = (GetNumber(minimumResult[2][panel], currentOrder, minimumCnt[2][panel]));
                    }
                }
                foreach (PanellingInfo.Panel panel in Status.Panels)
                {
                    Status.RigidTransformations[panel].AngleAroundX = Math.Tan(GetNumber(minimumResult[0][panel], currentOrder, minimumCnt[0][panel]));
                    Status.RigidTransformations[panel].AngleAroundY = Math.Tan(GetNumber(minimumResult[1][panel], currentOrder, minimumCnt[1][panel]));
                    Status.RigidTransformations[panel].AngleAroundZ = Math.Tan(GetNumber(minimumResult[2][panel], currentOrder, minimumCnt[2][panel]));
                }
            }

            private double GetNumber(double minresult, double currentorder, int Num)
            {
                //1->0.5,10->0.05
                return minresult + currentorder * Num;
            }
        }


        /// <summary>
        /// Interleaved Iterationを実行するクラスです。Panellingアルゴリズムの本体です。
        /// </summary>
        public class InterleavedIteration
        {
            /// <summary>
            /// 現在のパネリングにかかわる状態です。
            /// </summary>
            public CurrentStatus Status;
            /// <summary>
            /// 与えられた条件です。
            /// </summary>
            public InitialStatus Argument;

            /// <summary>
            /// Interleaved Iterationを実行します。
            /// </summary>
            /// <returns>なし。</returns>
            public void Solve()
            {
                ContinuousOptimizationStep cos = new ContinuousOptimizationStep();
                DiscreateOptimizationStep dos = new DiscreateOptimizationStep();
                InitializeTransformation initrt = new InitializeTransformation();
                InitializeTransformation initrtSetinit = new InitializeTransformation();

                //this.Status = new CurrentStatus();

                this.Status.StatusMessageManager.Show();
                Status.StatusMessageManager.WriteLine("Initialize");

                //以下は参照であることに注意。一一代入しなおす必要はない。はずであったが…。
                cos.Argument = this.Argument;
                cos.Status = this.Status;
                dos.SetInitStatus = new CurrentStatus();
                initrtSetinit.Status = dos.SetInitStatus;
                dos.Status = this.Status;
                dos.Argument = this.Argument;
                initrt.Status = this.Status;

                //逆引き辞書初期化
                Dictionary<PanellingInfo.Panel, List<PanellingInfo.CurveNetwork>> PanelCurveNetworkAssignment = new Dictionary<PanellingInfo.Panel, List<PanellingInfo.CurveNetwork>>();
                foreach (PanellingInfo.Panel panelforrev in Status.Panels)
                {
                    PanelCurveNetworkAssignment.Add(panelforrev, new List<PanellingInfo.CurveNetwork>());
                }

                //初期化
                //Status.CurveNetworks = new List<PanellingInfo.CurveNetwork>(Argument.CurveNetworks);
                Status.CurveNetworks = new List<PanellingInfo.CurveNetwork>();
                foreach (PanellingInfo.CurveNetwork oldcn in Argument.CurveNetworks)
                {
                    var newcn = new PanellingInfo.CurveNetwork(oldcn);
                    Status.CurveNetworks.Add(newcn);

                    PanelCurveNetworkAssignment[Status.CurveNetworkPanelAssignment1[oldcn]].Add(newcn);
                    PanelCurveNetworkAssignment[Status.CurveNetworkPanelAssignment2[oldcn]].Add(newcn);

                    Status.CurveNetworkPanelAssignment1.Add(newcn, Status.CurveNetworkPanelAssignment1[oldcn]);
                    Status.CurveNetworkPanelAssignment2.Add(newcn, Status.CurveNetworkPanelAssignment2[oldcn]);

                    Status.CurveNetworkPanelAssignment1.Remove(oldcn);
                    Status.CurveNetworkPanelAssignment2.Remove(oldcn);
                }
                //      Status.Molds = new List<PanellingInfo.Mold>();
                Status.Molds.Add(new PanellingInfo.Molds.PlaneMold());
                Status.Molds[0].SetParamater(0, Argument.MoldSize);
                Status.PanelMoldAssignment = new Dictionary<PanellingInfo.Panel, PanellingInfo.Mold>();
                Status.RigidTransformations = new Dictionary<PanellingInfo.Panel, PanellingInfo.RigidTransformation>();

                dos.SetInitStatus.CurveNetworks = Status.CurveNetworks;
                dos.SetInitStatus.CurveNetworkPanelAssignment1 = Status.CurveNetworkPanelAssignment1;
                dos.SetInitStatus.CurveNetworkPanelAssignment2 = Status.CurveNetworkPanelAssignment2;
                dos.SetInitStatus.Panels = Status.Panels;
                dos.SetInitStatus.Molds = new List<PanellingInfo.Mold>();
                dos.SetInitStatus.FixCurveNetwork = true;
                dos.SetInitStatus.StatusMessageManager = this.Status.StatusMessageManager;
                dos.SetInitStatus.MaxTransformPoint = this.Status.MaxTransformPoint;
                dos.SetInitStatus.MinTransformPoint = this.Status.MinTransformPoint;

                foreach (PanellingInfo.Panel p in Status.Panels)
                {
                    Status.PanelMoldAssignment.Add(p, Status.Molds[0]);
                    //Status.RigidTransformations.Add(p, new PanellingInfo.RigidTransformation());
                    //}
                    //foreach(PanellingInfo.Panel p in dos.SetInitStatus.Panels){
                    //        PanellingInfo.Mold mold = new PanellingInfo.Molds.PolynomialMold();
                    //PanellingInfo.Mold mold = new PanellingInfo.Molds.PolynomialMold();
                    //PanellingInfo.Mold mold = new PanellingInfo.Molds.PlaneMold();
                    PanellingInfo.Mold mold = new PanellingInfo.Molds.PolynomialMold();
                    mold.SetParamater(0, Argument.MoldSize);
                    dos.SetInitStatus.Molds.Add(mold);
                    dos.SetInitStatus.PanelMoldAssignment.Add(p, mold);


                    //Need to be fixed
                    /*
                    PanellingInfo.RigidTransformation panelRT = new PanellingInfo.RigidTransformation();

                    Point3d StartPoint = PanelCurveNetworkAssignment[p][0].ControlPoints[0];
                    Point3d EndPoint1 = PanelCurveNetworkAssignment[p][0].ControlPoints[PanelCurveNetworkAssignment[p][0].ControlPoints.Count() - 1];
                    Point3d EndPoint2 = PanelCurveNetworkAssignment[p][1].ControlPoints[PanelCurveNetworkAssignment[p][1].ControlPoints.Count() - 1];

                    panelRT.TargetPosition = new Point3d(
                    (EndPoint1 + StartPoint).X / 2.0,
                    (EndPoint1 + StartPoint).Y / 2.0,
                    (EndPoint2 + StartPoint).Z / 2.0);
                    panelRT.AngleAroundZ = Math.Atan((EndPoint1 - StartPoint).X / (EndPoint1 - StartPoint).Y);
                    panelRT.AngleAroundX = Math.Atan((EndPoint1 - StartPoint).Y / (EndPoint1 - StartPoint).Z);
                    */

                    dos.SetInitStatus.RigidTransformations.Add(p, new PanellingInfo.RigidTransformation());
                    Status.RigidTransformations.Add(p, new PanellingInfo.RigidTransformation());

                    //移動範囲枠設定
                    Status.MaxTransformPoint.Add(p, new Point3d(double.MinValue, double.MinValue, double.MinValue));
                    Status.MinTransformPoint.Add(p, new Point3d(double.MaxValue, double.MaxValue, double.MaxValue));
                }
                dos.SetInitStatus.DivergenceThreshold = Argument.DivergenceThreshold;
                dos.SetInitStatus.KinkAngleThreshold = Argument.KinkAngleThreshold;

                foreach (PanellingInfo.CurveNetwork tmpCN in Status.CurveNetworks)
                {
                    foreach (Point3d tmpCNCP in tmpCN.ControlPoints)
                    {
                        Status.MaxTransformPoint[Status.CurveNetworkPanelAssignment1[tmpCN]] =
                          new Point3d(
                          Math.Max(Status.MaxTransformPoint[Status.CurveNetworkPanelAssignment1[tmpCN]].X, tmpCNCP.X),
                          Math.Max(Status.MaxTransformPoint[Status.CurveNetworkPanelAssignment1[tmpCN]].Y, tmpCNCP.Y),
                          Math.Max(Status.MaxTransformPoint[Status.CurveNetworkPanelAssignment1[tmpCN]].Z, tmpCNCP.Z)
                          );
                        Status.MinTransformPoint[Status.CurveNetworkPanelAssignment1[tmpCN]] =
                          new Point3d(
                          Math.Min(Status.MinTransformPoint[Status.CurveNetworkPanelAssignment1[tmpCN]].X, tmpCNCP.X),
                          Math.Min(Status.MinTransformPoint[Status.CurveNetworkPanelAssignment1[tmpCN]].Y, tmpCNCP.Y),
                          Math.Min(Status.MinTransformPoint[Status.CurveNetworkPanelAssignment1[tmpCN]].Z, tmpCNCP.Z)
                          );
                        Status.MaxTransformPoint[Status.CurveNetworkPanelAssignment2[tmpCN]] =
                          new Point3d(
                          Math.Max(Status.MaxTransformPoint[Status.CurveNetworkPanelAssignment2[tmpCN]].X, tmpCNCP.X),
                          Math.Max(Status.MaxTransformPoint[Status.CurveNetworkPanelAssignment2[tmpCN]].Y, tmpCNCP.Y),
                          Math.Max(Status.MaxTransformPoint[Status.CurveNetworkPanelAssignment2[tmpCN]].Z, tmpCNCP.Z)
                          );
                        Status.MinTransformPoint[Status.CurveNetworkPanelAssignment2[tmpCN]] =
                          new Point3d(
                          Math.Min(Status.MinTransformPoint[Status.CurveNetworkPanelAssignment2[tmpCN]].X, tmpCNCP.X),
                          Math.Min(Status.MinTransformPoint[Status.CurveNetworkPanelAssignment2[tmpCN]].Y, tmpCNCP.Y),
                          Math.Min(Status.MinTransformPoint[Status.CurveNetworkPanelAssignment2[tmpCN]].Z, tmpCNCP.Z)
                          );
                    }
                }


                //実行
                int iterationCnt = 10;
                for (int i = iterationCnt; i >= 0; i = i - 1)
                {
                    //Status.StatusMessageManager.SetProgress((iterationCnt - i) * 100 / iterationCnt);
                    Status.StatusMessageManager.WriteLine("Iteration #" + (iterationCnt - i) + " start.");

                    Status.DivergenceThreshold = Argument.DivergenceThreshold + i * 1.0;
                    Status.KinkAngleThreshold = Argument.KinkAngleThreshold + i * Math.PI / 180.0 * 0.5;

                    Status.StatusMessageManager.WriteLine("Discreate Optimization Step start.");
                    initrtSetinit.Solve();
                    Status.StatusMessageManager.WriteLine("Init RT\r\n" + dos.SetInitStatus.GetCurrentInformation());
                    dos.Evaluate();
                    //this.Status = dos.SetInitStatus; throw new Exception();

                    cos.Status = dos.Status;
                    cos.UpdateAlphaKink();
                    Status.StatusMessageManager.WriteLine("Continuous Optimization Step start.");
                    initrt.Solve();
                    //throw new Exception();
                    Status.StatusMessageManager.WriteLine("Init RT\r\n" + Status.GetCurrentInformation());
                    cos.Evaluate(100);
                    this.Status = cos.Status;

                    Status.StatusMessageManager.WriteLine("Re-initialization start.");

                    //re-initialization
                    List<PanellingInfo.Panel> panels;
                    Dictionary<PanellingInfo.Panel, PanellingInfo.Mold> replacedMold = new Dictionary<PanellingInfo.Panel, PanellingInfo.Mold>();

                    List<PanellingInfo.Panel> oldPanels = Status.GetPanelsThatDoesNotSatisfyThreshold(Status.KinkAngleThreshold, Status.DivergenceThreshold);

                    if ((panels = Status.GetPanelsThatDoesNotSatisfyThreshold(Status.KinkAngleThreshold, Status.DivergenceThreshold)).Count() > 0)
                    {
                        List<PanellingInfo.Mold> DoneMold = new List<PanellingInfo.Mold>();
                        foreach (PanellingInfo.Panel panel in panels)
                        {
                            if ((!oldPanels.Contains(panel))) { continue; }
                            if (DoneMold.Contains(Status.PanelMoldAssignment[panel])) { continue; }

                            PanellingInfo.Mold m = new PanellingInfo.Molds.PlaneMold();

                            if (Status.PanelMoldAssignment[panel] is PanellingInfo.Molds.PlaneMold)
                            {
                                m = new PanellingInfo.Molds.CylinderMold();

                            }
                            else if (Status.PanelMoldAssignment[panel] is PanellingInfo.Molds.CylinderMold)
                            {
                                m = new PanellingInfo.Molds.ParaboloidMold();
                            }
                            else if (Status.PanelMoldAssignment[panel] is PanellingInfo.Molds.ParaboloidMold)
                            {
                                m = new PanellingInfo.Molds.TorusMold();
                            }
                            else if (Status.PanelMoldAssignment[panel] is PanellingInfo.Molds.TorusMold)
                            {
                                m = new PanellingInfo.Molds.PolynomialMold();
                            }
                            else
                            {
                                //m = Status.PanelMoldAssignment[panel];
                                oldPanels.Remove(panel);
                                continue;
                            }

                            m.SetParamater(0, Argument.MoldSize);

                            //Status.Molds.Remove(Status.PanelMoldAssignment[panel]);
                            Status.PanelMoldAssignment[panel] = m;
                            Status.Molds.Add(m);

                            DoneMold.Add(m);

                            //Status.RigidTransformations[panel] = new PanellingInfo.RigidTransformation();
                        }
                        initrt.Status = this.Status;
                        cos.Status = this.Status;

                        initrt.Solve();
                        cos.Status.FixCurveNetwork = true;
                        cos.Evaluate(100);
                        cos.Status.FixCurveNetwork = false;
                    }
                    dos.Status = this.Status;


                    /*
                    while((panels = Status.GetPanelsThatDoesNotSatisfyThreshold(Status.KinkAngleThreshold, Status.DivergenceThreshold)).Count() > 0 && oldPanels.Count() > 0){
                    Status.StatusMessageManager.WriteLine(panels.Count() + " panels has failed");
                    foreach(PanellingInfo.Panel panel in panels){
                    if((!oldPanels.Contains(panel))){continue;}

                    PanellingInfo.Mold m = new PanellingInfo.Molds.PlaneMold();
                    if(Status.PanelMoldAssignment[panel] is PanellingInfo.Molds.PlaneMold){
                    Status.StatusMessageManager.WriteLine("Plane mold");
                    m = new PanellingInfo.Molds.CylinderMold();
                    }else if(Status.PanelMoldAssignment[panel] is PanellingInfo.Molds.CylinderMold){
                    Status.StatusMessageManager.WriteLine("Cylindar mold");
                    m = new PanellingInfo.Molds.ParaboloidMold();
                    }else if(Status.PanelMoldAssignment[panel] is PanellingInfo.Molds.ParaboloidMold){
                    Status.StatusMessageManager.WriteLine("Paraboloid mold");
                    m = new PanellingInfo.Molds.TorusMold();
                    }else if(Status.PanelMoldAssignment[panel] is PanellingInfo.Molds.TorusMold){
                    Status.StatusMessageManager.WriteLine("Torus mold");
                    m = new PanellingInfo.Molds.PolynomialMold();
                    }else{
                    Status.StatusMessageManager.WriteLine("Poly mold");
                    //m = Status.PanelMoldAssignment[panel];
                    oldPanels.Remove(panel);
                    continue;
                    }

                    Status.PanelMoldAssignment[panel] = m;
                    Status.Molds.Add(m);
                    Status.RigidTransformations[panel] = new PanellingInfo.RigidTransformation();
                    if(replacedMold.ContainsKey(panel)){
                    Status.Molds.Remove(replacedMold[panel]);
                    replacedMold[panel] = m;
                    }else{
                    replacedMold.Add(panel, m);
                    }
                    }

                    Status.StatusMessageManager.WriteLine("ReInit/COS Start");
                    cos.Status = this.Status;
                    cos.Status.FixCurveNetwork = true;
                    cos.Evaluate(100);
                    cos.Status.FixCurveNetwork = false;
                    dos.Status = cos.Status;
                    this.Status = cos.Status;
                    Status.StatusMessageManager.WriteLine("ReInit/COS End");


                    foreach(PanellingInfo.Panel panel in oldPanels){
                    if(!panels.Contains(panel) && oldPanels.Count() > 1){
                    oldPanels.Remove(panel);
                    }
                    }

                    }

                    */
                    //if(MessageBox.Show("Sucess " + i, "Count", MessageBoxButtons.OKCancel) != DialogResult.OK){throw new Exception();}
                    Status.StatusMessageManager.WriteLine(Status.GetCurrentInformation());
                }
            }
        }



        /// <summary>
        /// 現在の状態を示すクラスです。
        /// </summary>
        public class CurrentStatus
        {
            /// <summary>
            /// カーブネットワークの配列です。各カーブネットワークはパネルの一辺に相当します。
            /// </summary>
            public List<PanellingInfo.CurveNetwork> CurveNetworks;

            public MessageManager StatusMessageManager = new MessageManager();
            /// <summary>
            /// どのパネルにどの型を用いるかを示します。
            /// </summary>
            public Dictionary<PanellingInfo.Panel, PanellingInfo.Mold> PanelMoldAssignment
              = new Dictionary<PanellingInfo.Panel, PanellingInfo.Mold>();
            /// <summary>
            /// 型の一覧です。
            /// </summary>
            public List<PanellingInfo.Mold> Molds = new List<PanellingInfo.Mold>();
            /// <summary>
            /// パネルの一覧です。
            /// </summary>
            public List<PanellingInfo.Panel> Panels = new List<PanellingInfo.Panel>();
            /// <summary>
            /// どのパネルをどのように移動・回転するかを示します。
            /// </summary>
            public Dictionary<PanellingInfo.Panel, PanellingInfo.RigidTransformation> RigidTransformations
              = new Dictionary<PanellingInfo.Panel, PanellingInfo.RigidTransformation>();
            /// <summary>
            /// どのカーブネットワークがどのパネルの一辺を成すかを示します。
            /// </summary>
            public Dictionary<PanellingInfo.CurveNetwork, PanellingInfo.Panel> CurveNetworkPanelAssignment1
              = new Dictionary<PanellingInfo.CurveNetwork, PanellingInfo.Panel>();
            /// <summary>
            /// どのカーブネットワークがどのパネルの一辺を成すかを示します。
            /// </summary>
            public Dictionary<PanellingInfo.CurveNetwork, PanellingInfo.Panel> CurveNetworkPanelAssignment2
              = new Dictionary<PanellingInfo.CurveNetwork, PanellingInfo.Panel>();

            /// <summary>
            /// 現在のKink angleの閾値を示します。
            /// </summary>
            public double KinkAngleThreshold;

            /// <summary>
            /// 現在のDivergenceの閾値を示します。
            /// </summary>
            public double DivergenceThreshold;

            /// <summary>
            /// カーブネットワークを固定します。連続最適化ステップで必要になります。Paramatersでカーブネットワークに相当する値を除外します。
            /// </summary>
            public bool FixCurveNetwork = false;

            public Dictionary<PanellingInfo.Panel, Rhino.Geometry.Point3d> MaxTransformPoint = new Dictionary<PanellingInfo.Panel, Rhino.Geometry.Point3d>();
            public Dictionary<PanellingInfo.Panel, Rhino.Geometry.Point3d> MinTransformPoint = new Dictionary<PanellingInfo.Panel, Rhino.Geometry.Point3d>();



            /// <summary>
            /// 現在の状況を取得します。
            /// </summary>
            /// <returns></returns>
            public string GetCurrentInformation()
            {
                string report = "";
                foreach (PanellingInfo.Panel panel in Panels)
                {
                    report += "Panel information.\r\n";
                    report += "Mold type " + PanelMoldAssignment[panel].GetMoldType + ".\r\n";

                    report += "Paramaters ";
                    foreach (double d in PanelMoldAssignment[panel].Paramaters)
                    {
                        report += d + ",";
                    }
                    report += "\r\n";

                    PanellingInfo.RigidTransformation rg = this.RigidTransformations[panel];
                    report += "Transformation {" + rg.TargetPosition.X + "," + rg.TargetPosition.Y + "," + rg.TargetPosition.Z + "} X:" + rg.AngleAroundX + " Y:" + rg.AngleAroundY + " Z:" + rg.AngleAroundZ + "\r\n";

                    report += "\r\n";
                }
                return report;
            }

            public Dictionary<PanellingInfo.Panel, Rhino.Geometry.Surface> GetActualPanelsNoTransform()
            {
                Dictionary<PanellingInfo.Panel, Rhino.Geometry.Surface> ret = new Dictionary<PanellingInfo.Panel, Rhino.Geometry.Surface>();
                foreach (PanellingInfo.Panel panel in this.Panels)
                {
                    Rhino.Geometry.Surface panel1 = this.PanelMoldAssignment[panel].GetPanel();
                    ret.Add(panel, panel1);
                }
                return ret;
            }


            /// <summary>
            /// 座標変換などを施した実際のパネルを取得します。
            /// </summary>
            /// <returns></returns>
            public Dictionary<PanellingInfo.Panel, Rhino.Geometry.Surface> GetActualPanels()
            {
                return GetActualPanels(GetActualPanelsNoTransform());
            }
            /// <summary>
            /// 座標変換などを施した実際のパネルを取得します。
            /// </summary>
            /// <returns></returns>
            public Dictionary<PanellingInfo.Panel, Rhino.Geometry.Surface> GetActualPanels(Dictionary<PanellingInfo.Panel, Rhino.Geometry.Surface> surfaces)
            {
                Dictionary<PanellingInfo.Panel, Rhino.Geometry.Surface> ret = new Dictionary<PanellingInfo.Panel, Rhino.Geometry.Surface>();
                foreach (PanellingInfo.Panel panel in this.Panels)
                {
                    Rhino.Geometry.Surface panel1 = (Rhino.Geometry.Surface)surfaces[panel].Duplicate();
                    PanellingInfo.RigidTransformation rg1 = this.RigidTransformations[panel];
                    panel1.Rotate(Math.Atan(rg1.AngleAroundX), new Vector3d(1, 0, 0), new Point3d(0, 0, 0));
                    panel1.Rotate(Math.Atan(rg1.AngleAroundY), new Vector3d(0, 1, 0), new Point3d(0, 0, 0));
                    panel1.Rotate(Math.Atan(rg1.AngleAroundZ), new Vector3d(0, 0, 1), new Point3d(0, 0, 0));
                    panel1.Translate(rg1.TargetPosition.X, rg1.TargetPosition.Y, rg1.TargetPosition.Z);

                    //ToDo:型の切り出し。
                    ret.Add(panel, panel1);
                }
                return ret;
            }

            public Dictionary<PanellingInfo.Panel, double> GetCurrentDivergenceSquare(Dictionary<PanellingInfo.Panel, Rhino.Geometry.Surface> surfaces)
            {
                var ret = new Dictionary<PanellingInfo.Panel, double>();

                foreach (PanellingInfo.Panel panel in this.Panels)
                {
                    ret.Add(panel, 0);
                }

                var ActualPanels = this.GetActualPanels(surfaces);
                Rhino.Geometry.Vector3d[] dummyv3d;
                for (int i = 0; i < this.CurveNetworks.Count(); i++)
                {
                    PanellingInfo.CurveNetwork cn = this.CurveNetworks[i];

                    var panel1 = ActualPanels[this.CurveNetworkPanelAssignment1[cn]];
                    var panel2 = ActualPanels[this.CurveNetworkPanelAssignment2[cn]];
                    for (int j = 0; j < cn.ControlPoints.Count(); j++)
                    {
                        Point3d point = cn.ControlPoints[j];

                        //Divergence
                        double u1, v1, u2, v2;
                        Point3d cp1, cp2;
                        panel1.ClosestPoint(point, out u1, out v1);
                        panel1.Evaluate(u1, v1, 0, out cp1, out dummyv3d);
                        ret[this.CurveNetworkPanelAssignment1[cn]] += Math.Pow(cp1.DistanceTo(point), 2.0);
                        //ret[this.CurveNetworkPanelAssignment1[cn]] = Math.Max(cp1.DistanceTo(point), ret[this.CurveNetworkPanelAssignment1[cn]]);
                        panel2.ClosestPoint(point, out u2, out v2);
                        panel2.Evaluate(u2, v2, 0, out cp2, out dummyv3d);
                        ret[this.CurveNetworkPanelAssignment2[cn]] += Math.Pow(cp2.DistanceTo(point), 2.0);
                        //ret[this.CurveNetworkPanelAssignment2[cn]] = Math.Max(cp2.DistanceTo(point), ret[this.CurveNetworkPanelAssignment2[cn]]);
                    }
                }
                return ret;
            }

            public Dictionary<PanellingInfo.Panel, double> GetCurrentDivergenceVariance(Dictionary<PanellingInfo.Panel, Rhino.Geometry.Surface> surfaces)
            {
                var ret = new Dictionary<PanellingInfo.Panel, double>();
                var maxdiv = new Dictionary<PanellingInfo.Panel, double>();
                var mindiv = new Dictionary<PanellingInfo.Panel, double>();

                foreach (PanellingInfo.Panel panel in this.Panels)
                {
                    maxdiv.Add(panel, 0);
                    mindiv.Add(panel, double.MaxValue);
                }

                var ActualPanels = this.GetActualPanels(surfaces);
                Rhino.Geometry.Vector3d[] dummyv3d;
                for (int i = 0; i < this.CurveNetworks.Count(); i++)
                {
                    PanellingInfo.CurveNetwork cn = this.CurveNetworks[i];

                    var panel1 = ActualPanels[this.CurveNetworkPanelAssignment1[cn]];
                    var panel2 = ActualPanels[this.CurveNetworkPanelAssignment2[cn]];
                    for (int j = 0; j < cn.ControlPoints.Count(); j++)
                    {
                        Point3d point = cn.ControlPoints[j];

                        //Divergence
                        double u1, v1, u2, v2;
                        Point3d cp1, cp2;
                        panel1.ClosestPoint(point, out u1, out v1);
                        panel1.Evaluate(u1, v1, 0, out cp1, out dummyv3d);
                        maxdiv[this.CurveNetworkPanelAssignment1[cn]] = Math.Max(cp1.DistanceTo(point), maxdiv[this.CurveNetworkPanelAssignment1[cn]]);
                        mindiv[this.CurveNetworkPanelAssignment1[cn]] = Math.Min(cp1.DistanceTo(point), mindiv[this.CurveNetworkPanelAssignment1[cn]]);

                        panel2.ClosestPoint(point, out u2, out v2);
                        panel2.Evaluate(u2, v2, 0, out cp2, out dummyv3d);
                        maxdiv[this.CurveNetworkPanelAssignment2[cn]] = Math.Max(cp2.DistanceTo(point), maxdiv[this.CurveNetworkPanelAssignment2[cn]]);
                        mindiv[this.CurveNetworkPanelAssignment2[cn]] = Math.Min(cp2.DistanceTo(point), mindiv[this.CurveNetworkPanelAssignment2[cn]]);
                    }
                }
                foreach (PanellingInfo.Panel panel in this.Panels)
                {
                    ret.Add(panel, maxdiv[panel] - mindiv[panel]);
                }
                return ret;
            }

            /*
                /// <summary>
                /// 離散最適化ステップで使うEfficiencyを求めます。
                /// </summary>
                /// <returns>型とEfficiencyの関係。辞書配列です。</returns>
                public Dictionary<PanellingInfo.Mold,double> GetEfficiency(){
                  Dictionary<PanellingInfo.Mold,double> cost = new Dictionary<PanellingInfo.Mold,double>();
                  Dictionary<PanellingInfo.Mold,double> count = new Dictionary<PanellingInfo.Mold,double>();
                  Dictionary<PanellingInfo.Mold,double> ret = new Dictionary<PanellingInfo.Mold,double>();
                  foreach(PanellingInfo.Mold m in this.Molds){
                    cost.Add(m, m.GetCost());
                    count.Add(m, 0);
                  }
                  foreach(KeyValuePair<PanellingInfo.Panel,PanellingInfo.Mold> kvp in this.PanelMoldAssignment){
                    cost[kvp.Value] += kvp.Value.GetPanelCost();
                    count[kvp.Value]++;
                  }
                  foreach(PanellingInfo.Mold m in this.Molds){
                    ret.Add(m, count[m] / cost[m]);
                  }
                  return ret;
                }
            */

            /// <summary>
            /// 現在の閾値を満たさないパネルの一覧を取得します。
            /// </summary>
            /// <returns>閾値を満たさないパネルの一覧。</returns>
            public List<PanellingInfo.Panel> GetPanelsThatDoesNotSatisfyThreshold(double KinkAngleTh, double DivergenceTh)
            {
                List<PanellingInfo.Panel> ret = new List<PanellingInfo.Panel>();
                Dictionary<PanellingInfo.Panel, bool> retOrig = new Dictionary<PanellingInfo.Panel, bool>();
                var ActualPanels = this.GetActualPanels();

                Rhino.Geometry.Vector3d[] dummyv3d;
                for (int i = 0; i < this.CurveNetworks.Count(); i++)
                {
                    PanellingInfo.CurveNetwork cn = this.CurveNetworks[i];

                    var panel1 = ActualPanels[this.CurveNetworkPanelAssignment1[cn]];
                    var panel2 = ActualPanels[this.CurveNetworkPanelAssignment2[cn]];
                    /*
                    Rhino.Geometry.Surface panel1 = this.PanelMoldAssignment[this.CurveNetworkPanelAssignment1[cn]].GetPanel();
                    PanellingInfo.RigidTransformation rg1 = this.RigidTransformations[this.CurveNetworkPanelAssignment1[cn]];
                    panel1.Rotate(rg1.AngleAroundX, new Vector3d(1, 0, 0), new Point3d(0, 0, 0));
                    panel1.Rotate(rg1.AngleAroundY, new Vector3d(0, 1, 0), new Point3d(0, 0, 0));
                    panel1.Rotate(rg1.AngleAroundZ, new Vector3d(0, 0, 1), new Point3d(0, 0, 0));
                    panel1.Translate(rg1.TargetPosition.X, rg1.TargetPosition.Y, rg1.TargetPosition.Z);

                    Rhino.Geometry.Surface panel2 = this.PanelMoldAssignment[this.CurveNetworkPanelAssignment2[cn]].GetPanel();
                    PanellingInfo.RigidTransformation rg2 = this.RigidTransformations[this.CurveNetworkPanelAssignment2[cn]];
                    panel2.Rotate(rg2.AngleAroundX, new Vector3d(1, 0, 0), new Point3d(0, 0, 0));
                    panel2.Rotate(rg2.AngleAroundY, new Vector3d(0, 1, 0), new Point3d(0, 0, 0));
                    panel2.Rotate(rg2.AngleAroundZ, new Vector3d(0, 0, 1), new Point3d(0, 0, 0));
                    panel2.Translate(rg2.TargetPosition.X, rg2.TargetPosition.Y, rg2.TargetPosition.Z);
                */

                    for (int j = 0; j < cn.ControlPoints.Count(); j++)
                    {
                        Point3d point = cn.ControlPoints[j];

                        //Divergence
                        double u1, v1, u2, v2;
                        Point3d cp1, cp2;
                        panel1.ClosestPoint(point, out u1, out v1);
                        panel1.Evaluate(u1, v1, 0, out cp1, out dummyv3d);
                        if (cp1.DistanceTo(point) > DivergenceTh)
                        {
                            //this.StatusMessageManager.WriteLine("Unsatisfied Div Th " + cp1.DistanceTo(point) + "/" + DivergenceTh);
                            if (!retOrig.ContainsKey(this.CurveNetworkPanelAssignment1[cn]))
                            {
                                retOrig.Add(this.CurveNetworkPanelAssignment1[cn], true);
                            }
                        }

                        panel2.ClosestPoint(point, out u2, out v2);
                        panel2.Evaluate(u2, v2, 0, out cp2, out dummyv3d);
                        if (cp2.DistanceTo(point) > DivergenceTh)
                        {
                            //this.StatusMessageManager.WriteLine("Unsatisfied Div Th " + cp2.DistanceTo(point) + "/" + DivergenceTh);
                            if (!retOrig.ContainsKey(this.CurveNetworkPanelAssignment2[cn]))
                            {
                                retOrig.Add(this.CurveNetworkPanelAssignment2[cn], true);
                            }
                        }

                        //Kink angles

                        if (Math.Acos(panel1.NormalAt(u1, v1) * panel2.NormalAt(u2, v2)) > KinkAngleTh)
                        {
                            //this.StatusMessageManager.WriteLine("Unsatisfied Kink Angle " + Math.Acos(panel1.NormalAt(u1, v1) * panel2.NormalAt(u2, v2)) + "/" + KinkAngleTh);
                            if (!retOrig.ContainsKey(this.CurveNetworkPanelAssignment1[cn]))
                            {
                                retOrig.Add(this.CurveNetworkPanelAssignment1[cn], true);
                            }
                            if (!retOrig.ContainsKey(this.CurveNetworkPanelAssignment2[cn]))
                            {
                                retOrig.Add(this.CurveNetworkPanelAssignment2[cn], true);
                            }
                        }
                    }
                }
                foreach (KeyValuePair<PanellingInfo.Panel, bool> kvp in retOrig)
                {
                    ret.Add(kvp.Key);
                }
                return ret;
            }

            public double[] GetMaximumParamatersDrift()
            {
                List<double> ret = new List<double>();

                //Add CN CP
                if (this.FixCurveNetwork == false)
                {
                    for (int i = 0; i < CurveNetworks.Count(); i++)
                    {
                        for (int j = 0; j < CurveNetworks[i].ControlPoints.Count(); j++)
                        {
                            ret.Add(0.5);
                            ret.Add(0.5);
                            ret.Add(0.5);
                        }
                    }
                }

                //Rigid Transformation
                for (int i = 0; i < Panels.Count(); i++)
                {
                    PanellingInfo.RigidTransformation tmpRT =
                      RigidTransformations[Panels[i]];
                    ret.Add(this.FixCurveNetwork ? 0.5 : 0.5);
                    ret.Add(this.FixCurveNetwork ? 0.5 : 0.5);
                    ret.Add(this.FixCurveNetwork ? 0.5 : 0.5);
                    ret.Add(Math.PI / 180.0 * 30.0);
                    ret.Add(Math.PI / 180.0 * 30.0);
                    ret.Add(Math.PI / 180.0 * 30.0);
                }

                //Mold Paramaters
                for (int i = 0; i < Molds.Count(); i++)
                {
                    PanellingInfo.Mold tmpm = Molds[i];

                    for (int j = 1; j < tmpm.Paramaters.Count(); j++)
                    {
                        ret.Add(tmpm.MaximumParamatersDrift[j]);
                    }
                }
                return ret.ToArray();
            }


            /// <summary>
            /// 現在の状態をベクトルとして取り出します。連続最適化ステップで必要になります。
            /// </summary>
            public double[] Paramaters
            {
                get
                {
                    List<double> ret = new List<double>();

                    //Add CN CP
                    if (this.FixCurveNetwork == false)
                    {
                        for (int i = 0; i < CurveNetworks.Count(); i++)
                        {
                            for (int j = 0; j < CurveNetworks[i].ControlPoints.Count(); j++)
                            {
                                ret.Add(CurveNetworks[i].ControlPoints[j].X);
                                ret.Add(CurveNetworks[i].ControlPoints[j].Y);
                                ret.Add(CurveNetworks[i].ControlPoints[j].Z);
                            }
                        }
                    }

                    //Rigid Transformation
                    for (int i = 0; i < Panels.Count(); i++)
                    {
                        PanellingInfo.RigidTransformation tmpRT =
                          RigidTransformations[Panels[i]];
                        ret.Add(tmpRT.TargetPosition.X);
                        ret.Add(tmpRT.TargetPosition.Y);
                        ret.Add(tmpRT.TargetPosition.Z);
                        ret.Add(tmpRT.AngleAroundX);
                        ret.Add(tmpRT.AngleAroundY);
                        ret.Add(tmpRT.AngleAroundZ);
                    }

                    //Mold Paramaters
                    for (int i = 0; i < Molds.Count(); i++)
                    {
                        PanellingInfo.Mold tmpm = Molds[i];

                        for (int j = 1; j < tmpm.Paramaters.Count(); j++)
                        {
                            ret.Add(tmpm.Paramaters[j]);
                        }
                    }
                    return ret.ToArray();
                }
                set
                {
                    //caution: no error check.

                    int count = 0;

                    //FixCurveNetworkはParamatersの読み込み時と書き込み時に同じ値になっていないといけません。
                    //ここで例外が発生した場合、FixCurveNetworkを不適切に書き換えていないか確認してください。
                    if (this.FixCurveNetwork == false)
                    {
                        for (int i = 0; i < CurveNetworks.Count(); i++)
                        {
                            for (int j = 0; j < CurveNetworks[i].ControlPoints.Count(); j++)
                            {
                                CurveNetworks[i].ControlPoints[j] = new Point3d(value[count], value[count + 1], value[count + 2]); count += 3;
                            }
                        }
                    }

                    //Rigid Transformation
                    for (int i = 0; i < Panels.Count(); i++)
                    {
                        PanellingInfo.RigidTransformation tmpRT =
                          RigidTransformations[Panels[i]];
                        /*
                        tmpRT.TargetPosition.X = Math.Max(Math.Min(value[count], this.MaxTransformPoint[Panels[i]].X), this.MinTransformPoint[Panels[i]].X);count++;//fixme
                        tmpRT.TargetPosition.Y = Math.Max(Math.Min(value[count], this.MaxTransformPoint[Panels[i]].Y), this.MinTransformPoint[Panels[i]].Y);count++;//fixme
                        tmpRT.TargetPosition.Z = Math.Max(Math.Min(value[count], this.MaxTransformPoint[Panels[i]].Z), this.MinTransformPoint[Panels[i]].Z);count++;//fixme
                  */
                        tmpRT.TargetPosition.X = value[count]; count++;//fixme
                        tmpRT.TargetPosition.Y = value[count]; count++;//fixme
                        tmpRT.TargetPosition.Z = value[count]; count++;//fixme
                        tmpRT.AngleAroundX = value[count]; count++;
                        tmpRT.AngleAroundY = value[count]; count++;
                        tmpRT.AngleAroundZ = value[count]; count++;
                    }

                    //Mold Paramaters
                    for (int i = 0; i < Molds.Count(); i++)
                    {
                        //this.StatusMessageManager.WriteLine("Mold " + Molds[i].Paramaters.Count() + " type" + Molds[i].GetType());
                        for (int j = 1; j < Molds[i].Paramaters.Count(); j++)
                        {
                            Molds[i].SetParamater(j, value[count]); count++;
                        }
                    }
                }
            }

            /// <summary>
            /// 現在の全てのカーブネットワークに含まれる制御点を取り出します。
            /// </summary>
            /// <returns>制御点の配列。</returns>
            public List<Point3d> GetAllPoints()
            {
                List<Point3d> ret = new List<Point3d>();
                for (int i = 0; i < CurveNetworks.Count(); i++)
                {
                    ret.AddRange(CurveNetworks[i].ControlPoints);
                }
                return ret;
            }

            /// <summary>
            /// 現在の全てのカーブネットワークに含まれる制御点の数を取得します。
            /// </summary>
            /// <returns>全制御点数。</returns>
            public int GetPointCount()
            {
                int cnt = 0;
                for (int i = 0; i < CurveNetworks.Count(); i++)
                {
                    cnt += CurveNetworks[i].ControlPoints.Count();
                }
                return cnt;
            }
        }

        /// <summary>
        /// プログラムに与えられた条件を示します。
        /// </summary>
        public class InitialStatus
        {
            /// <summary>
            /// 与えられたカーブネットワークを示します。
            /// </summary>
            public List<PanellingInfo.CurveNetwork> CurveNetworks;
            /// <summary>
            /// 目的とする形を示します。
            /// </summary>
            public Rhino.Geometry.Surface TargetSurface;

            /// <summary>
            /// 最終的なKink angleの閾値を示します。
            /// </summary>
            public double KinkAngleThreshold;

            /// <summary>
            /// 最終的なDivergenceの閾値を示します。
            /// </summary>
            public double DivergenceThreshold;

            /// <summary>
            /// 型の大きさを示します。
            /// </summary>
            public double MoldSize = 1000;

            /// <summary>
            /// 現在の全てのカーブネットワークに含まれる制御点を取り出します。
            /// </summary>
            /// <returns>制御点の配列。</returns>
            public List<Point3d> GetAllPoints()
            {
                List<Point3d> ret = new List<Point3d>();
                for (int i = 0; i < CurveNetworks.Count(); i++)
                {
                    ret.AddRange(CurveNetworks[i].ControlPoints);
                }
                return ret;
            }
        }

        /// <summary>
        /// パネリングにかかわる情報を示します。
        /// </summary>
        public class PanellingInfo
        {

            /// <summary>
            /// 各パネルの一辺に相当するカーブネットワークを示します。
            /// </summary>
            public class CurveNetwork
            {
                /// <summary>
                /// 制御点です。
                /// </summary>
                public List<Point3d> ControlPoints = new List<Point3d>();

                /// <summary>
                /// カーブネットワークを作成します。別途制御点を追加してください。
                /// </summary>
                public CurveNetwork()
                {
                }
                /// <summary>
                /// カーブネットワークを既存のカーブネットワークから作成します。
                /// </summary>
                /// <param name="orig">元となるカーブネットワーク。</param>
                public CurveNetwork(CurveNetwork orig)
                {
                    foreach (Point3d p3d in orig.ControlPoints)
                    {
                        ControlPoints.Add(p3d);
                    }
                }
            }

            /// <summary>
            /// パネルを示します。現在このクラスでは特に情報を保管していません。識別の為のみに用います。
            /// </summary>
            public class Panel
            {
            }

            /// <summary>
            /// 平行移動およびx軸回り・y軸回りの回転からなる移動を示します。
            /// </summary>
            public class RigidTransformation
            {
                /// <summary>
                /// 平行移動する距離を示します。
                /// </summary>
                public Point3d TargetPosition;
                /// <summary>
                /// x軸回りの回転角(ラジアン)を示します。
                /// </summary>
                public double AngleAroundX { get { return _AngleAroundX; } set { _AngleAroundX = value; } }
                private double _AngleAroundX;
                /// <summary>
                /// y軸回りの回転角(ラジアン)を示します。
                /// </summary>
                public double AngleAroundY { get { return _AngleAroundY; } set { _AngleAroundY = value; } }
                private double _AngleAroundY;
                /// <summary>
                /// y軸回りの回転角(ラジアン)を示します。
                /// </summary>
                public double AngleAroundZ { get { return _AngleAroundZ; } set { _AngleAroundZ = value; } }
                private double _AngleAroundZ;
                /// <summary>
                /// 初期化します。
                /// </summary>
                public RigidTransformation() { _AngleAroundX = 0.01; _AngleAroundY = 0.01; }
                /// <summary>
                /// 既存のクラスをコピーして初期化します。
                /// </summary>
                /// <param name="orig">コピー元</param>
                public RigidTransformation(RigidTransformation orig)
                {
                    this.TargetPosition = orig.TargetPosition;
                    this.AngleAroundX = orig.AngleAroundX;
                    this.AngleAroundY = orig.AngleAroundY;
                }
            }

            /// <summary>
            /// 型を示します。インターフェースです。
            /// </summary>
            public interface Mold
            {
                /// <summary>
                /// 別の型との違いを示します。
                /// </summary>
                /// <param name="target">比較対象。</param>
                /// <returns>違いの大きさ。</returns>
                double DifferenceFrom(Mold target);
                /// <summary>
                /// 型をそれと近似するPolinomialな型として評価します。
                /// </summary>
                /// <returns>相当するMolds.PolynomialMold</returns>
                Molds.PolynomialMold GetPolynomialRepresentation();
                /// <summary>
                /// 型の種類を番号で取得します。
                /// </summary>
                int GetMoldType { get; }
                /// <summary>
                /// 型に関わるパラメータを取得します。
                /// </summary>
                List<double> Paramaters { get; }

                List<double> MaximumParamatersDrift { get; }
                /// <summary>
                /// パラメータに代入します。
                /// </summary>
                /// <param name="i">番号</param>
                /// <param name="d">値</param>
                /// <returns>なし</returns>
                void SetParamater(int i, double d);
                /// <summary>
                /// 実際の型の形を取得します。
                /// </summary>
                /// <returns>実際のサーフェス。</returns>
                Rhino.Geometry.Surface GetPanel();
                /// <summary>
                /// 製造コストを示します。
                /// </summary>
                /// <returns>コスト。</returns>
                double GetCost();
                /// <summary>
                /// 型を使ったパネルの製造コストを示します。
                /// </summary>
                /// <returns>コスト。</returns>
                double GetPanelCost();
            }

            /// <summary>
            /// 具体的な型の一覧です。
            /// </summary>
            public class Molds
            {
                /// <summary>
                /// Polynomial型です。P(x,y)=ax^2+by^2+cx^3+dx^2y+exy^2+fy^3で表現できます。
                /// </summary>
                public class PolynomialMold : Mold
                {
                    private List<double> PolynomialParamaters = new List<double> { 10000, 0, 0, 0, 0, 0, 0 };
                    public List<double> ParamatersMax = new List<double> { 1e6, 1e-2, 1e-2, 1e-4, 1e-4, 1e-4, 1e-4 };
                    public List<double> Paramaters { get { return PolynomialParamaters; } }
                    public List<double> MaximumParamatersDrift { get { return new List<double> { 1e6, 1e-3, 1e-3, 1e-5, 1e-5, 1e-5, 1e-5 }; } }

                    public void SetParamater(int i, double d)
                    {
                        //PolynomialParamaters[i] = d % ParamatersMax[i];
                        PolynomialParamaters[i] = d;
                    }

                    public void SetPolynomialParamaters(double L, double a, double b, double c, double d, double e, double f)
                    {
                        PolynomialParamaters[0] = L;
                        PolynomialParamaters[1] = a;
                        PolynomialParamaters[2] = b;
                        PolynomialParamaters[3] = c;
                        PolynomialParamaters[4] = d;
                        PolynomialParamaters[5] = e;
                        PolynomialParamaters[6] = f;
                    }

                    public double DifferenceFrom(Mold origMold)
                    {
                        List<PolynomialMold> origs = new List<PolynomialMold>();
                        PolynomialMold orig = origMold.GetPolynomialRepresentation();
                        for (int i = 0; i < 8; i++)
                        {
                            origs.Add(new PolynomialMold());
                        }
                        origs[0].SetPolynomialParamaters(orig.Paramaters[0], orig.Paramaters[1], orig.Paramaters[2], orig.Paramaters[3], orig.Paramaters[4], orig.Paramaters[5], orig.Paramaters[6]);
                        origs[1].SetPolynomialParamaters(orig.Paramaters[0], orig.Paramaters[2], orig.Paramaters[1], orig.Paramaters[6], -orig.Paramaters[5], orig.Paramaters[4], -orig.Paramaters[3]);
                        origs[2].SetPolynomialParamaters(orig.Paramaters[0], orig.Paramaters[1], orig.Paramaters[2], -orig.Paramaters[3], -orig.Paramaters[4], -orig.Paramaters[5], -orig.Paramaters[6]);
                        origs[3].SetPolynomialParamaters(orig.Paramaters[0], orig.Paramaters[2], orig.Paramaters[1], -orig.Paramaters[6], orig.Paramaters[5], -orig.Paramaters[4], orig.Paramaters[3]);
                        origs[4].SetPolynomialParamaters(orig.Paramaters[0], -orig.Paramaters[1], -orig.Paramaters[2], -orig.Paramaters[3], orig.Paramaters[4], -orig.Paramaters[5], orig.Paramaters[6]);
                        origs[5].SetPolynomialParamaters(orig.Paramaters[0], -orig.Paramaters[2], -orig.Paramaters[1], -orig.Paramaters[6], -orig.Paramaters[5], -orig.Paramaters[4], -orig.Paramaters[3]);
                        origs[6].SetPolynomialParamaters(orig.Paramaters[0], -orig.Paramaters[1], -orig.Paramaters[2], orig.Paramaters[3], -orig.Paramaters[4], orig.Paramaters[5], -orig.Paramaters[6]);
                        origs[7].SetPolynomialParamaters(orig.Paramaters[0], -orig.Paramaters[2], -orig.Paramaters[1], orig.Paramaters[6], orig.Paramaters[5], orig.Paramaters[4], orig.Paramaters[3]);

                        double ret = double.MaxValue;
                        foreach (PolynomialMold pm in origs)
                        {
                            ret = Math.Min(ret, DifferenceFromSingle(pm));
                        }
                        return ret;
                    }

                    private double DifferenceFromSingle(Mold orig)
                    {
                        double ret = 0;
                        double[] pos6d1 = this.Get6dPosition();
                        double[] pos6d2 = orig.GetPolynomialRepresentation().Get6dPosition();
                        for (int i = 0; i < 6; i++)
                        {
                            ret += Math.Pow(pos6d1[i] - pos6d2[i], 2);
                        }
                        return Math.Sqrt(ret);
                    }

                    public PolynomialMold GetPolynomialRepresentation() { return this; }

                    public double[] Get6dPosition()
                    {
                        double[] ret = {2 * Math.Pow(Paramaters[0], 2) / 3.0 / Math.Sqrt(5) * Paramaters[1],
            2 * Math.Pow(Paramaters[0], 2) / 3.0 / Math.Sqrt(5) * Paramaters[2],
            Math.Pow(Paramaters[0], 3) / Math.Sqrt(15) * (Paramaters[4] + Paramaters[6]),
            Math.Pow(Paramaters[0], 3) / Math.Sqrt(15) * (Paramaters[3] + Paramaters[5]),
            Math.Sqrt(8) * Math.Pow(Paramaters[0], 3) / Math.Sqrt(105) * Paramaters[3],
            Math.Sqrt(8) * Math.Pow(Paramaters[0], 3) / Math.Sqrt(105) * Paramaters[6]};
                        return ret;
                    }

                    public int GetMoldType { get { return 0; } }

                    public double GetCost()
                    {
                        return 10;
                    }
                    public double GetPanelCost() { return 1; }

                    public Rhino.Geometry.Surface GetPanel()
                    {
                        int ucount = 20;
                        int vcount = 20;
                        double L = Paramaters[0];
                        Point3d[] p3ds = new Point3d[ucount * vcount];

                        for (int i = 0; i < ucount; i++)
                        {
                            for (int j = 0; j < vcount; j++)
                            {
                                double x = L * ((double)i / (double)ucount - 0.5);
                                double y = L * ((double)j / (double)vcount - 0.5);

                                p3ds[j * ucount + i] = new Point3d(x, y, GetZAtPoint2d(x, y));
                            }
                        }
                        return Rhino.Geometry.NurbsSurface.CreateFromPoints(p3ds, ucount, vcount, 0, 0);
                    }

                    private double GetZAtPoint2d(double x, double y)
                    {
                        return ((Paramaters[1]) * x * x + (Paramaters[2]) * y * y + (Paramaters[3]) * x * x * x
                          + (Paramaters[4]) * x * x * y + (Paramaters[5]) * x * y * y + (Paramaters[6]) * y * y * y);
                    }
                }

                public class PlaneMold : Mold
                {
                    public List<double> Paramaters { get { return paramaters; } }
                    public List<double> paramaters = new List<double> { 10000 };//L
                    public List<double> MaximumParamatersDrift { get { return new List<double> { 1e6 }; } }

                    public double DifferenceFrom(Mold target)
                    {
                        return this.GetPolynomialRepresentation().DifferenceFrom(target);
                    }

                    public PolynomialMold GetPolynomialRepresentation()
                    {
                        PolynomialMold ret = new PolynomialMold();
                        ret.SetPolynomialParamaters(this.Paramaters[0], 0, 0, 0, 0, 0, 0);
                        return ret;
                    }

                    public int GetMoldType { get { return 1; } }

                    public double GetCost()
                    {
                        return 1;
                    }
                    public double GetPanelCost() { return 1; }

                    public void SetParamater(int i, double d)
                    {
                        paramaters[i] = d;
                    }

                    public Rhino.Geometry.Surface GetPanel()
                    {
                        double L = Paramaters[0];
                        Rhino.Geometry.PolylineCurve plc1 = new Rhino.Geometry.PolylineCurve();
                        Rhino.Geometry.Polyline pl1 = new Rhino.Geometry.Polyline(2);
                        Rhino.Geometry.Polyline pl2 = new Rhino.Geometry.Polyline(2);

                        pl1.Add(-L / 2.0, -L / 2.0, 0);
                        pl1.Add(-L / 2.0, L / 2.0, 0);
                        pl2.Add(L / 2.0, -L / 2.0, 0);
                        pl2.Add(L / 2.0, L / 2.0, 0);

                        return
                          Rhino.Geometry.NurbsSurface.CreateRuledSurface(pl1.ToNurbsCurve(), pl2.ToNurbsCurve());
                    }
                }

                public class ParaboloidMold : Mold
                {
                    //Paraboloid is defined as z=x^2/a^2+y^2/b^2
                    public List<double> Paramaters { get { return paramaters; } }
                    public List<double> paramaters = new List<double> { 10000, 1e4, 1e4 };//L,a,b
                    public List<double> MaximumParamatersDrift { get { return new List<double> { 1e6, double.MaxValue, double.MaxValue }; } }

                    public double DifferenceFrom(Mold target)
                    {
                        return this.GetPolynomialRepresentation().DifferenceFrom(target);
                    }

                    public PolynomialMold GetPolynomialRepresentation()
                    {
                        PolynomialMold ret = new PolynomialMold();
                        ret.SetPolynomialParamaters(this.Paramaters[0], Math.Pow(Paramaters[1], -2.0), Math.Pow(Paramaters[2], -2.0), 0, 0, 0, 0);
                        return ret;
                    }

                    public int GetMoldType { get { return 2; } }

                    public double GetCost()
                    {
                        return 5;
                    }
                    public double GetPanelCost() { return 1; }

                    public void SetParamater(int i, double d)
                    {
                        paramaters[i] = d;
                    }

                    public Rhino.Geometry.Surface GetPanel()
                    {
                        return this.GetPolynomialRepresentation().GetPanel();
                    }
                }

                public class CylinderMold : Mold
                {
                    public List<double> Paramaters { get { return paramaters; } }
                    public List<double> paramaters = new List<double> { 10000, 0, Math.PI * 2.0 };
                    public List<double> MaximumParamatersDrift { get { return new List<double> { 1e6, double.MaxValue, Math.PI / 2.0 }; } }

                    public double DifferenceFrom(Mold target)
                    {
                        return this.GetPolynomialRepresentation().DifferenceFrom(target);
                    }

                    public PolynomialMold GetPolynomialRepresentation()
                    {
                        PolynomialMold ret = new PolynomialMold();
                        ret.SetPolynomialParamaters(this.Paramaters[0], Paramaters[1], 0, 0, 0, 0, 0);
                        return ret;
                    }

                    public int GetMoldType { get { return 3; } }

                    public void SetParamater(int i, double d)
                    {
                        paramaters[i] = d;
                    }

                    public double GetCost()
                    {
                        return 2;
                    }
                    public double GetPanelCost() { return 1; }

                    public Rhino.Geometry.Surface GetPanel()
                    {
                        Rhino.Geometry.Surface ret = this.GetPolynomialRepresentation().GetPanel();
                        ret.Rotate(Paramaters[2], new Vector3d(0, 0, 1), new Point3d(0, 0, 0));
                        return ret;
                        /*
                        Rhino.Geometry.Surface ret = Rhino.Geometry.NurbsSurface.CreateFromCylinder(new Rhino.Geometry.Cylinder(
                          new Rhino.Geometry.Circle(new Rhino.Geometry.Plane(new Point3d(-this.Paramaters[0] / 2.0, 0, -Paramaters[1]), new Vector3d(1, 0, 0)), Math.Abs(Paramaters[1])), this.paramaters[0]));
                        ret.Rotate(Paramaters[2], new Vector3d(0, 0, 1), new Point3d(0, 0, 0));
                        return ret;
                        */
                    }
                }

                public class TorusMold : Mold
                {
                    public List<double> Paramaters { get { return paramaters; } }
                    public List<double> paramaters = new List<double> { 10000, 2e3, 1e3, 0 };
                    public List<double> MaximumParamatersDrift { get { return new List<double> { 1e6, double.MaxValue, double.MaxValue, double.MaxValue }; } }

                    public double DifferenceFrom(Mold target)
                    {
                        return this.GetPolynomialRepresentation().DifferenceFrom(target);
                    }

                    public void SetParamater(int i, double d)
                    {
                        paramaters[i] = d;
                    }

                    public PolynomialMold GetPolynomialRepresentation()
                    {
                        PolynomialMold ret = new PolynomialMold();
                        ret.SetPolynomialParamaters(this.Paramaters[0], -1.0 / 2.0 / Paramaters[1], -1.0 / 2.0 * Math.Cos(Paramaters[3]) / (Paramaters[1] + Paramaters[2]), 0, 0, 0, -1.0 / 2.0 * Paramaters[2] * Math.Sin(Paramaters[3]) / Math.Pow(Paramaters[1] + Paramaters[2], 2));
                        return ret;
                    }

                    public int GetMoldType { get { return 4; } }

                    public double GetCost()
                    {
                        return 4;
                    }
                    public double GetPanelCost() { return 1; }

                    public Rhino.Geometry.Surface GetPanel()
                    {
                        return this.GetPolynomialRepresentation().GetPanel();
                        /*
                        Rhino.Geometry.Surface ret = Rhino.Geometry.NurbsSurface.CreateFromTorus(new Rhino.Geometry.Torus(
                          new Rhino.Geometry.Plane(new Point3d(0, 0, 0), new Vector3d(Math.Cos(Paramaters[3]), 0, Math.Sin(Paramaters[3]))), Math.Abs(Paramaters[1]) + Math.Abs(Paramaters[2]), Math.Abs(Paramaters[2])));
                        ret.Rotate(Paramaters[3], new Vector3d(0, 1, 0), new Point3d(Paramaters[1], 0, 0));
                        ret.Translate(new Vector3d(-Paramaters[1], -Paramaters[2], 0));
                        return ret;
                      */
                    }
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

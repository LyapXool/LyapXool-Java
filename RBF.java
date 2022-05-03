package com.LyapXool;

import java.io.IOException;
import java.time.Duration;
import java.time.Instant;
import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;
class RBF {
    public RBF(int Iode_dimension, int Ic, int Imaxpos, int Imaxneg, double Ialpha, int Ipoints_directional, double Iradius, double Icart_grid_density, double[] Imin_geometric_limits, double[] Imax_geometric_limits, boolean Inormal, boolean Iprinting)
    {
        ode_dimension = Iode_dimension;
        c = Ic;
        maxpos = Imaxpos;
        maxneg = Imaxneg;
        alpha = Ialpha;
        points_directional = Ipoints_directional;
        radius = Iradius;
        cart_grid_density = Icart_grid_density;
        min_geometric_limits = Imin_geometric_limits;
        max_geometric_limits = Imax_geometric_limits;
        normal = Inormal;
        printing = Iprinting;
    }

    public void wbase()
    {
        Instant start = Instant.now();

        rbfbasis=ArrayOperations.ResizeArray(rbfbasis,ode_dimension,ode_dimension);

        for (int i = 1; i <= ode_dimension; ++i)
        {
            for (int j = 1; j <= i; ++j)
            {
                double sqrt = Math.sqrt(1.0 / (2.0 * j * (j + 1.0)));
                if (j == i)
                    rbfbasis[i - 1][j - 1] = (i + 1) * sqrt;
                else
                    rbfbasis[i - 1][j - 1] = sqrt;
            }
        }
        Instant end = Instant.now();
        System.out.println("=====Obtener la base vectorial tomó: "+ Duration.between(start, end).getSeconds() + " segundos=====");


    }
    public void grid()
    {        Instant start = Instant.now();

        int i = 0;
        int j = 0;
        int k = 0;

        int elements;
        int mn = 1;
        elements = Math.abs(maxneg) + Math.abs(maxpos) + 1;
        int[] m = new int[ode_dimension];

        List<Integer> w = new ArrayList<Integer>();
        List<List<Integer>> v = new ArrayList<List<Integer>>(ode_dimension);

        {
            ///ESTO ES v[i].resize(elements)
            for (i = 0; i < ode_dimension; ++i) {
                w.clear();
                for (j = 0; j < elements; ++j) {
                    w.add(maxneg + j);
                }
                v.add(w);
            }
            /// Termina

            for (i = 0; i < ode_dimension; ++i)
            {
                m[i] = v.get(i).size();
                mn *= m[i];
            }


            coord=ArrayOperations.ResizeArray(coord, mn, ode_dimension);
            for (i = 0; i < mn; ++i)
            {
                k = i;
                for (j = ode_dimension - 1; j >= 0; --j)
                {
                    coord[i][j] = v.get(j).get(k % m[j]);
                    k /= m[j];
                }
            }
        }
        Instant end = Instant.now();
        System.out.println("=====Obtener la malla cartesiana en ZxZ tomó: "+ Duration.between(start, end).getSeconds() + " segundos=====");

    }
    public void evaluatinggrid()
    {        Instant start = Instant.now();

        int elements;
        int mn = 1;
        double maxmax = Double.NEGATIVE_INFINITY;
        double minmin = Double.POSITIVE_INFINITY;

        for (int jc = 0; jc < ode_dimension; ++jc)
        {
            if (min_geometric_limits[jc] <= minmin)
            {
                minmin = min_geometric_limits[jc];
            }
            if (max_geometric_limits[jc] >= maxmax)
            {
                maxmax = max_geometric_limits[jc];
            }
        }

        int dimension = 1 + (int)((Math.abs(maxmax) + Math.abs(minmin)) / cart_grid_density);
        double[] evaluatingpoints = new double[dimension];
        evaluatingpoints=Arrays.copyOf(evaluatingpoints, dimension);
        Arrays.fill(evaluatingpoints, 0.0);
        for (int i = 0; i < dimension; ++i)
        {
            evaluatingpoints[i] = minmin + i * cart_grid_density;
        }
        elements = (int)evaluatingpoints.length;

        int[] m = new int[ode_dimension];

        List<Double> w = new ArrayList<Double>();
        List<List<Double>> v = new ArrayList<List<Double>>(ode_dimension);


        for (int i = 0; i != ode_dimension; ++i)
        {
            w.clear();
            int kj = 0;
            for (int j = 0; j != elements; ++j)
            {
                w.add(evaluatingpoints[kj]);
                ++kj;
            }
            v.add(new ArrayList(w));
        }

        for (int i = 0; i < ode_dimension; ++i)
        {
            m[i] = (int)v.get(i).size();
            mn *= m[i];
        }

        cuadricula=ArrayOperations.ResizeArray(cuadricula,(int)mn, (int)ode_dimension);

        for (int i = 0; i < mn; ++i)
        {
            int k = i;
            for (int j = ode_dimension - 1; j >= 0; --j)
            {
                cuadricula[i][j] = v.get(j).get(k % m[j]);
                k /= m[j];
            }
        }
        Instant end = Instant.now();
        System.out.println("=====Obtener la malla de evalucacion cartesiana tomó: "+ Duration.between(start, end).getSeconds() + " segundos=====");
    }
    public double[][] effectivegrid( double[][] gridtobeclean)
    {        Instant start = Instant.now();

        for (int jc = 0; jc < ode_dimension; ++jc)
        {
            if (max_geometric_limits[jc] <= min_geometric_limits[jc])
            {
                System.out.println("ERROR: Maximum should be larger than Minimum");
                System.out.println("Entry: " + jc + " value: " + max_geometric_limits[jc]);

                System.exit(9);
            }
        }
        int dim1 = (int)gridtobeclean.length;//longitud
        int dim2 = (int)gridtobeclean[0].length;//anchura
        List<Integer> counter = new ArrayList<Integer>();
        for (int i = 0; i < dim1; ++i)
        {
            int inside = 0;
            for (int jc = 0; jc < ode_dimension; ++jc)
            {
                if ((gridtobeclean[i][jc] <= max_geometric_limits[jc]
                        && gridtobeclean[i][jc] >= min_geometric_limits[jc]))
                {
                    ++inside;
                }
                else
                {
                    break;
                }
            }
            if (inside == ode_dimension)
            {
                counter.add(i);
            }
        }

        int fin = (int)counter.size();
        double[][] cleanedgrid = new double[fin][dim2];
        int n = 0;
        for (var i = 0; i < fin; ++i)
        {
            for (int j = 0; j < dim2; ++j)
            {
                cleanedgrid[n][j] = gridtobeclean[counter.get(i)][j];
            }
            ++n;
        }
        counter.clear();
        Instant end = Instant.now();
        System.out.println("=====Limpiar la malla en los limities tomó: "+ Duration.between(start, end).getSeconds() + " segundos=====");

        return cleanedgrid;
    }

    public double[][] rbfgrid()
    {        Instant start = Instant.now();

        int i = 0, j = 0;
        int length = (int)Math.pow((Math.abs(maxneg) + Math.abs(maxpos) + 1), ode_dimension);

        gridrbf=ArrayOperations.ResizeArray(gridrbf, length, ode_dimension);
        ArrayOperations.FillMatrixConst(gridrbf, 0.0);

        int dimgpoint1 = (int)coord.length;
        int dimgpoint2 = (int)coord[0].length;
        for (i = 0; i < dimgpoint1; ++i)
        {
            for (j = 0; j < dimgpoint2; ++j)
            {
                for (int All = 0; All < ode_dimension; ++All)
                {
                    gridrbf[i][All] += alpha * coord[i][j] * rbfbasis[j][All];
                }
            }
        }
        {
            for (i = 0; i < dimgpoint1; ++i)
            {
                for (j = 0; j < dimgpoint2; ++j)
                {
                    for (int All = 0; All < ode_dimension; ++All)
                    {
                        gridrbf[i][All] += alpha * coord[i][j] * rbfbasis[j][All];
                    }
                }
                if (dimgpoint2 == 2)
                {
                    for (int All = 0; All < ode_dimension; ++All)
                    {
                        gridrbf[i][All] += 0.5 * alpha * (rbfbasis[0][All] + rbfbasis[1][All]);
                    }
                }
                if (dimgpoint2 == 3)
                {
                    for (int All = 0; All < ode_dimension; ++All)
                    {
                        gridrbf[i][All] += 0.5 * alpha * (rbfbasis[0][All] + rbfbasis[1][All] + rbfbasis[2][All]);
                    }
                }
            }

        }
        Instant end = Instant.now();
        System.out.println("=====Obtener la malla RBF: "+ Duration.between(start, end).getSeconds() + " segundos=====");

        return gridrbf;
    }
    public void alphafunction()
    {        Instant start = Instant.now();

        int N = collocationpoints.length;
        alphavector=Arrays.copyOf(alphavector, N);
        Arrays.fill(alphavector, -1.0);
    }
    public double[][] interpolationmatrixA()
    {
        Instant start = Instant.now();

        int j = 0, k = 0;


        double atzero;

        int dimA = (int)collocationpoints.length;
        int dimAc = (int)collocationpoints[0].length;
        double[][] Amat = new double[dimA][dimA];


        ArrayOperations.FillMatrixConst(Amat, 0.0);


        atzero = Wendland.WndlndFnctnFirst(0.0, c);


        {
            double[] diffsave = new double[dimAc];
            double[] savingcallj = new double[dimAc];
            double[] savingcallk = new double[dimAc];
            double[] resultj = new double[dimAc];
            double[] resultk = new double[dimAc];

            Arrays.fill(diffsave, 0.0);
            Arrays.fill(savingcallj, 0.0);
            Arrays.fill(savingcallk, 0.0);
            Arrays.fill(resultj, 0.0);
            Arrays.fill(resultk, 0.0);

            double twopointsdistance = 0.0;

            double wdlfvalue1 = 0.0;
            double wdlfvalue2 = 0.0;
            double checking = 0.0;


            for (j = 0; j < dimA; ++j)
            {
                for (int all = 0; all < dimAc; ++all)
                    savingcallj[all] = collocationpoints[j][all];


                EqDiff.eqdff(normal, savingcallj, resultj);
                for (k = 0; k < dimA; ++k)
                {
                    for (int all = 0; all < dimAc; ++all)
                        savingcallk[all] = collocationpoints[k][all];

                    diffsave = ArrayOperations.vdiff(savingcallj, savingcallk);
                    if (k == j)
                    {
                        Amat[j][k] = -atzero * ArrayOperations.dot(resultj, resultj);
                    }
                    else
                    {
                        twopointsdistance = Math.sqrt(ArrayOperations.dot(diffsave, diffsave));
                        checking = 1.0 - c * twopointsdistance;
                        if (checking > 0.0)
                        {
                            EqDiff.eqdff(normal, savingcallk, resultk);
                            wdlfvalue1 = Wendland.WndlndFnctnFirst(twopointsdistance, c);
                            wdlfvalue2 = Wendland.WndlndFnctnSecond(twopointsdistance, c);
                            Amat[j][k] = -wdlfvalue2 * ArrayOperations.dot(diffsave, resultj) * ArrayOperations.dot(diffsave, resultk) - wdlfvalue1 * ArrayOperations.dot(resultj, resultk);
                        }
                    }
                }
            }
        }
        Instant end = Instant.now();
        System.out.println("=====Construir la matriz de colocación tomó: "+ Duration.between(start, end).getSeconds() + " segundos=====");

        return Amat;
    }

    public void direcgrid() throws IOException {
        Instant start = Instant.now();

        int lrows = (int)collocationpoints.length;
        int lcols = (int)collocationpoints[0].length;
        stride = points_directional * 2 + 1;
        double[][] coldirectgrid = new double[lrows * stride][lcols];

        int newlenght = (int)(points_directional * 2 * lrows);

        int j, jd;

        double norm;

        directgrid=ArrayOperations.ResizeArray(directgrid, newlenght, lcols);


        double[][] domain = new double[newlenght][lcols];
        double[] savingdomain = new double[lcols];
        double[] evaldfunction = new double[lcols];
        {
            for (int i = 0; i < lrows; ++i)
            {
                j = stride * i;
                jd = (stride - 1) * i;

                for (int All = 0; All < lcols; ++All)
                {
                    coldirectgrid[j][All] = collocationpoints[i][All];
                    savingdomain[All] = collocationpoints[i][All];
                }
                EqDiff.eqdff(normal, savingdomain, evaldfunction);
                norm = Math.sqrt(ArrayOperations.dot(evaldfunction, evaldfunction));
                int kp = 0;
                for (int kd = 0; kd < points_directional; kd += 1)
                {
                    for (int All = 0; All < lcols; ++All)
                    {
                        directgrid[jd + kp][All] = collocationpoints[i][All] + (radius / points_directional) * (kd + 1) * alpha * (evaldfunction[All] / norm);
                        directgrid[jd + kp + 1][All] = collocationpoints[i][All] - (radius / points_directional) * (kd + 1) * alpha * (evaldfunction[All] / norm);
                        coldirectgrid[j + kp + 1][All] = directgrid[jd + kp][All];
                        coldirectgrid[j + kp + 2][All] = directgrid[jd + kp + 1][All];
                    }
                    kp += 2;
                }
            }
        }

        List<Integer> counter = new ArrayList<Integer>();
        List<Integer> counterf = new ArrayList<Integer>();
        {
            int cdrows = (int)coldirectgrid.length;
            int drows = (int)directgrid.length;
            for (int i = 0; i < cdrows; ++i)
            {
                for (int jc = 0; jc < ode_dimension; ++jc)
                {
                    if ((coldirectgrid[i][jc] <= max_geometric_limits[jc]) && (coldirectgrid[i][jc] >= min_geometric_limits[jc]))
                    {
                        counter.add(i);
                    }

                }
            }
            for (int ii = 0; ii < drows; ++ii)
            {
                for (int jc = 0; jc < ode_dimension; ++jc)
                {
                    if ((directgrid[ii][jc] <= max_geometric_limits[jc]) && (directgrid[ii][jc] >= min_geometric_limits[jc]))
                    {
                        counterf.add(ii);
                    }
                }
            }
        }

        int ana = (int)counter.size();
        int flo = (int)counterf.size();

        boolean[] boolcoldirectgrid = new boolean[(int)stride * lrows];
        boolean[] booldirectgrid = new boolean[(int)(lrows * (stride - 1))];

        double[][] cleanbigag = new double[][] { };
        double[][] cleanbigfg = new double[][] { };
        Arrays.fill(boolcoldirectgrid,false);
        Arrays.fill(booldirectgrid,false);


        cleanbigag=ArrayOperations.ResizeArray(cleanbigag, ana, lcols);
        cleanbigfg=ArrayOperations.ResizeArray(cleanbigfg, flo, lcols);

        int n = 0;
        int m = 0;
        int end = counter.size();
        for (int i = 0; i < end; ++i)
        {
            boolcoldirectgrid[counter.get(i)] = true;
            for (int All = 0; All < lcols; ++All)
                cleanbigag[n][All] = coldirectgrid[counter.get(i)][All];
            n++;
        }
        end = counterf.size();
        for (int i = 0; i < end; ++i)
        {
            booldirectgrid[counterf.get(i)] = true;
            for (int All = 0; All < lcols; ++All)
                cleanbigfg[m][All] = directgrid[counterf.get(i)][All];
            m++;
        }

        counter.clear();
        counterf.clear();

        if (printing)
        {
            Generalities.printColumnsToFile("direcgrid", 0, directgrid);
        }
        Instant endt = Instant.now();
        System.out.println("=====Construir la malla direccional tomó: "+ Duration.between(start, endt).getSeconds() + " segundos=====");

    }
    public Matrix makeCholesky(double[][] Amat){
        Instant start = Instant.now();

        Amat=Cholesky.cholesky(Amat);
        Matrix copiedCholesky = new Matrix(Amat);
        Instant end = Instant.now();
        System.out.println("=====Obtener la factorización de Cholesky tomó: "+ Duration.between(start, end).getSeconds() + " segundos=====");

        return copiedCholesky;
    }

    static public int ode_dimension;
    public int c;
    public int maxpos;
    public int maxneg;
    public double alpha;
    public int points_directional;
    public double radius;
    double cart_grid_density;
    public double[] min_geometric_limits;
    public double[] max_geometric_limits;
    public boolean normal;
    public boolean printing;
    public int stride;


    public double[][] coldirectgrid = new double[][] { };
    public double[][] directgrid = new double[][] { };
    public double[][] rbfbasis = new double[][] { };
    public double[][] coord = new double[][] { };
    public double[][] cuadricula = new double[][] { };
    public double[][] gridrbf = new double[][] { };

    public double[] alphavector = new double[] { };
    double[][] collocationpoints = new double[][] { };
    double[][] Amat = new double[][] { };

}

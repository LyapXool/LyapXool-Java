package com.LyapXool;

import java.io.IOException;
import java.time.Duration;
import java.time.Instant;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
public class Lyapunov {
    public Lyapunov(int Iode_dimension, int Ic, int Ipoints_directional, double Icritval, boolean Inormal, boolean Iprinting)
    {
        ode_dimension = Iode_dimension;
        c = Ic;
        points_directional = Ipoints_directional;
        critval = Icritval;
        normal = Inormal;
        printing = Iprinting;
    }

    public double[] lyapequation( double[] alphavector, Matrix copiedCholesky)
    {
        Instant start = Instant.now();
        int maxbet = (int)alphavector.length;
        double[] betaod = new double[maxbet];
        Matrix temp = new Matrix(maxbet, 1);
        temp = copiedCholesky.transpose().solve(copiedCholesky.solve(Matrix.copyVector(alphavector)));
        betaod=Matrix.copyVector(temp,maxbet); //WORKS!!!
        Instant end = Instant.now();
        System.out.println("=====Resolver la ecuación de Liapunov tomó: "+ Duration.between(start, end).getSeconds() + " segundos=====");
        return betaod;
    }
    public double[][] lyapunovfunctions(int currentiteration, boolean type_of_grid, double[] betaodv, double[][] evalcoordinates, double[][] collocationpoints) throws IOException {
        Instant start = Instant.now();
        int i = 0, k = 0;
        int maxite = (int) evalcoordinates.length;
        int maxbet = (int) collocationpoints.length;
        int pointdim = (int) collocationpoints[0].length;
        double[][] lyapFuncOrbDer = new double[maxite][2];


        {

            double[] diffpoints = new double[pointdim];
            double[] diffpointski = new double[pointdim];
            double[] diffpointskineg = new double[pointdim];
            double[] resulti = new double[pointdim];
            double[] resultk = new double[pointdim];
            double[] saving = new double[pointdim];
            double[] savingdomain = new double[pointdim];

            Arrays.fill(diffpoints, 0.0);
            Arrays.fill(diffpointski, 0.0);
            Arrays.fill(diffpointskineg, 0.0);
            Arrays.fill(resulti, 0.0);
            Arrays.fill(resultk, 0.0);
            Arrays.fill(saving, 0.0);
            Arrays.fill(savingdomain, 0.0);

            double proctk = 0.0;
            double producting = 0.0;
            double twopointsdistance = 0.0;
            double wdlfvalue1 = 0.0;
            double wdlfvalue2 = 0.0;
            double checking = 0.0;
            for (i = 0; i < maxite; ++i) {
                for (int all = 0; all < pointdim; ++all)
                    savingdomain[all] = evalcoordinates[i][all];

                EqDiff.eqdff(normal, savingdomain, resulti);
                for (k = 0; k < maxbet; ++k) {
                    for (int all = 0; all < pointdim; ++all)
                        saving[all] = collocationpoints[k][all];

                    diffpoints = ArrayOperations.vdiff(savingdomain, saving);
                    twopointsdistance = Math.sqrt(ArrayOperations.dot(diffpoints, diffpoints));
                    checking = 1.0 - c * twopointsdistance;
                    if (checking > 0.0) {
                        EqDiff.eqdff(normal, saving, resultk);
                        wdlfvalue1 = Wendland.WndlndFnctnFirst(twopointsdistance, c);
                        wdlfvalue2 = Wendland.WndlndFnctnSecond(twopointsdistance, c);
                        diffpointski = ArrayOperations.vdiff(saving, savingdomain);
                        proctk = ArrayOperations.dot(diffpointski, resultk);
                        producting = betaodv[k] * proctk;
                        lyapFuncOrbDer[i][0] += producting * wdlfvalue1;
                        lyapFuncOrbDer[i][1] += -wdlfvalue2 * producting * ArrayOperations.dot(diffpointski, resulti) -
                                betaodv[k] * wdlfvalue1 * ArrayOperations.dot(resulti, resultk);
                    }
                }
            }
        }
        double[] lyapFunc = new double[maxite];
        double[] orbDer = new double[maxite];
        for (int j = 0; j < maxite; ++j)
        {
            lyapFunc[j] = lyapFuncOrbDer[j][0];
            orbDer[j]  =lyapFuncOrbDer[j][1];
        }


        if(type_of_grid){
            Generalities.printVectorToFile("lyapFuncd",currentiteration,lyapFunc);
            Generalities.printVectorToFile("OrbDerd",currentiteration,orbDer);

        }else{
            Generalities.printVectorToFile("lyapFuncc",currentiteration,lyapFunc);
            Generalities.printVectorToFile("orbDerc",currentiteration,orbDer);
        }
        Instant end = Instant.now();
        System.out.println("=====Obtener la función de Liapunov tomó: "+ Duration.between(start, end).getSeconds() + " segundos=====");

        return  lyapFuncOrbDer;
    }
    public void chainrecurrentset(int currentiteration, boolean type_of_grid, boolean with_orbder, double[] lyapFunc, double[] orbDer, double[][] evalcoordinates, double[] vector) throws IOException {
        Instant start = Instant.now();
        List<Integer> counterzero = new ArrayList<Integer>();
        int maxlength = (int)evalcoordinates.length;
        int maxwidth = (int)evalcoordinates[0].length;

        for (int j = 0; j < maxlength; ++j)
        {
            if (with_orbder)
            {
                if (vector[j] > critval)
                {
                    counterzero.add(j);
                }
            }
            else
            {
                if (-vector[j] > -critval)
                {
                    counterzero.add(j);
                }
            }
        }
        int faillength = (int)counterzero.size();

        double[] crslyapun = new double[faillength];
        double[] crsorbder = new double[faillength];
        double[][] failinggrid = new double[faillength][maxwidth];
        double[] failinglyapunov = new double[faillength];
        double[] failingorbder = new double[faillength];
        int end = counterzero.size();
        int m = 0;
        {
            for (int ii = 0; ii < end; ++ii)
            {
                crslyapun[m] = lyapFunc[counterzero.get(ii)];
                crsorbder[m] = orbDer[counterzero.get(ii)];
                failinglyapunov[m] = lyapFunc[counterzero.get(ii)];
                failingorbder[m] = orbDer[counterzero.get(ii)];

                for (int All = 0; All < maxwidth; ++All)
                    failinggrid[m][All] = evalcoordinates[counterzero.get(ii)][All];

                m++;
            }
        }
        counterzero.size();
        if (printing)
        {
            if (type_of_grid)
            {
                Generalities.printColumnsToFile("fdirecgrid", currentiteration, failinggrid);
                Generalities.printVectorToFile("flfdirecgrid", currentiteration,  failinglyapunov);
                Generalities.printVectorToFile("flfpdirecgrid", currentiteration,  failingorbder);
            }
            else
            {
                Generalities.printColumnsToFile("fcartesian", currentiteration, failinggrid);
                Generalities.printVectorToFile("flfcartesian", currentiteration, failinglyapunov);
                Generalities.printVectorToFile("flfpcartesian", currentiteration, failingorbder);
            }
        }
        Instant endt = Instant.now();
        System.out.println("=====Obtener las recurrencias tomó: "+ Duration.between(start, endt).getSeconds() + " segundos=====");

    }
    public double[] firstderivative(int currentiteration, boolean normal, boolean type_of_grid,  double[][] evalcoordinates, RBF rbf) throws IOException {
        Instant start = Instant.now();
        int i = 0, j = 0, k = 0;
        int evaldim = (int)evalcoordinates.length;
        double[][] fdvector = new double[evaldim][ode_dimension];
        double[] normed = new double[evaldim];

        ArrayOperations.FillMatrixConst(fdvector, 0.0);
        double checking = 0.0;
        double twopointsdistance = 0.0;
        double wdlfvalue1 = 0.0;
        double wdlfvalue2 = 0.0;
        {
            double[] saving = new double[ode_dimension];
            double[] savingdomain = new double[ode_dimension];
            double[] diffpoints = new double[ode_dimension];
            double[] resultk = new double[ode_dimension];

            int maxite = (int)betaod.length;

            for (j = 0; j < evaldim; ++j)
            {
                for (int All = 0; All < ode_dimension; ++All)
                    saving[All] = evalcoordinates[j][All];

                for (i = 0; i < ode_dimension; ++i)
                {
                    for (k = 0; k < maxite; ++k)
                    {
                        for (int All = 0; All < ode_dimension; ++All)
                            savingdomain[All] = rbf.collocationpoints[k][All];

                        EqDiff.eqdff(normal, savingdomain, resultk);
                        diffpoints = ArrayOperations.vdiff(saving, savingdomain);
                        twopointsdistance = Math.sqrt(ArrayOperations.dot(diffpoints, diffpoints));
                        checking = 1.0 - c * twopointsdistance;
                        if (checking > 0.0)
                        {
                            wdlfvalue1 = Wendland.WndlndFnctnFirst(twopointsdistance, c);
                            wdlfvalue2 = Wendland.WndlndFnctnSecond(twopointsdistance, c);
                            fdvector[j][i] += betaod[k] * (-resultk[i] * wdlfvalue1
                                    - diffpoints[i]
                                    * ArrayOperations.dot(diffpoints, resultk)
                                    * wdlfvalue2);
                        }
                    }
                }
            }
        }
        double numbernormsquare = 0.0;
        normed=Arrays.copyOf(normed, evaldim);
        double[] temp1 = new double[ode_dimension];
        for (int p = 0; p < evaldim; ++p)
        {
            for (int All = 0; All < ode_dimension; ++All)
            {
                temp1[All] = fdvector[p][All];
            }
            numbernormsquare = ArrayOperations.dot(temp1, temp1);
            normed[p] = Math.sqrt(numbernormsquare);
        }
        if (printing)
        {
            if (type_of_grid)
            {
                Generalities.printColumnsToFile("lyapprimexdir", 0, fdvector);
                Generalities.printVectorToFile("normeddire", 0, normed);
            }
            else
            {
                Generalities.printColumnsToFile("lyapprimexcar", 0, fdvector);
                Generalities.printVectorToFile("normedcar", 0,  normed);
            }
        }
        Instant end = Instant.now();
        System.out.println("=====Obtener el gradiente y su norma: "+ Duration.between(start, end).getSeconds() + " segundos=====");

        return normed;
    }
    public double[] getnewalpha(int currentiteration, double[] avector, double[] orbDer) throws IOException {
        Instant start = Instant.now();
        double summing=0.0;
        double normalizationfactor=0.0;
        int N=avector.length;
        avector=Arrays.copyOf(avector, N);
        Arrays.fill(avector,0.0);
        for(int iii=0; iii<N; ++iii)
        {
            summing=0.0;
            for(int j=0; j<2*points_directional;++j)
            {
                summing+=orbDer[(2*points_directional)*(iii)+j];
            }
            if(summing>0.0)
            {
                summing=0.0;
            }
            avector[iii]=summing/((double)(2*points_directional));
            normalizationfactor+=avector[iii];
        }
        avector=ArrayOperations.vectorMultiply((double)Math.abs(N/normalizationfactor),avector);
        if(printing)
        {
            Generalities.printVectorToFile("alphavector",  currentiteration, avector);
        }
        Instant end = Instant.now();
        System.out.println("=====Obtener el nuevo vector alfa tomó: "+ Duration.between(start, end).getSeconds() + " segundos=====");

        return avector;
    }

    public double[] betaod = new double[] { };
    public int ode_dimension;
    public int c;
    public int points_directional;
    public double critval;
    public boolean normal;
    public boolean printing;

}

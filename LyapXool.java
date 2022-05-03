package com.LyapXool;


import java.io.IOException;
import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.util.Calendar;
import java.util.Locale;

public class LyapXool {
    public static void main(String[] args) throws IOException {
        LocalDate localDate= LocalDate.now();
        Locale spanishLocale=new Locale("es", "ES");
        Calendar now = Calendar.getInstance();
        String dateInSpanish=localDate.format(DateTimeFormatter.ofPattern("EEEE dd MMMM de yyyy",spanishLocale));
        System.out.println("El cálculo de la función de Liapunov comenzó el día " + dateInSpanish + " a las " + now.get(Calendar.HOUR_OF_DAY)+":"+ now.get(Calendar.MINUTE)+":"+now.get(Calendar.SECOND)+"."+now.get(Calendar.MILLISECOND)+" horas.");
        Generalities.PrintGnlInfo();

        RBF rbf = new RBF(Instructions.ode_dimension, Instructions.c, Instructions.maxmax, Instructions.minmin, Instructions.alpha, Instructions.points_directional, Instructions.radius, Instructions.cart_grid_density, Instructions.min_geometric_limits, Instructions.max_geometric_limits, Instructions.normal, Instructions.printing);
        Lyapunov lpv = new Lyapunov(Instructions.ode_dimension, Instructions.c, Instructions.points_directional, Instructions.critval, Instructions.normal, Instructions.printing);

        double[][] cartesianevalgrid = new double[][] { };
        double[][] lyapFuncOrbDer  = new double[][] { };
        rbf.wbase();
        rbf.grid();
        rbf.evaluatinggrid();
        cartesianevalgrid=rbf.effectivegrid( rbf.cuadricula);
        Generalities.printColumnsToFile("cg", cartesianevalgrid);
        rbf.rbfgrid();
        rbf.collocationpoints=rbf.effectivegrid(rbf.gridrbf);
        Generalities.printColumnsToFile("collocation", rbf.collocationpoints);
        rbf.alphafunction();
        rbf.direcgrid();
        rbf.Amat=rbf.interpolationmatrixA();
        Matrix copiedCholesky=rbf.makeCholesky(rbf.Amat);


        for(int i=1; i<=Instructions.totaliterations; ++i){
            System.out.println("\t\t\t===== Iteración "+ i +"=====");
            lpv.betaod=lpv.lyapequation(rbf.alphavector, copiedCholesky);
            lyapFuncOrbDer=lpv.lyapunovfunctions(i, true, lpv.betaod, rbf.directgrid, rbf.collocationpoints);

            int maxite =  rbf.directgrid.length;
            double[] lyapFunc = new double[maxite];
            double[] orbDer = new double[maxite];

            for (int j = 0; j < maxite; ++j)
            {
                lyapFunc[j] = lyapFuncOrbDer[j][0];
                orbDer[j]  =lyapFuncOrbDer[j][1];
            }

            lpv.chainrecurrentset(i, true, true, lyapFunc, orbDer, rbf.directgrid, orbDer);

            rbf.alphavector=lpv.getnewalpha(i, rbf.alphavector, orbDer);

            lyapFuncOrbDer=lpv.lyapunovfunctions(i, false, lpv.betaod, cartesianevalgrid, rbf.collocationpoints);
            lpv.firstderivative(i, Instructions.normal, false,  cartesianevalgrid, rbf);

            maxite = cartesianevalgrid.length;
            double[] lyapFuncc = new double[maxite];
            double[] orbDerc = new double[maxite];
            for (int j = 0; j < maxite; ++j)
            {
                lyapFuncc[j] = lyapFuncOrbDer[j][0];
                orbDerc[j]  =lyapFuncOrbDer[j][1];
            }

            lpv.chainrecurrentset(i, true, true, lyapFuncc, orbDerc, rbf.directgrid, orbDer);

        }
        System.out.println("El cálculo de la función de Liapunov terminó el día " + dateInSpanish + " a las " + now.get(Calendar.HOUR_OF_DAY)+":"+ now.get(Calendar.MINUTE)+":"+now.get(Calendar.SECOND)+"."+now.get(Calendar.MILLISECOND)+" horas.");

    }
}

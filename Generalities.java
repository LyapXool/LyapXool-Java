package com.LyapXool;

import java.io.FileWriter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.util.Calendar;
import java.util.Date;
import java.util.Locale;

public class Generalities {
    static void printMAtrix(String nameVariable, int interationcount, double[][] matrixToPrint){
        int rows=matrixToPrint.length;
        int cols=matrixToPrint[0].length;
        System.out.print(nameVariable + "ite" + interationcount + " = [ ");
        for (int i = 0; i < rows; i++)
        {
            for (int j=0; j<cols;++j)
            {
                System.out.print(matrixToPrint[i][j]+" ");
                System.out.println();
            }
        }
        System.out.println("];");
    }
    static public void printVector(String nameVariable, int interationcount, double[] vectorToPrint){
        int end=vectorToPrint.length;
        System.out.print(nameVariable + "ite" + interationcount + " = [ ");
        for(int i=0; i<end;++i)
        {
            System.out.print(vectorToPrint[i]+" ");
        }
        System.out.println("];");
    }
    static public void printColumnsToFile(String nameVariable, int interationcount, double[][] matrixToPrint) throws IOException {
        int minRows = matrixToPrint.length;
        int minCols = matrixToPrint.length;

        String fileName1 = "s" + nameVariable + "ite" + interationcount + "x.m";

        FileWriter sw1 = new FileWriter(fileName1);

        sw1.write(nameVariable + "x = [ ");
        for (int i = 0; i < minRows; ++i)
        {
            sw1.write(matrixToPrint[i][0] + " ");
        }
        sw1.write("];");
        sw1.close();

        String fileName2 = "s" + nameVariable + "ite" + interationcount + "y.m";
        FileWriter sw2 = new FileWriter(fileName2);

        sw2.write(nameVariable + "y = [ ");
        for (int i = 0; i < minRows; ++i)
        {
            sw2.write(matrixToPrint[i][1] + " ");
        }
        sw2.write("];");
        sw2.close();
    }
    static public void printColumnsToFile(String nameVariable, double[][] matrixToPrint) throws IOException {
        int minRows = matrixToPrint.length;
        int minCols = matrixToPrint.length;

        String fileName1 = "s" + nameVariable + "ite"  + "x.m";

        FileWriter sw1 = new FileWriter(fileName1);

        sw1.write(nameVariable + "x = [ ");
        for (int i = 0; i < minRows; ++i)
        {
            sw1.write(matrixToPrint[i][0] + " ");
        }
        sw1.write("];");
        sw1.close();

        String fileName2 = "s" + nameVariable + "ite" + "y.m";
        FileWriter sw2 = new FileWriter(fileName2);

        sw2.write(nameVariable + "y = [ ");
        for (int i = 0; i < minRows; ++i)
        {
            sw2.write(matrixToPrint[i][1] + " ");
        }
        sw2.write("];");
        sw2.close();
    }
    static public void printVectorToFile(String nameVariable, int interationcount, double[] vectorToPrint) throws IOException {
        int minRows = vectorToPrint.length;
        String fileName1 = "s" + nameVariable + "ite" + interationcount+ ".m";
        FileWriter sw1 = new FileWriter(fileName1);

        sw1.write(nameVariable + " = [ ");
        for (int i = 0; i < minRows; ++i)
        {
            sw1.write(vectorToPrint[i] + " ");
        }
        sw1.write("];");
        sw1.close();
    }



    static  public void PrintGnlInfo() throws IOException  {
        String fileName1 = "data.lpx";

        FileWriter sw1 = new FileWriter(fileName1);

        sw1.write("The parameters used to produced these computations are the following: \n");
        sw1.write("\n");
        sw1.write("ode_dimension: " + Instructions.ode_dimension + "\n");
        sw1.write("c: " + Instructions.c + "\n");
        sw1.write("maxmax: " + Instructions.maxmax + "\n");
        sw1.write("minmin: " + Instructions.minmin + "\n");
        sw1.write("alpha: " + Instructions.alpha + "\n");
        sw1.write("normal: " + Instructions.normal + "\n");
        sw1.write("cart_grid_density: " + Instructions.cart_grid_density + "\n");
        sw1.write("radius: " + Instructions.radius + "\n");
        sw1.write("critval: " + Instructions.critval + "\n");
        sw1.write("min_geometric_limits: { ");
        for (int i = 0; i < Instructions.ode_dimension; ++i)
        {
            sw1.write(Instructions.min_geometric_limits[i] + " ");
        }
        sw1.write("}" + "\n");
        sw1.write("max_geometric_limits: { ");
        for (int i = 0; i < Instructions.ode_dimension; ++i)
        {
            sw1.write(Instructions.max_geometric_limits[i] + " ");
        }
        sw1.write("}" + "\n");
        sw1.write("points_directional: " + Instructions.points_directional + "\n");
        sw1.write("printing: " + Instructions.printing + "\n");
        sw1.write("\n");
        sw1.write("\n");
        sw1.write("\n");
        sw1.close();
        PrintDisc();
    }

    static public void PrintDisc() throws IOException {

        LocalDate localDate= LocalDate.now();
        Locale spanishLocale=new Locale("es", "ES");
        Calendar now = Calendar.getInstance();
        String dateInSpanish=localDate.format(DateTimeFormatter.ofPattern("EEEE dd MMMM de yyyy",spanishLocale));

        String fileName1 = "gnldisc.lpx";

        FileWriter sw1 = new FileWriter(fileName1);

        sw1.write("El cálculo de la función de Liapunov comenzó el día " + dateInSpanish + " a las " + now.get(Calendar.HOUR_OF_DAY)+":"+ now.get(Calendar.MINUTE)+":"+now.get(Calendar.SECOND)+"."+now.get(Calendar.MILLISECOND)+" horas.\n\n" +
                "LyapXool# has been written by Carlos Argáez from the Marine and Freshwater Research Institute, Iceland.\n" +
                "The code is primarily based on C. Argáez's own code: LyapXool, Argáez et al. DOI: 10.1016/j.softx.2020.100616\n" +
                "However, it does use a great deal of public sources: \n" +
                "\t 1)  textbooks for a first course in computer science for the next generation of scientists and engineers.\n" +
                "\t https://introcs.cs.princeton.edu/java/home/.\n"+
                "\t In particular, the Java Class on Matrix https://introcs.cs.princeton.edu/java/95linear/Matrix.java.html\n" +
                "\t https://jamesmccaffrey.wordpress.com/2015/03/06/inverting-a-matrix-using-c/ \n" +
                "\t 2) textbooks for a first course in computer science for the next generation of scientists and engineers.\n" +
                "\t https://introcs.cs.princeton.edu/java/home/.\n"+
                "\t In particular, the Java Class on Cholesky https://introcs.cs.princeton.edu/java/95linear/Cholesky.java.html\n");
        sw1.write("\n");
        sw1.write("\n");
        sw1.write("\n");
        sw1.write("Everybody can download, use and modify this code.");
        sw1.close();
    }

}

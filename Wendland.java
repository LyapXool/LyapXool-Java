package com.LyapXool;

public class Wendland {
    static public double WndlndFnctn(double r, int c)
    {
        double psi = 0.0;
        double crs = c * r;
        double maximo = Math.max(1.0 - crs, 0.0);

        double crscuadrado = crs * crs;
        double maximocuatro = maximo * maximo * maximo * maximo;
        double maximoocho = maximocuatro * maximocuatro;
        return psi = maximoocho * (32.0 * crscuadrado * crs + 25.0 * crscuadrado + 8.0 * crs + 1.0);
    }

    static public double WndlndFnctnFirst(double r, int c)
    {
        double psifirst = 0.0;
        double crs = c * r;
        double maximo = Math.max(1.0 - crs, 0);
        double maximocubo = maximo * maximo * maximo;
        double maximosiete = maximocubo * maximocubo * maximo;
        psifirst = (-22 * c * c) * (maximosiete) * (16.0 * crs * crs + 7.0 * crs + 1.0);
        return psifirst;
    }
    static public double WndlndFnctnSecond(double r, int c)
    {
        double psisecond = 0.0;
        double crs = c * r;
        double maximo = Math.max(1.0 - crs, 0);
        double maximocubo = maximo * maximo * maximo;
        double maximoseis = maximocubo * maximocubo;
        double ccuadrado = c * c;
        psisecond = (528 * ccuadrado * ccuadrado) * (maximoseis) * (6.0 * crs + 1.0);
        return psisecond;
    }
    static public double WndlndFnctnThird(double r, int c)
    {
        //To Be CHECKED!!!
        double psithird = 0.0;
        double crs = c * r;
        double maximo = Math.max(1.0 - crs, 0);
        double maximocubo = maximo * maximo * maximo;
        double maximoseis = maximocubo * maximocubo;
        double ccubo = c * c * c;
        double cseis = ccubo * ccubo;
        psithird = (-22176.0 * cseis) * (Math.sqrt(0.5 * (Math.max(1.0 - crs, 0))));

        return psithird;
    }
}

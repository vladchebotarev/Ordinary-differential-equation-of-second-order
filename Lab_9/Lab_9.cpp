// Lab_9.cpp: определяет точку входа для консольного приложения.
//


#include "stdafx.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include <Windows.h>


using namespace std;

const double p = 1.0, q = 0.0, r = -1.0 / 4.0, s = 0.0;
const double alfa = 0.0, beta = 1.0, gamma = 0.0;  //2-obojetny 1-nie ma pochodnej
const double fi = 0.0, psi = 1.0, teta = 0.0;
const double x_poczatek = 0.0, x_koniec = 2.0;

double wzor_analityczny(double x)
{
	return (exp(x / 2.0) - exp(2.0 - x / 2.0)) / (1.0 - exp(2.0));
}

void algorytm_thomasa(double *l, double *d, double *u, double *b, double *x, int N)
{
	double *temp_b = new double[N];
	double *temp_d = new double[N];

	temp_d[0] = d[0];
	temp_b[0] = b[0];

	for (int i = 1; i < N; i++)
	{
		temp_d[i] = d[i] - l[i - 1] * (u[i - 1] / temp_d[i - 1]);
	}

	for (int i = 1; i < N; i++)
		temp_b[i] = b[i] - l[i - 1] * temp_b[i - 1] / temp_d[i - 1];

	x[N - 1] = temp_b[N - 1] / temp_d[N - 1];

	for (int i = N - 2; i >= 0; i--)
	{
		x[i] = (temp_b[i] - u[i] * x[i + 1]) / temp_d[i];
	}

	delete[] temp_b;
	delete[] temp_d;
}

double max_blad(double *blad, int N)
{
	double max = fabs(blad[0]);

	for (int i = 0; i < N; i++)
	if (fabs(blad[i]) > max)
		max = fabs(blad[i]);

	return max;
}

double dyskretyzacja_Numerowa(double h, int N)
{
	double *l, *d, *u, *b, *x, *blad, xpb = x_poczatek, xp = x_poczatek;

	l = new double[N];
	d = new double[N];
	u = new double[N];
	b = new double[N];
	x = new double[N];
	blad = new double[N];

	u[0] = alfa / h;
	d[0] = beta - alfa / h;
	b[0] = -gamma;

	for (int i = 1; i < N - 1; i++)
	{
		l[i - 1] = 1.0 / (h * h) + (-1.0 / 4.0) / 12.0;
		d[i] = (-2.0) / (h *h) + (-1.0 / 4.0) * 10.0 / 12.0;
		u[i] = 1.0 / (h * h) + (-1.0 / 4.0) / 12.0;
		b[i] = -s;
	}

	l[N - 2] = -fi / h;
	d[N - 1] = -fi / h + psi;
	b[N - 1] = -teta;

	algorytm_thomasa(l, d, u, b, x, N);

	for (int i = 0; i < N; i++)
	{
		blad[i] = fabs(x[i] - wzor_analityczny(xpb));
		xpb += h;
	}

	if (N == 20)
	{
		fstream plik;
		plik.open("wyniki_numerowa.txt", fstream::out);

		cout << endl << "Dyskretyzacja_Numerowa punkty" << endl;
		cout << "  i |         punkt |          x[n] |           U(x)| " << endl;
		cout << "-----------------------------------------------------" << endl;

		for (int i = 0; i < N; i++)
		{
			plik << xp << " " << x[i] << " " << wzor_analityczny(xp) << " " << endl;
			cout.width(4);
			cout << i << "|";
			cout.width(15);
			cout << xp << "|";
			cout.width(15);
			cout << x[i] << "|";
			cout.width(15);
			cout << wzor_analityczny(xp) << "|" << endl;

			xp += h;

		}

		plik.close();
	}

	delete[] l;
	delete[] d;
	delete[] u;
	delete[] x;
	delete[] b;


	return max_blad(blad, N);
}

double dyskretyzacja_konwencjonalna_trzypunktowa(double h, int N)
{
	double *l, *d, *u, *b, *x, *blad, xpb = x_poczatek, xp = x_poczatek;

	l = new double[N];
	d = new double[N];
	u = new double[N];
	b = new double[N];
	x = new double[N];
	blad = new double[N];

	u[0] = alfa / h;
	d[0] = beta - alfa / h;
	b[0] = -gamma;

	for (int i = 1; i < N - 1; i++)
	{
		l[i - 1] = p / (h * h) - q / (2.0 * h);;
		d[i] = (-2.0 * p) / (h *h) + r;
		u[i] = p / (h * h) - q / (2.0 * h);;
		b[i] = -s;
	}

	l[N - 2] = -fi / h;
	d[N - 1] = -fi / h + psi;
	b[N - 1] = -teta;

	algorytm_thomasa(l, d, u, b, x, N);

	for (int i = 0; i < N; i++)
	{
		blad[i] = fabs(x[i] - wzor_analityczny(xpb));
		xpb += h;
	}

	if (N == 20)
	{
		fstream plik;
		plik.open("wyniki_konwencjonalna.txt", fstream::out);

		cout << endl << "Dyskretyzacja Konwencjonalna Trzypunktowa punkty" << endl;
		cout << "  i |         punkt |          x[n] |           U(x)| " << endl;
		cout << "-----------------------------------------------------" << endl;

		for (int i = 0; i < N; i++)
		{
			blad[i] = fabs(x[i] - wzor_analityczny(xp));
			plik << xp << " " << x[i] << " " << wzor_analityczny(xp) << " " << endl;
			cout.width(4);
			cout << i << "|";
			cout.width(15);
			cout << xp << "|";
			cout.width(15);
			cout << x[i] << "|";
			cout.width(15);
			cout << wzor_analityczny(xp) << "|" << endl;

			xp += h;

		}

		plik.close();
	}

	delete[] l;
	delete[] d;
	delete[] u;
	delete[] x;
	delete[] b;

	return max_blad(blad, N);

}
/*
void dyskretyzacja_Numerowa_rys( double x_poczatek, double h, int N )
{


double *l, *d, *u, *b, *x, *blad;

l = new double[N];
d = new double[N];
u = new double[N];
b = new double[N];
x = new double[N];
blad = new double[N];

u[0] = 0.0;
d[0] = 1.0;
b[0] = 0.0;

for ( int i = 1; i < N - 1; i++ )
{
l[i - 1] = 1.0 + (h * h) / 12.0;
d[i] = - 2.0 * ( 1 - (5.0 *  h * h ) / 12.0 );
u[i] = 1.0 + ( h * h ) / 12.0;
b[i] = -( 2.0 * sin( ( i + 1 ) * h ) + 10.0 * 2.0 * sin( i * h ) + 2.0 * sin( ( i - 1 ) * h ) ) * ( h * h ) / 12.0;
}

l[N - 2] = 0.0;
d[N - 1] = 1.0;
u[N - 1] = 0.0;
b[N - 1] = -M_PI;

algorytm_thomasa( l, d, u, b, x, N );



delete[] l;
delete[] d;
delete[] u;
delete[] x;
delete[] b;
delete[] blad;

//return max_blad(blad);
}

void dyskretyzacja_konwencjonalna_trzypunktowa_rys( double x_poczatek, double h, int N)
{
fstream plik;
string nazwa = "wyniki_konwencjonalna.txt";

double *l, *d, *u, *b, *x, *blad;

l = new double[N];
d = new double[N];
u = new double[N];
b = new double[N];
x = new double[N];
blad = new double[N];

u[0] = 0.0;
d[0] = 1.0;
b[0] = 0.0;

for ( int i = 1; i < N - 1; i++ )
{
l[i - 1] = 1.0 / ( h * h );
d[i] = (-2.0 / (h * h )) + 1.0;
u[i] = 1.0 / ( h * h );
b[i] = - (2 * sin( i * h));
}

l[N - 2] = 0.0;
d[N - 1] = 1.0;
u[N - 1] = 0.0;
b[N - 1] =  -M_PI;

algorytm_thomasa( l, d, u, b, x, N );

plik.open( nazwa.c_str( ), fstream::out );

cout << endl << "Dyskretyzacja_konwencjonalna_trzypunktowa" << endl;
cout << "  i |         punkt |          x[n] |           U(x)|          blad |" << endl;
cout << "-------------------------------------------------------------------" << endl;
for ( int i = 0; i < N; i++ )
{
blad[i] = fabs( x[i] - wzor_analityczny( x_poczatek ) );

plik << x_poczatek  << " " << x[i] << " " << wzor_analityczny( x_poczatek ) << " " <<  blad[i]  << endl;
cout.width( 4 );
cout << i << "|";
cout.width( 15 );
cout << x_poczatek << "|";
cout.width( 15 );
cout << x[i] << "|";
cout.width( 15 );
cout << wzor_analityczny( x_poczatek ) << "|";
cout.width( 15 );
cout << blad[i] << "|" << endl;
//cout << x_poczatek << " " << x[i] << " " << wzor_analityczny( x_poczatek ) << " " << blad[i] << endl;
x_poczatek += h;
}

plik.close( );

delete[] l;
delete[] d;
delete[] u;
delete[] x;
delete[] b;
delete[] blad;

}

*/
void dyskretyzacja()
{
	double h, dkt, dn;
	int N;

	fstream bledy, dane;


	bledy.open("bledy.txt", fstream::out);

	//cout << endl << "\Bledy" << endl;
	//cout << "  n |         h |   trzypunktowa |        Numerow|" << endl;
	//cout << "-----------------------------------------------------" << endl;


	for (N = 10; N < 20000; N += 10)
	{
		h = (x_koniec - x_poczatek) / (N - 1);
		dkt = log10(dyskretyzacja_konwencjonalna_trzypunktowa(h, N));
		dn = log10(dyskretyzacja_Numerowa(h, N));

		bledy << log10(h) << " " << dkt << " " << dn << endl;

		/*cout.width( 4 );
		cout << N << "|";
		cout.width( 15 );
		cout << h  << "|";
		cout.width( 15 );
		cout << dkt << "|";
		cout.width( 15 );
		cout << dn << "|" << endl;*/
	}

	bledy.close();

}

int _tmain(int argc, _TCHAR* argv[])
{
	dyskretyzacja();
	system("Pause");
	return 0;
}



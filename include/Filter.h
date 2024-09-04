#ifndef FILTER_H_
#define FILTER_H_
#include "tool.h"
inline void movingAverage2D(double** grad, int nx, int ny, int wx, int wy, int times)
{
	//expand
	double** grad_ex;
	int bond = 30;
	int nx1 = nx + 2 * bond;
	int ny1 = ny + 2 * bond;
	grad_ex = alloc_double_2d(nx1, ny1);
	exmodel(grad, grad_ex, nx, ny, bond);
	for (int it = 0; it < times; it++)
	{
		for (int i = 0; i < nx1; ++i)
		{
			for (int j = 0; j < ny1; ++j)
			{
				double sum = 0.0;
				int count = 0;
				for (int ii = max(0, i - wx / 2); ii <= min(nx1 - 1, i + wx / 2); ++ii)
				{
					for (int jj = max(0, j - wy / 2); jj <= min(ny1 - 1, j + wy / 2); ++jj)
					{
						sum += grad_ex[ii][jj];
						count++;
					}
				}
				grad_ex[i][j] = sum / count;
			}
		}
	}
	for (int ix = bond; ix < nx + bond; ix++)
		for (int iy = bond; iy < ny + bond; iy++)
		{
			grad[ix - bond][iy - bond] = grad_ex[ix][iy];
		}
	free_2d(grad_ex, nx1);
}
inline void Gauss2D(double** mod, int nx, int nz, int h_size, double sigma)
{
	int bond, ix, iz, i, j;
	bond = (h_size - 1) / 2;
	double** exmod;
	exmod = alloc_double_2d(nx + 2 * bond, nz + 2 * bond);
	exmodel(mod, exmod, nx, nz, bond);
	double** gau, tmp;
	gau = alloc_double_2d(h_size, h_size);
	MakeGauss(gau, h_size, sigma);
	for (ix = bond; ix < nx + bond; ix++)
		for (iz = bond; iz < nz + bond; iz++)
		{
			tmp = 0;
			for (i = 0; i < h_size; i++)
				for (j = 0; j < h_size; j++)
					tmp += gau[i][j] * exmod[ix - bond + i][iz - bond + j];
			mod[ix - bond][iz - bond] = tmp;
		}
	free_2d(exmod, nx + 2 * bond);
	free_2d(gau, h_size);
}
inline void smoothgauss2dtopo(double** v, int nz, int nx, double dz, double freq, vector<int>& itopo, int maxn, double fracx, double fracz, double vref, int nptt)
{
	double* beta1, * beta2, ** vf;
	double wl, taux, tauz, taux2, tauz2, h, xl2, xl1, d, betatot;
	int ix, iz, il2, il1, k1, k2, j1, j2, jj1, jj2;

	beta1 = new double[maxn];
	beta2 = new double[maxn];
	vf = new double* [nx];
	for (ix = 0; ix < nx; ix++)vf[ix] = new double[nz];

	// Correlation lengths of Gaussian filter are defined as fraction of the wavelength 
	// as defined by the inverted frequency component and a reference velocity

	wl = vref / freq;
	taux = fracx * wl;
	tauz = fracz * wl;
	h = dz;
	if (taux == 0 && tauz == 0)
		for (ix = 0; ix < nx; ix++)
			for (iz = 0; iz < nz; iz++)
				vf[ix][iz] = v[ix][iz];
	else
	{
		for (ix = 0; ix < nx; ix++)
			for (iz = 0; iz < itopo[ix] + nptt; iz++)
				vf[ix][iz] = v[ix][iz];
		xl2 = 3. * taux;
		xl1 = 3. * tauz;
		il2 = int(xl2 / h) + 1;
		il1 = int(xl1 / h + 1);
		if (2 * il1 + 1 > maxn)
		{
			cout << "*******increase maxn in smoothgauss*******" << endl;
			exit(0);
		}
		if (2 * il2 + 1 > maxn)
		{
			cout << "*******increase maxn in smoothgauss*******" << endl;
			exit(0);
		}
		tauz2 = tauz * tauz;
		taux2 = taux * taux;
		k2 = 0;
		for (ix = -il2; ix <= il2; ix++)
		{
			d = ix * h;
			beta2[k2] = exp(-d * d / taux2);
			k2++;
		}
		k1 = 0;
		for (iz = -il1; iz <= il1; iz++)
		{
			d = iz * h;
			beta1[k1] = exp(-d * d / tauz2);
			k1++;
		}

		for (ix = 0; ix < nx; ix++)
			for (iz = itopo[ix] + nptt; iz < nz; iz++)
			{
				vf[ix][iz] = 0;
				betatot = 0;
				jj2 = -1;
				for (j2 = -il2; j2 <= il2; j2++)
				{
					k2 = j2 + ix;
					jj2++;
					if (k2 < 0 || k2 >= nx)continue;
					jj1 = -1;
					for (j1 = -il1; j1 <= il1; j1++)
					{
						k1 = j1 + iz;
						jj1++;
						if (k1 < itopo[ix] + nptt || k1 >= nz)continue;
						vf[ix][iz] += beta1[jj1] * beta2[jj2] * v[k2][k1];
						betatot += beta1[jj1] * beta2[jj2];
					}
				}
				vf[ix][iz] /= betatot;
			}
	}

	for (ix = 0; ix < nx; ix++)
		for (iz = 0; iz < nz; iz++)
			v[ix][iz] = vf[ix][iz];

	delete[] beta1;
	delete[] beta2;
	for (ix = 0; ix < nx; ix++)delete[] vf[ix];
	delete[] vf;
}
#endif
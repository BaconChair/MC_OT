// multi-specular-k.cpp : 定义控制台应用程序的入口点。
//

// multi_trywtire.cpp : 定义控制台应用程序的入口点。
//
#include "stdafx.h"
#include<stdio.h> 
#include <stdlib.h>
#include<time.h>
#include <math.h>



int main()
{
	srand(time(NULL));


	int mi, mr;
	int N = 2;//层数
	double kk[19] = { 0 };
	for (int h = 0; h <= 9; h++) { kk[h] = 0.0001*(h+1); }
	for (int h = 0; h <= 9; h++) { kk[9 + h] = 0.001*(h + 1); }
	double k[5] = { 0,0.001,0.001,0.001,0.001 };//N+1
	double n[5] = { 1.0,1.5,1.7,1.7,2.0 };//N+1
	double r[4] = { 50,30,40.5,10 };//N
	double lam[28];//波长
	for (int h = 0; h <= 9; h++) { lam[h] = 314.1592653 / (h + 1); }
	for (int h = 0; h <= 9; h++) { lam[9 + h] = 31.41592653 / (h + 1); }
	for (int h = 0; h <= 9; h++) { lam[18 + h] = 3.141592653 / (h + 1); }
	for (int kkj = 18; kkj <= 18; kkj++)
	{
		
		//k[2] = kk[kkj];
		

		for (int h = 18; h <= 18; h++)
		{
			double M_PI = 3.141592653;//π
			int num[181] = { 0 };//分为180份，计数
			int count, out = 0, absorb = 0;
			int not_in = 0;
			for (count = 0; count < 1000000; count++)
			{ //0

				double rou;
				double fai;//散射角
				double lv;//光程
				double pp[3];//入射点坐标
				double nor[3] = { 0,0,0 };//法线
				double ray[3] = { 0,0,-1 };//方向向量
				double nri = 0;
				double x, y, z;
				double refl(double incray[3], double normal[3], double nri);
				void reflray(double ray[3], double normal[3]); void refrray(double ray[3], double normal[3], double nri);
				int tracing(double pp[3], double ray[3], int mi, double *lv, int N, double br, double sr);
				double angle(double Ma[3]); void scatter(double ang, int num[181]);

				double Rx = rand() / (RAND_MAX + 1.0), Ry = rand() / (RAND_MAX + 1.0);
				x = 2 * r[0] * Rx - r[0]; y = 2 * r[0] * Ry - r[0];
				if (x*x + y*y > r[0] * r[0]) { not_in++; }
				else//与球相遇
				{//1
					pp[0] = x; pp[1] = y; pp[2] = sqrt(r[0] * r[0] - x*x - y*y);
					for (int i1 = 0; i1 < 3; i1++)
					{
						nor[i1] = pp[i1] / r[0];
					}//noraml

					nri = n[1] / n[0];
					rou = refl(ray, nor, nri);
					double R1 = rand() / (RAND_MAX + 1.0);
					if (R1 <= rou)
					{//2
						out++;
						reflray(ray, nor);
						fai = angle(ray); scatter(fai, num);

					}//2
					else//enter into sphere
					{//2_1
						mi = 2;
						nri = n[1] / n[0];
						refrray(ray, nor, nri);
						while (mi > 1)
						{//2_2

							mr = tracing(pp, ray, mi, &lv, N, r[mi - 2], r[mi - 1]);
							//	printf("mr2=%d\n", mr);
							if (mi < (N + 1))//不在最内球
							{

								double Ro = rand() / (RAND_MAX + 1.0);
								double beta = 4 * M_PI*k[mi - 1] / lam[h];
								double s = -log(Ro) / beta;
								if (s < lv) { absorb++; mi = 1; }
								else
								{//2_5
									if (mr == (mi - 1))//如果下一层介质向外
									{//2_6
										for (int i1 = 0; i1 < 3; i1++) { nor[i1] = -pp[i1] / r[mi - 2]; }
										//double a = nor[0] * nor[0] + nor[1] * nor[1] + nor[2] * nor[2];
										//printf("a1=%f\t", a); printf("r=%f\n", r[mi - 2]);
										nri = n[mi - 2] / n[mi - 1];
										rou = refl(ray, nor, nri);
										double R1 = rand() / (RAND_MAX + 1.0);
										if (R1 <= rou)
										{
											reflray(ray, nor);

										}
										else
										{
											mi = mr;
											refrray(ray, nor, nri);
											if (mr == 1)
											{
												out++;
												fai = angle(ray); scatter(fai, num);
												mi = 1;
											}
										}
									}//2_6
									else if (mr == (mi + 1))//如果下一层介质向内
									{//2_6
										for (int i1 = 0; i1 < 3; i1++) { nor[i1] = pp[i1] / r[mi - 1]; }

										//double a = nor[0] * nor[0] + nor[1] * nor[1] + nor[2] * nor[2];
										//printf("a2=%f\n", a);
										nri = n[mi] / n[mi - 1];
										rou = refl(ray, nor, nri);
										double R1 = rand() / (RAND_MAX + 1.0);
										if (R1 <= rou)
										{
											reflray(ray, nor);
										}
										else
										{
											mi = mr;
											refrray(ray, nor, nri);
										}
									}//2_6
								}//2_5

							}
							else if (mi == (N + 1))//2_4在最内球
							{
								double Ro = rand() / (RAND_MAX + 1.0);
								double beta = 4 * M_PI*k[mi - 1] / lam[h];
								double s = -log(Ro) / beta;
								if (s < lv) { absorb++; mi = 1; }
								else {
									for (int i1 = 0; i1 < 3; i1++) { nor[i1] = -pp[i1] / r[mi - 2]; }

									//double a = nor[0] * nor[0] + nor[1] * nor[1] + nor[2] * nor[2];
									//printf("a3=%f\n", a);
									nri = n[mi - 2] / n[mi - 1];
									rou = refl(ray, nor, nri);
									double R1 = rand() / (RAND_MAX + 1.0);
									if (R1 <= rou)
									{
										reflray(ray, nor);
									}
									else
									{
										mi = mr;
										refrray(ray, nor, nri);
									}


								}//2_4
							}
						}//2_2

					}//2_1

				}//1
			}//0

			int total = count - not_in;
			double alpha = 1.0*absorb / total;
			double chi = 2 * 3.14159*r[0] / lam[h];


			//printf("%f\t", chi);
			//printf("%.2f\t", n2);
			//printf("%f\n", alpha);
			//		printf("吸收条数=%d\t", absorb);
			//		printf("散射条数=%d\t", out);
			//		printf("%总光线数=%d\n", total);

					float phi1 = 2 * num[0] / sin(0.25*M_PI / 180) / (0.5*M_PI / 180) / out;
			//printf("%d\t%f\n", num[0], phi1);
				printf("%f\n", phi1);
					for (int u = 1; u < 180; u++)
				{
						float phi = 2 * num[u] / sin(u*M_PI / 180) / (1 * M_PI / 180) / out;
			//			//printf("%d\t%f\n", num[u], phi);
						printf("%f\n", phi);
					}
					float phi2 = 2 * num[180] / sin(179.75*M_PI / 180) / (0.5 * M_PI / 180) / out;
			//printf("%d\t%f\n", num[180], phi2);
					printf("%f\n", phi2);
		}
	}
	return 0;

}

double refl(double incray[3], double normal[3], double nri)//反射率函数
{
	double thetai, thetar, a, b, R;
	a = -incray[0] * normal[0] - incray[1] * normal[1] - incray[2] * normal[2];

	thetai = acos(a);
	if (sin(thetai) >= nri)
	{
		R = 1;
	}
	else {
		thetar = asin(sin(thetai) / nri);
		R = (pow((nri*cos(thetai) - cos(thetar)), 2) / pow((nri*cos(thetai) + cos(thetar)), 2) + pow((cos(thetai) - nri*cos(thetar)), 2) / pow((cos(thetai) + nri*cos(thetar)), 2)) / 2;
		if (R > 1)
		{
			R = 1;
		}
	}
	return R;
}

void reflray(double ray[3], double normal[3])//反射方向
{
	double ray1[3];
	ray1[0] = ray[0] - 2 * (normal[0] * ray[0] + normal[1] * ray[1] + normal[2] * ray[2])*normal[0];
	ray1[1] = ray[1] - 2 * (normal[0] * ray[0] + normal[1] * ray[1] + normal[2] * ray[2])*normal[1];
	ray1[2] = ray[2] - 2 * (normal[0] * ray[0] + normal[1] * ray[1] + normal[2] * ray[2])*normal[2];
	ray[0] = ray1[0] / sqrt(ray1[0] * ray1[0] + ray1[1] * ray1[1] + ray1[2] * ray1[2]);
	ray[1] = ray1[1] / sqrt(ray1[0] * ray1[0] + ray1[1] * ray1[1] + ray1[2] * ray1[2]);
	ray[2] = ray1[2] / sqrt(ray1[0] * ray1[0] + ray1[1] * ray1[1] + ray1[2] * ray1[2]);
}

void refrray(double ray[3], double normal[3], double nri)//折射方向
{
	double ray2[3];
	double thetai, thetar, a, b;
	a = -ray[0] * normal[0] - ray[1] * normal[1] - ray[2] * normal[2];

	thetai = acos(a);
	thetar = asin(sin(thetai) / nri);
	//ray2[0] = nri*ray[0] + (nri*cos(thetai) - cos(thetar))*normal[0];
	//ray2[1] = nri*ray[1] + (nri*cos(thetai) - cos(thetar))*normal[1];
	//ray2[2] = nri*ray[2] + (nri*cos(thetai) - cos(thetar))*normal[2];
	ray2[0] = ray[0] + (cos(thetai) - cos(thetar)*nri)*normal[0];
	ray2[1] = ray[1] + (cos(thetai) - cos(thetar)*nri)*normal[1];
	ray2[2] = ray[2] + (cos(thetai) - cos(thetar)*nri)*normal[2];
	ray[0] = ray2[0] / sqrt(ray2[0] * ray2[0] + ray2[1] * ray2[1] + ray2[2] * ray2[2]);
	ray[1] = ray2[1] / sqrt(ray2[0] * ray2[0] + ray2[1] * ray2[1] + ray2[2] * ray2[2]);
	ray[2] = ray2[2] / sqrt(ray2[0] * ray2[0] + ray2[1] * ray2[1] + ray2[2] * ray2[2]);

}
int tracing(double pp[3], double ray[3], int mi, double *lv, int N, double br, double sr)//追踪函数
{
	int mr;
	double bb, cc, tt;
	double ratio = sr / br;
	double pp0;
	if (mi < (N + 1))
	{
		pp0 = (pp[0] * pp[0] + pp[1] * pp[1] + pp[2] * pp[2]) / (br*br);
		if (pp0 >((1 + ratio*ratio) / 2))
		{
			bb = ray[0] * pp[0] + ray[1] * pp[1] + ray[2] * pp[2];
			cc = (pp[0] * pp[0] + pp[1] * pp[1] + pp[2] * pp[2]) - sr*sr;
			if ((bb*bb - cc) > 0)
			{
				tt = -sqrt(bb*bb - cc) - bb;
				mr = mi + 1;
				//		printf("tt=%f\n", tt);
				//		printf("mr=%d\n", mr);
			}
			else { tt = fabs(2 * bb); mr = mi - 1; }
			*lv = fabs(tt);
			pp[0] = pp[0] + ray[0] * tt;
			pp[1] = pp[1] + ray[1] * tt;
			pp[2] = pp[2] + ray[2] * tt;

			//	printf("mr1=%d\n", mr);

		}
		else
		{
			bb = ray[0] * pp[0] + ray[1] * pp[1] + ray[2] * pp[2];
			cc = (pp[0] * pp[0] + pp[1] * pp[1] + pp[2] * pp[2]) - br*br;
			tt = fabs(sqrt(fabs(bb*bb - cc)) - bb);
			*lv = fabs(tt);
			pp[0] = pp[0] + ray[0] * tt;
			pp[1] = pp[1] + ray[1] * tt;
			pp[2] = pp[2] + ray[2] * tt;
			mr = mi - 1;
		}
	}
	else if (mi == (N + 1))
	{
		bb = ray[0] * pp[0] + ray[1] * pp[1] + ray[2] * pp[2];
		tt = fabs(2 * bb);
		*lv = fabs(tt);
		pp[0] = pp[0] + ray[0] * tt;
		pp[1] = pp[1] + ray[1] * tt;
		pp[2] = pp[2] + ray[2] * tt;
		mr = mi - 1;

	}
	return mr;
}

double angle(double Ma[3])//求散射角
{
	double a2 = -Ma[2];

	double sita = acos(a2);
	return sita;


}

void scatter(double sita, int num[181])//散射光线计数
{

	double ang = sita;
	if (ang - (0.5*3.1415926 / 180)>0)
	{
		int j = (ang - (0.5*3.1415926 / 180)) / (3.1415926 / 180);
		num[j + 1]++;
	}
	else { num[0]++; }


}

               
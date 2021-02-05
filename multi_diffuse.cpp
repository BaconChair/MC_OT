// multi_diffuse.cpp : 定义控制台应用程序的入口点。
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
	int N = 4;//层数
	double k[5] = { 0,0.001,0.001,0.001,0.001 };//N+1
	double n[5] = { 1.0,1.3,1.5,1.7,1.5 };//N+1
	double rou0[10] = { 0 }, rou1[10] = { 0 };
	void refl(double rou0[10], double rou1[10],double n[5],int N);
	refl(rou0, rou1, n, N);//每个界面的两侧的反射率
	//printf("rou0=%f\trou1=%f\t", rou0[0], rou1[1]);
	double r[4] = { 50,45,40.5,36.45 };//N
	double lam[28];//波长
	for (int h = 0; h <= 9; h++) { lam[h] = 314.1592653 / (h + 1); }
	for (int h = 0; h <= 9; h++) { lam[9 + h] = 31.41592653 / (h + 1); }
	for (int h = 0; h <= 9; h++) { lam[18 + h] = 3.141592653 / (h + 1); }

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
			
			void reflray(double normal[3],double ray[3]); void refrray(double normal[3], double ray[3]);
			int tracing(double pp[3], double ray[3], int mi, double *lv, int N, double br, double sr);
			double angle(double Ma[3]); void scatter(double ang, int num[181]);

			double Rx = rand() / (RAND_MAX + 1.0), Ry = rand() / (RAND_MAX + 1.0);
			x = 2 * r[0] * Rx - r[0]; y = 2 * r[0] * Ry - r[0];
			if (x*x + y*y > r[0] * r[0]) { not_in++; }
			else//与球相遇
			{//1
				pp[0] = x; pp[1] = y; pp[2] = sqrt(r[0] * r[0] - x*x - y*y);
				
				for (int i1 = 0; i1 < 3; i1++) { nor[i1] = pp[i1] / r[0]; }
				
				
				double R1 = rand() / (RAND_MAX + 1.0);
				if (R1 <= rou0[0])
				{//2
					out++;
					reflray(nor,ray);
					//printf("ray1=%f\t%f\t%f\n", ray[0],ray[1],ray[2]);
					fai = angle(ray); 
					//printf("ang1=%f\n", fai);
					scatter(fai, num);

				}//2
				else//enter into sphere
				{//2_1
					mi = 2;
					nri = n[1] / n[0];
					refrray(nor,ray);
					while (mi > 1)
					{//2_2

						mr = tracing(pp, ray, mi, &lv, N, r[mi - 2], r[mi - 1]);
							//printf("mi=%d\t", mi);
							//printf("mr=%d\t", mr);
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
									//printf("a1=%f\t", a); printf("r1=%f\n", r[mi - 2]);
									//nri = n[mi - 2] / n[mi - 1];
									double R1 = rand() / (RAND_MAX + 1.0);
									if (R1 <= rou1[mi-2])
									{
										reflray(nor,ray);

									}
									else
									{
										mi = mr;
										refrray(nor,ray);
										//printf("ray2=%f\t%f\t%f\n", ray[0], ray[1], ray[2]);
										if (mr == 1)
										{
											out++;
											fai = angle(ray); 
											//printf("ang2=%f\n", fai);
											scatter(fai, num);
											mi = 1;
										}
									}
								}//2_6
								else if (mr == (mi + 1))//如果下一层介质向内
								{//2_6
									for (int i1 = 0; i1 < 3; i1++) { nor[i1] = pp[i1] / r[mi - 1]; }

									//double a = nor[0] * nor[0] + nor[1] * nor[1] + nor[2] * nor[2];
									//printf("a2=%f\t", a); 
									//printf("r2=%f\n", r[mi - 1]);
									//nri = n[mi] / n[mi - 1];
									double R1 = rand() / (RAND_MAX + 1.0);
									if (R1 <= rou0[mi-1])
									{
										reflray(nor,ray);
									}
									else
									{
										mi = mr;
										refrray(nor,ray);
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
								//printf("a3=%f\t", a); printf("r3=%f\n", r[mi - 2]);
								//nri = n[mi - 2] / n[mi - 1];
								
								double R1 = rand() / (RAND_MAX + 1.0);
								if (R1 <= rou1[mi-2])
								{
									reflray(nor,ray);
								}
								else
								{
									mi = mr;
									refrray(nor,ray);
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
	return 0;

}

void refl(double rou0[10], double rou1[10],double n[5],int N)//反射率函数
{
	for (int i = 0; i < N;i++)
	{
		
		if (n[i] < n[i+1])
		{
			for (float j = 0; j < 1.57; j += 0.001)
			{
				float f1 = pow((1 - pow(sin(j)*n[i] / n[i+1], 2)), 0.5);

				float rou2 = pow((n[i]*f1 - n[i+1]*cos(j)) / (n[i]*f1 + n[i+1]*cos(j)), 2) + pow((n[i+1]*f1 - n[i]*cos(j)) / (n[i+1]*f1 + n[i]*cos(j)), 2);
				rou0[i] = rou2*sin(j)*cos(j)*0.001 + rou0[i];
			}
			rou1[i] = rou0[i] *n[i]*n[i]/ (n[i+1]*n[i+1]) + (1 - n[i] * n[i] / (n[i + 1] * n[i + 1]));
		}
		else
		{
			for (float j = 0; j < 1.57; j += 0.001)
			{
				float f1 = pow((1 - pow(sin(j)*n[i+1] / n[i], 2)), 0.5);

				float rou2 = pow((n[i] * f1 - n[i + 1] * cos(j)) / (n[i] * f1 + n[i + 1] * cos(j)), 2) + pow((n[i + 1] * f1 - n[i] * cos(j)) / (n[i + 1] * f1 + n[i] * cos(j)), 2);
				rou1[i] = rou2*sin(j)*cos(j)*0.001 + rou1[i];
			}
			rou0[i] = rou1[i]*n[i+1]*n[i+1] / (n[i]*n[i]) + (1 - n[i+1]*n[i+1] / (n[i]*n[i]));
		}
	}
	
}

void reflray(double normal[3],double ray[3])//反射方向
{
	double xx[3],yy[3],zz[3];
	zz[0] = normal[0]; zz[1] = normal[1]; zz[2] = normal[2];
	if (zz[0] * zz[0] + zz[1] * zz[1] == 0) { xx[0] = 1; xx[1] = 0; xx[2] = 0; }
	else
	{
		xx[0] = zz[1]; xx[1] = -zz[0]; xx[2] = 0;
		double a = sqrt(xx[0] * xx[0] + xx[1] * xx[1] + xx[2] * xx[2]);
		xx[0] = xx[0] / a; xx[1] = xx[1] / a; xx[2] = xx[2] / a;
	}
	yy[0] = zz[1] * xx[2] - zz[2] * xx[1]; 
	yy[1] = zz[2] * xx[0] - zz[0] * xx[2];
	yy[2] = zz[0] * xx[1] - zz[1] * xx[0];
	double b = sqrt(yy[0] * yy[0] + yy[1] * yy[1] +yy[2] * yy[2]);
	yy[0] = yy[0] / b; yy[1] = yy[1] / b; yy[2] = yy[2] / b;
	double Ro = rand() / (RAND_MAX + 1.0); double Rp = rand() / (RAND_MAX + 1.0);
	double theta1 = acos(sqrt(1 - Ro)), phi1 = 2 * 3.141592653*Rp;
	double i1 = sin(theta1)*cos(phi1), j1 = sin(theta1)*sin(phi1), k1 = cos(theta1);
	ray[0] = xx[0] * i1 + yy[0] * j1 + zz[0] * k1;
	ray[1] = xx[1] * i1 + yy[1] * j1 + zz[1] * k1;
	ray[2] = xx[2] * i1 + yy[2] * j1 + zz[2] * k1;
	double c = sqrt(ray[0] * ray[0] + ray[1] * ray[1] + ray[2] * ray[2]);
	ray[0] = ray[0] / c; ray[1] = ray[1] / c; ray[2] = ray[2] / c;
	//printf("ray3=%f\t%f\t%f\n", ray[0], ray[1], ray[2]);

}

void refrray(double normal[3], double ray[3])//zhe射方向
{
	double xx[3], yy[3], zz[3];
	zz[0] = normal[0]; zz[1] = normal[1]; zz[2] = normal[2];
	if (zz[0] * zz[0] + zz[1] * zz[1] == 0) { xx[0] = 1; xx[1] = 0; xx[2] = 0; }
	else
	{
		xx[0] = zz[1]; xx[1] = -zz[0]; xx[2] = 0;
		double a = sqrt(xx[0] * xx[0] + xx[1] * xx[1] + xx[2] * xx[2]);
		xx[0] = xx[0] / a; xx[1] = xx[1] / a; xx[2] = xx[2] / a;
	}
	yy[0] = zz[1] * xx[2] - zz[2] * xx[1];
	yy[1] = zz[2] * xx[0] - zz[0] * xx[2];
	yy[2] = zz[0] * xx[1] - zz[1] * xx[0];
	double b = sqrt(yy[0] * yy[0] + yy[1] * yy[1] + yy[2] * yy[2]);
	yy[0] = yy[0] / b; yy[1] = yy[1] / b; yy[2] = yy[2] / b;
	double Ro = rand() / (RAND_MAX + 1.0); double Rp = rand() / (RAND_MAX + 1.0);
	double theta1 = acos(sqrt(1 - Ro)), phi1 = 2 * 3.141592653*Rp;
	double i1 = sin(theta1)*cos(phi1), j1 = sin(theta1)*sin(phi1), k1 = cos(theta1);
	ray[0] = xx[0] * i1 + yy[0] * j1 + zz[0] * k1;
	ray[1] = xx[1] * i1 + yy[1] * j1 + zz[1] * k1;
	ray[2] = xx[2] * i1 + yy[2] * j1 + zz[2] * k1;
	double c = sqrt(ray[0] * ray[0] + ray[1] * ray[1] + ray[2] * ray[2]);
	ray[0] = -ray[0] / c; ray[1] = -ray[1] / c; ray[2] = -ray[2] / c;


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
			else 
			{ tt = fabs(2 * bb); mr = mi - 1; }
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
	else //if (mi == (N + 1))
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


}         
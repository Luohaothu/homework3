#include <iostream>
#include <stdio.h>
using namespace std;
struct datas
{
	int i, j;
	float t1, t2, t3, p, m1, m2, s, w;
};

int main()
{
	datas d[49][19];
	FILE *fp;
	fp = fopen("./out", "r");
	if (fp!= NULL)
	{
		for (int i = 0; i < 49; i++)
		{
			for (int j = 0; j < 19; j++)
			{
				fscanf(fp, "%d%d\n", d[i][j].i, d[i][j].j);
				fscanf(fp, "%f\n", d[i][j].t1);
				fscanf(fp, "%f\n", d[i][j].t2);
				fscanf(fp, "%f\n", d[i][j].t3);
				fscanf(fp, "%f\n", d[i][j].p);
				fscanf(fp, "%f\n", d[i][j].m1);
				fscanf(fp, "%f\n", d[i][j].m2);
				fscanf(fp, "%f\n", d[i][j].s);
				fscanf(fp, "%f\n", d[i][j].w);
			}
		}
		printf("%d   %d   %f\n", d[1][1].i, d[1][1].j, d[1][1].t1);
	}
	return 0;
}

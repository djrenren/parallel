//
// Example: summation
//

#include <xmtc.h>
#include <xmtio.h>

#define n 8
#define log_n 3

void print_row(int* row)
{
	int i;
	for(i=1;i<=n;i++)
	{
	   if(row[i]<10)
		   printf("  ");
	   else if(row[i]<100)
		   printf(" ");
		printf("%d",row[i]);
		printf(" ");
	}
	printf("\n");
}
int main()
{
	int h;               // index for moving bottom-to-top
	int p;               // for any h, p=2^h
	int A[n+1];          // the original data
	int B[log_n+1][n+1]; // the data at each level of the tree
	int C[log_n+1][n+1];
	A[1]=1; A[2]=1; A[3]=2; A[4]=3; A[5]=5; A[6]=2; A[7]=1; A[8]=2;

	spawn(1,n)           // for i, 1<=i<=n pardo
	{
		int i;            // index for moving left-to-right
		i=$;              // XMT-C uses dollar sign
		B[0][i]=A[i];
	}

	h=log_n;
	for(p=n;p>1;p/=2)   // move up the tree
	{
		spawn(1,n/p)
		{
			int i=$;
			if(i == 1)
				C[h][1] = B[h][1];
			else if(i%2)
				C[h][i] = C[h+1][i/2];
			else
				C[h][i] = C[h+1][(i-1)/2 * b[h][i]];
		}
		h-=1;             // h goes from 1 to log_n
	}

	printf("\n");
	print_row(A);
	printf("\n");

	for(h=log_n;h>=0;h--)
	{
		print_row(C[h]);
	}
	printf("\n");

	printf("sum=%d\n",B[log_n][1]);
	printf("\n");
}
//
// end of file
//

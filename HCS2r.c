#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define PI 3.141592654
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS (1.2E-07)
#define RNMX (1.0-EPS)

double bnldev(double pp, long n, long *idum); double gammln(double xx); double ran1(long *idum); double randn (double mu, double sigma);

int main(void){

long N1, N2, TIME, REP, PHASE; double NMU, MIGS, sel, REC, pla, epis;

FILE * par, * inp, *input, * out1; 
par = fopen("pars1", "r"); 
inp = fopen("trajectory", "r");
out1 = fopen("outst", "w");

fscanf(par, "%ld", &PHASE);
fscanf(par, "%lf", &MIGS);

REP = 2000000, N1 = 10000, N2 = 10000, TIME = 100*(N1+N2), pla = 1.0; REC = 0.5; NMU = 0.1;

epis = -0.1;

long seed = -2345;

double D1, D2, r1, r2, s1, s2;  int one, two;  long R1, R2, S1, S2; 
double avefit1, avefit2, select = 0.0, rab1, raB1, rAb1, rAB1, rab2, raB2, rAb2, rAB2;
double ab1, aB1, Ab1, AB1, ab2, aB2, Ab2, AB2, fit_ab1, fit_aB1, fit_Ab1, fit_ab2, fit_aB2, fit_Ab2;
double newN1, newN2, ranpick1, ranpick2, m,m1,m2, ab1_migs=0.0, ab1_stay, aB1_migs=0.0, aB1_stay, Ab1_migs=0.0, Ab1_stay, AB1_migs=0.0, AB1_stay, ab2_migs=0.0, ab2_stay, aB2_migs=0.0, aB2_stay, Ab2_migs=0.0, Ab2_stay, AB2_migs=0.0, AB2_stay, ab1new, Ab1new, aB1new, AB1new, ab2new, Ab2new, aB2new, AB2new;
double hl_first_locus = 0.0, hl_second_locus = 0.0; 
long cc, k, i, l;  double a,b, ax, bx, tra[3*PHASE], smax; int rstart, rstart1, xx; long long gentime = 0, cosegtime = 0;
double pop1 = (double)(N1)/(double)(N1+N2);

	for (l = 0; l < (3*PHASE); l++){
		fscanf(inp, "%lf", &tra[l]);
		}

		smax = tra[PHASE/4];

for (i = 0 ; i < REP; i++){

  ab1 = 1.000; ab2 = 1.000;
  aB1 = 0.00; Ab1 = 0.00; AB1 = 0.00;
  aB2 = 0.00; Ab2 = 0.00; AB2 = 0.00;

   one = 0; two = 0;

	  if(ran1(&seed) <= 0.5)
	  {
	  	if(ran1(&seed) < pop1)
	      {
	    	  Ab1 = 1.00/(double)N1;
	    	  ab1 = 1.00 - Ab1;
	    	}
	  	else
	      {
	    	  Ab2 = 1.00/(double)N2;
	    	  ab2 = 1.00 - Ab2;
	  	  }
	
	  }

	  else
	  {
	    if(ran1(&seed) < pop1)
	      {
	      	aB1 = 1.00/(double)N1;
	      	ab1 = 1.00-aB1;
	  	  }
	  	else
	     {
	      	aB2 = 1.00/(double)N2;
	      	ab2 = 1.00-aB2;
	  	 }

	  }

	  rstart = (int)(ran1(&seed)*PHASE);
	  rstart1 = 0; 

  for( k = 0 ; k< TIME ; k++)
  {

  	/////////////////  MIGRATION ////////////////////

  	if (MIGS < 100.0)
    {
        newN1 = (double)N1;
        newN2 = (double)N2;

        ab1_stay = ab1 * newN1 ;
        aB1_stay = aB1 * newN1 ;
        Ab1_stay = Ab1 * newN1 ;
        AB1_stay = AB1 * newN1 ;

        ab2_stay = ab2 * newN2 ;
        aB2_stay = aB2 * newN2 ;
        Ab2_stay = Ab2 * newN2 ;
        AB2_stay = AB2 * newN2 ;

        ab1_migs = 0.0; aB1_migs = 0.0; Ab1_migs = 0.0; AB1_migs = 0.0;
        ab2_migs = 0.0; aB2_migs = 0.0; Ab2_migs = 0.0; AB2_migs = 0.0;

      for (xx = 0; xx< MIGS; xx++)
      {
          ranpick1 = ran1(&seed);
          ranpick2 = ran1(&seed);

          newN1 = newN1 - 1.0;
          newN2 = newN2 - 1.0;

          ab1 = ab1_stay/newN1 ;
          aB1 = aB1_stay/newN1 ;
          Ab1 = Ab1_stay/newN1 ;

          ab2 = ab2_stay/newN2 ;
          aB2 = aB2_stay/newN2 ;
          Ab2 = Ab2_stay/newN2 ;

          if (ranpick1 <= ab1)
          {
            ab1_migs = ab1_migs + 1.0;
            ab1_stay = ab1_stay - 1.0;
          }

          else if (ranpick1 > ab1 && ranpick1 <= (ab1 + aB1))
          {
            aB1_migs = aB1_migs + 1.0;
            aB1_stay = aB1_stay - 1.0;
          }

          else if (ranpick1 > (ab1 + aB1) && ranpick1 <= (ab1 + aB1 + Ab1))
          {
            Ab1_migs = Ab1_migs + 1.0;
            Ab1_stay = Ab1_stay - 1.0;
          }

          else
          {
            AB1_migs = AB1_migs + 1.0;
            AB1_stay = AB1_stay - 1.0;
          }

          if (ranpick2 <= ab2)
          {
            ab2_migs = ab2_migs + 1.0;
            ab2_stay = ab2_stay - 1.0;
          }

          else if (ranpick2 > ab2 && ranpick2 <= (ab2 + aB2))
          {
            aB2_migs = aB2_migs + 1.0;
            aB2_stay = aB2_stay - 1.0;
          }

          else if (ranpick2 > (ab2 + aB2) && ranpick2 <= (ab2 + aB2 + Ab2))
          {
            Ab2_migs = Ab2_migs + 1.0;
            Ab2_stay = Ab2_stay - 1.0;
          }

          else
          {
            AB2_migs = AB2_migs + 1.0;
            AB2_stay = AB2_stay - 1.0;
          }

      } // end for loop under migs<100

      ab1new = (ab1_stay + ab2_migs)/(double)N1;
      aB1new = (aB1_stay + aB2_migs)/(double)N1;
      Ab1new = (Ab1_stay + Ab2_migs)/(double)N1;
      AB1new = (AB1_stay + AB2_migs)/(double)N1;

      ab2new = (ab2_stay + ab1_migs)/(double)N2;
      aB2new = (aB2_stay + aB1_migs)/(double)N2;
      Ab2new = (Ab2_stay + Ab1_migs)/(double)N2;
      AB2new = (AB2_stay + AB1_migs)/(double)N2;

      ab1 = ab1new;
      Ab1 = Ab1new;
      aB1 = aB1new;
      AB1 = AB1new;
      ab2 = ab2new;
      Ab2 = Ab2new;
      aB2 = aB2new;
      AB2 = AB2new;

    } //end if migs<100


    else if (MIGS >= 100.0)
    {
      m1 = MIGS/(double)N1;
      m2 = MIGS/(double)N2;

      ab1new = ab1*(1.0-m1)+ab2*m1;
      Ab1new = Ab1*(1.0-m1)+Ab2*m1;

      aB1new = aB1*(1.0-m1)+aB2*m1;
      AB1new = AB1*(1.0-m1)+AB2*m1;

      ab2new = ab2*(1.0-m2)+ab1*m2;
      Ab2new = Ab2*(1.0-m2)+Ab1*m2;

      aB2new = aB2*(1.0-m2)+aB1*m2;
      AB2new = AB2*(1.0-m2)+AB1*m2;


    	ab1 = ab1new;
    	Ab1 = Ab1new;
    	aB1 = aB1new;
    	AB1 = AB1new;
    	ab2 = ab2new;
    	Ab2 = Ab2new;
    	aB2 = aB2new;
    	AB2 = AB2new;

    }   // end else if

			///////////////////   SELECTION     ///////////////////
  if (rstart1 == PHASE) rstart1 = 0;

////////////// additive and epistasis except sign

      avefit1 = ab1*(1.0 - (tra[rstart1+rstart]*(1.0+epis)))*(1.0 -(tra[rstart1+rstart]*(1.0+epis))) + Ab1*(1.0 + tra[rstart1+rstart])*(1.0 - tra[rstart1+rstart]) + aB1*(1.0 - tra[rstart1+rstart])*(1.0 + tra[rstart1+rstart]) +  AB1*(1.0 + (tra[rstart1+rstart]*(1.0+epis)))*(1.0 + (tra[rstart1+rstart]*(1.0+epis)));

        fit_ab1= ((1.0 - (tra[rstart1+rstart]*(1.0+epis)))*(1.0 - (tra[rstart1+rstart]*(1.0+epis))))/avefit1;
        fit_Ab1= ((1.0 + tra[rstart1+rstart])*(1.0 - tra[rstart1+rstart]))/avefit1;
        fit_aB1= ((1.0 - tra[rstart1+rstart])*(1.0 + tra[rstart1+rstart]))/avefit1;

////// plasticity = sign epistasis
/*
epis = 0.1234;

        avefit1 = ab1*(1.0 -2.0*tra[rstart1+rstart]) + Ab1*(1.0 + 2.0*tra[rstart1+rstart]) + aB1*(1.0 -2.0*tra[rstart1+rstart]*(1-pla)) +  AB1*(1.0 +2.0*tra[rstart1+rstart]*(1-pla));

        fit_ab1=(1.0 - 2.0*tra[rstart1+rstart])/avefit1;
        fit_Ab1= (1.0 + 2.0*tra[rstart1+rstart])/avefit1;
        fit_aB1= (1.0 - 2.0*tra[rstart1+rstart]*(1-pla))/avefit1;
*/
////////common to models:

        if(ab1>0.0) ab1 = fit_ab1*ab1;
        if(aB1>0.0) aB1 = fit_aB1*aB1;
        if(Ab1>0.0) Ab1 = fit_Ab1*Ab1;
        if(AB1>0.0) AB1 = 1.0 - ab1 - aB1 - Ab1;

                     ////////////////    population 2    --- drift (C = 0) under all scenarios

        rstart1+=1;

   ////////////// RECCURENT MUATION  ////////////////

                        rab1 = (aB1 + Ab1)*NMU + (ab1 * (1 - 2*NMU)); 
                        raB1 = (ab1 + AB1)*NMU + (aB1 * (1 - 2*NMU)); 
                        rAb1 = (AB1 + ab1)*NMU + (Ab1 * (1 - 2*NMU)); 
                        rAB1 = (Ab1 + aB1)*NMU + (AB1 * (1 - 2*NMU)); 

                        rab2 = (aB2 + Ab2)*NMU + (ab2 * (1 - 2*NMU));
                        raB2 = (ab2 + AB2)*NMU + (aB2 * (1 - 2*NMU));
                        rAb2 = (AB2 + ab2)*NMU + (Ab2 * (1 - 2*NMU));
                        rAB2 = (Ab2 + aB2)*NMU + (AB2 * (1 - 2*NMU));


                        ab1 = rab1;
                        aB1 = raB1;
                        Ab1 = rAb1;
                        AB1 = rAB1;

                        ab2 = rab2;
                        aB2 = raB2;
                        Ab2 = rAb2;
                        AB2 = rAB2;


   /////////////////  RECOMBINATION   ///////////////

	/////////////  population 1  ///////////////

	if((ab1+aB1)*(Ab1+AB1)*(ab1+Ab1)*(aB1+AB1) > 0.0)
	{

	  D1 = ((ab1*AB1)-(aB1*Ab1))*REC;

	  ab1 = ab1 - D1;
	  aB1 = aB1 + D1;
	  Ab1 = Ab1 + D1;
	  AB1 = AB1 - D1;

	}

	/////////////  population 2  ///////////////

	if((ab2+aB2)*(Ab2+AB2)*(ab2+Ab2)*(aB2+AB2) > 0.0)
	{

	  D2 = ((ab2*AB2)-(aB2*Ab2))*REC;

	  ab2 = ab2 - D2;
	  aB2 = aB2 + D2;
	  Ab2 = Ab2 + D2;
	  AB2 = AB2 - D2;

	}

	/////////////  REPRODUCTION  ///////////////

    //////////// population 1  //////////////

    if(Ab1 > 0.0) r1 = Ab1/(Ab1+AB1);
    else r1 = 0.0;

    if(aB1 > 0.0) r2 = aB1/(Ab1+aB1+AB1);
    else r2 = 0.0;

    ab1 = bnldev(ab1, N1, &seed);
    R1 = N1 - (long)ab1;

    aB1 = bnldev(r2, R1, &seed);
    R2 = R1 - (long)aB1;

    Ab1 = bnldev(r1, R2, &seed);

    AB1 = ((double)N1 - ab1 - aB1 - Ab1)/(double)N1;
    ab1 = ab1/(double)N1;
    aB1 = aB1/(double)N1;
    Ab1 = Ab1/(double)N1;

    //////////// population 2  //////////////

    if(Ab2 > 0.0) s1 = Ab2/(Ab2+AB2);
    else s1 = 0.0;

    if(aB2 > 0.0) s2 = aB2/(Ab2+aB2+AB2);
    else s2 = 0.0;

    ab2 = bnldev(ab2, N2, &seed );
    S1 = N2 - (long)ab2;

    aB2 = bnldev(s2, S1, &seed );
    S2 = S1 - (long)aB2;

    Ab2 = bnldev(s1, S2, &seed );

    AB2 = ((double)N2 - ab2 - aB2 - Ab2)/(double)N2;
    ab2 = ab2/(double)N2;
    aB2 = aB2/(double)N2;
    Ab2 = Ab2/(double)N2;

    //////////////// SUMMARIES ////////////////////

	hl_first_locus  += 2*( (ab1+aB1)*(1.0-(ab1+aB1))*pop1+(ab2+aB2)*(1.0-(ab2+aB2))*(1.0-pop1));
	hl_second_locus += 2*( (ab1+Ab1)*(1.0-(ab1+Ab1))*pop1+(ab2+Ab2)*(1.0-(ab2+Ab2))*(1.0-pop1));

	ax = (((ab1+aB1)*(double)N1 + (ab2+aB2)*(double)N2))/((double)(N1+N2));
	bx = (((ab1+Ab1)*(double)N1 + (ab2+Ab2)*(double)N2))/((double)(N1+N2));

if(ax*bx*(1.0-ax)*(1.0-bx) > 0.000000001) cosegtime +=1;

  }} /// end time and reps

fprintf(out1, "epis %.4f Nmu %.3f P %ld smax %.3f N1 %ld N2 %ld MIGRANTS %.1f hl_first_locus %.4f hl_second_locus %.4f totalcosegtime %lld cosegregation_time %.3Lf TIME %ld REPS %ld\n", epis, NMU, PHASE, smax, N1, N2, MIGS, hl_first_locus/(double)REP ,hl_second_locus/(double)REP, cosegtime, (long double)cosegtime/(long double)REP, TIME, REP);

fclose(par);
fclose(inp);
fclose(out1);

return 0;

} // END MAIN //

double randn (double mu, double sigma){
  double U1, U2, W, mult;
  static double X1, X2;
  static int call = 0;

  if (call == 1)
    {
      call = !call;
      return (mu + sigma * (double) X2);
    }

  do
    {
      U1 = -1 + ((double) rand () / RAND_MAX) * 2;
      U2 = -1 + ((double) rand () / RAND_MAX) * 2;
      W = pow (U1, 2) + pow (U2, 2);
    }
  while (W >= 1 || W == 0);

  mult = sqrt ((-2 * log (W)) / W);
  X1 = U1 * mult;
  X2 = U2 * mult;

  call = !call;

  return (mu + sigma * (double) X1);
}


double bnldev(double pp, long n, long *idum)
{

long j;

static long nold =(-1);

double am, em, g, angle, p, bnl, sq, t, y;
static double pold=(-1.0), pc, plog, pclog, en, oldg;

p=(pp <= 0.5 ? pp : 1.0-pp);

am = n*p;

if(n<25){
	bnl = 0.0;

	for(j=1; j<=n; j++)
		if(ran1(idum) < p) ++bnl;
}else if (am < 1.0) {
  	g=exp(-am);
  	t = 1.0;
	for (j=0; j<=n; j++){
		t *= ran1(idum);
		if (t < g) break;
}
bnl=(j<=n ? j:n);
}else {


if (n != nold) {
	en = n;
	oldg = gammln(en+1.0);
	nold = n;
}

if (p != pold) {
	pc = 1.0-p;
	plog = log(p);
	pclog = log(pc);
	pold = p;
}

sq = sqrt(2.0*am*pc);
do{
	do{
		angle = PI * ran1(idum);
		y = tan(angle);
		em = sq*y+am;
	} while (em < 0.0 || em >= (en + 1.0));

	em = floor(em);
	t = 1.2 * sq*(1.0+y*y)*exp(oldg-gammln(em+1.0)-gammln(en-em+1.0)+em*plog+(en-em)*pclog);
}while (ran1(idum) > t);


bnl = em;

}

if (p != pp) bnl = n - bnl;
return bnl;

}

double gammln(double xx)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	long j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}

double ran1(long *idum)
{
	long j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	double temp;

	if (*idum <= 0 || !iy) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k;
			if (*idum < 0) *idum += IM;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if (*idum < 0) *idum += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = *idum;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}


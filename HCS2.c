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

long N1, N2, TIME, REP, PHASE; double NMU, MIGS, sel, REC, pla, epis, C;

FILE * par, * inp, * out1, * out2; 
par = fopen("pars1", "r"); 
inp = fopen("trajectory", "r");
out1 = fopen("outst1", "w");
out2 = fopen("outst2", "w");

fscanf(par, "%ld", &PHASE);
fscanf(par, "%lf", &MIGS);

C = 1, N1 = 5000, N2 = 5000, TIME = 100*(N1+N2), REP = 10*TIME; REC = 0.5; NMU = 0.1;

epis = 0.1;     ////// = 0.0 additivitiy ,   < 0.0 diminishing (negative)   , > 0.0 reinforcing (positive)    ///// can define epistasis here, except sign epi with is defined by commenting and uncomenting code in selection section
pla = epis;

long seed = -2345;

double D1, D2, r1, r2, s1, s2;  int one, two;  long R1, R2, S1, S2; long double Ds1 = 0.0, Ds2 = 0.0, dt1 = 0.0, dt2 = 0.0, DA1 = 0.0, DA2 = 0.0;
int mode, SE;
long start_time_first = 0, start_time_second = 0;
long fixation_time_first = 0, fixation_time_second = 0, number_of_fixed_first = 0, number_of_fixed_second = 0, segtime_first = 0, segtime_second = 0;
long loss_time_first = 0, loss_time_second = 0, number_of_loss_first = 0, number_of_loss_second = 0;
double avg_fix_time_first = 0.0, avg_fix_time_second = 0.0, avg_loss_time_first = 0.0, avg_loss_time_second = 0.0, prob_fix_first = 0.0, prob_fix_second = 0.0, ave_segtime_first = 0.0, ave_segtime_second = 0.0;
double avefit1, avefit2, select = 0.0 ;
double ab1, aB1, Ab1, AB1, ab2, aB2, Ab2, AB2, fit_ab1, fit_aB1, fit_Ab1, fit_ab2, fit_aB2, fit_Ab2;   
double newN1, newN2, ranpick1, ranpick2, m,m1,m2, ab1_migs=0.0, ab1_stay, aB1_migs=0.0, aB1_stay, Ab1_migs=0.0, Ab1_stay, AB1_migs=0.0, AB1_stay, ab2_migs=0.0, ab2_stay, aB2_migs=0.0, aB2_stay, Ab2_migs=0.0, Ab2_stay, AB2_migs=0.0, AB2_stay, ab1new, Ab1new, aB1new, AB1new, ab2new, Ab2new, aB2new, AB2new;
double hl_first_locus = 0.0, hl_second_locus = 0.0; long seg100A = 0, seg4A = 0, seg40A = 0, seg10A = 0, seg100B = 0, seg4B = 0, seg40B = 0, seg10B = 0,  seg100AB = 0, seg4AB = 0, seg40AB = 0, seg10AB = 0;
long cc, k, i, l;  double a,b, ax, bx, tra[3*PHASE], smax; int rstart, rstart1, xx; long long gentime = 0, cosegtime = 0;
double pop1 = (double)(N1)/(double)(N1+N2);

SE = 1;

	for (l = 0; l < (3*PHASE); l++)
		{
		fscanf(inp, "%lf", &tra[l]);
		}

 smax = tra[PHASE/4];

for (i = 0 ; i < REP; i++)
{
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
	    	one = 1;
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
	    	two = 1;
	  }
	  start_time_first = 0;
	  start_time_second = 0;

	  rstart = (int)(ran1(&seed)*PHASE);
	  rstart1 = 0; 

  for( k = 0 ; k< TIME ; k++)
  {

////////////////  MUTATION ////////////////////

/////////////  first locus  /////////////

   if (ran1(&seed) < NMU && one == 0)
     {

       if(ran1(&seed) < pop1)
         {
             if(ran1(&seed) < ab1)
             {
               ab1 = ab1 - 1.00/(double)N1;
               Ab1 = Ab1 + 1.00/(double)N1;
             }
             else
             {
               aB1 = aB1 - 1.00/(double)N1;
               AB1 = AB1 + 1.00/(double)N1;
             }

         }
       else
         {
             if(ran1(&seed) < ab2)
             {
               ab2 = ab2 - 1.00/(double)N2;
               Ab2 = Ab2 + 1.00/(double)N2;
             }

             else
             {
               aB2 = aB2 - 1.00/(double)N2;
               AB2 = AB2 + 1.00/(double)N2;
             }
         }
       one = 1;
       start_time_first = k;
       }

///////////////  second locus  /////////////

			   if (ran1(&seed) < NMU && two == 0)
			      {

			        if (ran1(&seed) < pop1)
			           {
			             if(ran1(&seed) < ab1)
			             {
			               ab1 = ab1 - 1.00/(double)N1;
			               aB1 = aB1 + 1.00/(double)N1;
			             }
			             else
			             {
			               Ab1 = Ab1 - 1.00/(double)N1;
			               AB1 = AB1 + 1.00/(double)N1;
			             }
			           }
			        else
			           {
			             if(ran1(&seed) < ab2)
			             {
			               ab2 = ab2 - 1.00/(double)N2;
			               aB2 = aB2 + 1.00/(double)N2;
			             }
			             else
			             {
			               Ab2 = Ab2 - 1.00/(double)N2;
			               AB2 = AB2 + 1.00/(double)N2;
			             }
			          }
          two = 1;
	  start_time_second = k;
        }

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

          newN1 = newN1 - 1.0;
          newN2 = newN2 - 1.0;

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
  sel = tra[rstart1+rstart];

	////////////////////////////////////////   POPULATION 1 ///////////////////////////////////////////

/////////////// reinforcing and diminishing

mode = 1;
      avefit1 = ab1*(1.0 - sel*(1.0+epis))*(1.0 -sel*(1.0+epis)) + Ab1*(1.0 + sel)*(1.0 - sel) + aB1*(1.0 - sel)*(1.0 + sel) +  AB1*(1.0 + sel*(1.0+epis))*(1.0 + sel*(1.0+epis));

        fit_ab1= ((1.0 - sel*(1.0+epis))*(1.0 -sel*(1.0+epis)))/avefit1;
        fit_Ab1= ((1.0 + sel)*(1.0 - sel))/avefit1;
        fit_aB1= ((1.0 - sel)*(1.0 + sel))/avefit1;

///////////// plasticity 
/*
mode = 2;
        avefit1 = ab1*(1.0 - sel) + Ab1*(1.0 + sel) + aB1*(1.0 -sel*(1-pla)) +  AB1*(1.0 + sel*(1-pla));

        fit_ab1=(1.0 - sel)/avefit1;
        fit_Ab1= (1.0 + sel)/avefit1;
        fit_aB1= (1.0 - sel*(1-pla))/avefit1;
*/
////////common to models:

        if(ab1>0.0) ab1 = fit_ab1*ab1;
        if(aB1>0.0) aB1 = fit_aB1*aB1;
        if(Ab1>0.0) Ab1 = fit_Ab1*Ab1;
        if(AB1>0.0) AB1 = 1.0 - ab1 - aB1 - Ab1;

///////////////////////////////    POPULATION 2    --- uncomment the same portion of the code including common part as above for uniform (no storage effect) selection //////////////////////////////////////////////////////////////////

      sel = C * sel;

      avefit2 = ab2*(1.0 - sel*(1.0+epis))*(1.0 -sel*(1.0+epis)) + Ab2*(1.0 + sel)*(1.0 - sel) + aB2*(1.0 - sel)*(1.0 + sel) +  AB2*(1.0 + sel*(1.0+epis))*(1.0 + sel*(1.0+epis));

        fit_ab2= ((1.0 - sel*(1.0+epis))*(1.0 -sel*(1.0+epis)))/avefit2;
        fit_Ab2= ((1.0 + sel)*(1.0 - sel))/avefit2;
        fit_aB2= ((1.0 - sel)*(1.0 + sel))/avefit2;

///////////// plasticity
/*
      sel = C * sel;

        avefit2 = ab2*(1.0 - sel) + Ab2*(1.0 + sel) + aB2*(1.0 -sel*(1-pla)) +  AB2*(1.0 + sel*(1-pla));

        fit_ab2=(1.0 - sel)/avefit2;
        fit_Ab2= (1.0 + sel)/avefit2;
        fit_aB2= (1.0 - sel*(1-pla))/avefit2;
*/
////////common to models in pop 2:

SE = 0;
        if(ab2>0.0) ab2 = fit_ab2*ab2;
        if(aB2>0.0) aB2 = fit_aB2*aB2;
        if(Ab2>0.0) Ab2 = fit_Ab2*Ab2;
        if(AB2>0.0) AB2 = 1.0 - ab2 - aB2 - Ab2;

        //////////////////////  end pop2

        rstart1+=1;

///////////////////////////////////////////////////////////  RECOMBINATION   /////////////////////////////////////////////////////////////////////

/////////////  population 1  ///////////////

	if((ab1+aB1)*(Ab1+AB1)*(ab1+Ab1)*(aB1+AB1) > 0.0)
	{

	  D1 = ((ab1*AB1)-(aB1*Ab1))*REC;

	  ab1 = ab1 - D1;
	  aB1 = aB1 + D1;
	  Ab1 = Ab1 + D1;
	  AB1 = AB1 - D1;

	if (D1 < 0){
          	if (-1.0*(Ab1+AB1)*(aB1+AB1) > -1.0*(ab1+aB1)*(ab1+Ab1)) Ds1 += D1/(REC*-1.0*(Ab1+AB1)*(aB1+AB1));
          	else Ds1 += D1/(REC*-1.0*(ab1+aB1)*(ab1+Ab1));
        } 
	else if (D1 > 0){
          	if ((ab1+aB1)*(aB1+AB1) < (Ab1+AB1)*(ab1+Ab1)) Ds1 += D1/(REC*(ab1+aB1)*(aB1+AB1));
          	else Ds1 += D1/(REC*(Ab1+AB1)*(ab1+Ab1));
        }
	DA1 += D1/REC;
	dt1 += 1.0;
	
	}

	/////////////  population 2  ///////////////

	if((ab2+aB2)*(Ab2+AB2)*(ab2+Ab2)*(aB2+AB2) > 0.0)
	{

	  D2 = ((ab2*AB2)-(aB2*Ab2))*REC;

	  ab2 = ab2 - D2;
	  aB2 = aB2 + D2;
	  Ab2 = Ab2 + D2;
	  AB2 = AB2 - D2;

      	if (D2 < 0){
                if (-1.0*(Ab2+AB2)*(aB2+AB2) > -1.0*(ab2+aB2)*(ab2+Ab2)) Ds2 += D2/(REC*-1.0*(Ab2+AB2)*(aB2+AB2));
                else Ds2 += D2/(REC*-1.0*(ab2+aB2)*(ab2+Ab2));
        }
        else if (D2 > 0){
                if ((ab2+aB2)*(aB2+AB2) < (Ab2+AB2)*(ab2+Ab2)) Ds2 += D2/(REC*(ab2+aB2)*(aB2+AB2));
                else Ds2 += D2/(REC*(Ab2+AB2)*(ab2+Ab2));
        }
        DA2 += D2/REC;
       	dt2 += 1.0;

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

	if((ax*(1-ax) == 0.0) && one == 1)
	  {
	  if(ax == 0.0)
	    {
	      fixation_time_first += k-start_time_first;
	      number_of_fixed_first +=1;
	    }

	  else if(ax == 1.0)
	    {
	      loss_time_first += k-start_time_first;
	      number_of_loss_first +=1;
	    }
	  segtime_first += k-start_time_first;
	  one = 2;
	  }
	    if((bx*(1-bx) == 0.0) && two == 1)
	  {

	  if(bx == 0.0)
	    {
	      fixation_time_second += k-start_time_second;
	      number_of_fixed_second +=1;
	    }

	  else if(bx == 1.0)
	    {
	      loss_time_second += k-start_time_second;
	      number_of_loss_second +=1;
	    }
	  segtime_second +=  k-start_time_second;
	  two = 2;
	  }

	  if(one > 1 && two > 1)
	  {
	  break;
	  }

if(k == (4*(N1+N2))){
if(one == 1) seg4A +=1;
if(two == 1) seg4B +=1;
if(one == 1 && two == 1) seg4AB +=1;
}

if(k == (40*(N1+N2))){
if(one == 1) seg40A +=1;
if(two == 1) seg40B +=1;
if(one == 1 && two == 1) seg40AB +=1;
}

if(k == (10*(N1+N2))){
if(one == 1) seg10A +=1;
if(two == 1) seg10B +=1;
if(one == 1 && two == 1) seg10AB +=1;
}

if(one == 1 && two == 1) cosegtime +=1;

  } // end TIME

	gentime += k;

if(k == TIME){	
	if(one == 1) seg100A +=1;
	if(two == 1) seg100B +=1;
	if(one == 1 && two ==1) seg100AB +=1;
}

	avg_fix_time_first = (double)fixation_time_first/(double)number_of_fixed_first;
	avg_fix_time_second = (double)fixation_time_second/(double)number_of_fixed_second;

	avg_loss_time_first = (double)loss_time_first/(double)number_of_loss_first;
	avg_loss_time_second = (double)loss_time_second/(double)number_of_loss_second;

	prob_fix_first = (double)number_of_fixed_first/(double)(REP);
	prob_fix_second = (double)number_of_fixed_second/(double)(REP);

	ave_segtime_first = (double)segtime_first/(double)(REP);
	ave_segtime_second = (double)segtime_second/(double)(REP);

if( i == 299999 && hl_first_locus > 10000000.0 && hl_second_locus > 10000000.0){REP= 500000; break;} 

} // end REP

fprintf(out1, "mode %d SE %d epis %.4f C %.4f Nmu %.3f P %ld smax %.3f N1 %ld N2 %ld MIGRANTS %.1f hl_ave_locus %.4f hl_drift %.4f avg_fix_time %.4f avg_loss_time %.1f REPS %ld\n",mode, SE, epis, C, NMU, PHASE, smax, N1, N2, MIGS, (hl_first_locus+hl_second_locus)/(2.0*(double)REP), ((hl_first_locus+hl_second_locus)/(2.0*(double)REP))/2, (double)(fixation_time_first+fixation_time_second)/(double)(number_of_fixed_second+number_of_fixed_first), (double)(loss_time_first+loss_time_second)/(double)(number_of_loss_second+number_of_loss_first), REP);
fprintf(out2, "mode %d SE %d epis %.4f Nmu %.3f P %ld smax(assuming sine fn) %.3f N1 %ld N2 %ld MIGRANTS %.1f numberseg100Nfirst %ld  numberseg4Nfirst %ld  numberseg40Nfirst %ld  numberseg10Nfirst %ld numberseg100Nsecond %ld  numberseg4Nsecond %ld  numberseg40Nsecond %ld  numberseg10Nsecond %ld ", mode, SE, epis, NMU, PHASE, smax, N1, N2, MIGS, seg100A, seg4A, seg40A, seg10A, seg100B, seg4B, seg40B, seg10B);
fprintf(out2, "numberseg100Nboth %ld  numberseg4Nboth %ld  numberseg40Nboth %ld  numberseg10Nboth %ld hl_first_locus  %.4f hl_second_locus   %.4f avg_fix_time_first   %.4f avg_fix_time_second  %.4f avg_loss_time_first   %.4f avg_loss_time_second  %.4f prob_fix_first  %.8f prob_fix_second  %.8f ", seg100AB, seg4AB, seg40AB, seg10AB, hl_first_locus/(double)REP ,hl_second_locus/(double)REP, avg_fix_time_first, avg_fix_time_second, avg_loss_time_first, avg_loss_time_second, prob_fix_first, prob_fix_second);
fprintf(out2, "D1standard %.6Lf D1ave %.6Lf D2standard %.6Lf D2ave %.6Lf dt1 %Lf  dt2 %Lf", Ds1, DA1, Ds2, DA2, dt1, dt2);
fprintf(out2, "genomesegtime %.3Lf totalcosegtime %lld cosegregation_time %.3Lf TIME %ld REPS %ld\n", (long double)gentime/(long double)REP, cosegtime, (long double)cosegtime/(long double)REP, TIME, REP);

fclose(par);
fclose(inp);
fclose(out1);
fclose(out2);

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
} else {

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
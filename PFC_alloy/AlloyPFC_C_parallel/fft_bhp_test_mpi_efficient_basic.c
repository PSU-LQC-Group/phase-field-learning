/*Binary Alloy Phase Field Crystal modelling code using 
semi implicit Fourier space PFC algorithm v2
Written by Jonathan Stolle Oct. 2010 based on Ken Elder's EXD.f

-periodic boundary conditions
-initialize states:
ntype:
2 - restart from previous data
5 - basic initial condition (solid seed in the middle of a liquid pool)

corrections:
-to be more efficient, a fourier transform was neglected
-solves equation (9.50) in the textbook, assuming lattice parameter is
constant; that is, it excludes effects due to Vegard's Law (eta = 0) 
*/

/* libraries -- fftw libraries might have 
different names e.g., dfftw, sfftw */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <fftw.h>
#include <rfftw.h>
#include <rfftw_mpi.h>
#include <fftw_mpi.h>

/*Constant Declarations*/
#define PI 3.14159265 /* Approximation */
#define DEG_TO_RAD PI/180 /* Convert Degrees to Radians */
#define ranseed 300 /*set this manually for now*/

/*global variables -- those which change rarely*/
long int Nx,Ny,N2,Ny2; /*needed for in place transform!*/ 
	/*size of y-dimension in Fourier space in r2c transform*/
	
/*other variables introduced because of MPI*/
int local_nlast, local_last_start, local_nlast2_after_trans,
	local_last2_start_after_trans, total_local_size;
int comm3d,ierr,irc,stride,USE_WORK;
int myid,numprocs;
int curr_proc;
	
/*input parameters*/
char densf[5], concf[5]; /*make sure length > length string read in */
char stemp[20]; /*make sure length > length string read in */

double dx, dt;	//grid spacing/timestep
double theta;	//angle
double BL, BL2, BX, RK, w, t, u, v; //PFC parameters
double rno, co, cl;	//average initial field values
long int nnend, nout, nstart; //number of iterations/output files
int ntype;	//initialization condition type
double qo;	//magnitude of reciprocal lattice vector
double dLy, dLx; //determine seed size of seeds
double noise; //determine amplitude of noise

/*variables first defined in the fourier factors part*/
double facx, facy, rli, enfac,qfac, qcon,fac,ef;
double *qsf;	//operations for Vegard's law term
double *ql,*qlc,*qnc,*qn; //operations on linear/nonlinear parts

/*function prototypes */
int output_to_file(double *x, char *name, long int ns);
int initialize_fields(double *rn, double *c);
int calculate_factors();

/*main function*/
int main(int argc, char *argv[]){
	/*important overall variables*/
	//field variables
	double *rn,*c; /*density and concentration*/
	
	/*fourier transform variables*/
	fftw_complex *data,*datac; /*density and concentration data*/
	fftw_complex *cc,*crn;	//complex pointer to field variable arrays
	/*real space equivalents of above variables*/
	fftw_real *tdata,*tdatac;	
	fftw_real *work;	//temporary storage variable
	
	/*various counter variables*/
	long int ix, iy, ic, ns,kk,ixh; 

	/*input and output files*/	
	FILE *in,*out;

	/*MPI_variables -- use if variables passed through window
	do not work*/
	/*int argc;
	char **argv;*/

	/*mpi transform*/
	rfftwnd_mpi_plan planf,plani;

	/*start mpi */
	MPI_Init(&argc,&argv); 	
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD,&numprocs);

	printf("\nmyid: %d, numprocs: %d",myid,numprocs);
	
	/*read in parameter file*/ 
	in=fopen("alloy.in","rt");

	fscanf(in,"%s %s",densf, stemp);
	fscanf(in,"%s %s",concf, stemp);
	fscanf(in,"%ld %s %ld %s",&Nx,stemp,&Ny,stemp);
	fscanf(in,"%lf %s %lf %s",&dx,stemp,&dt,stemp);
	fscanf(in,"%lf%s %lf%s %lf%s",&BL,stemp,&BL2,stemp,&BX,stemp);
	fscanf(in,"%lf%s %lf%s %lf%s",&RK,stemp,&w,stemp,&t,stemp);
	fscanf(in,"%lf%s %lf%s",&u,stemp,&v,stemp);
	fscanf(in,"%lf%s %lf%s %lf%s",&rno,stemp,&co,stemp,&cl,stemp);
	fscanf(in,"%ld%s %ld%s %ld%s",&nnend,stemp,&nout,stemp,&nstart,stemp);
	fscanf(in,"%d%s %lf%s",&ntype,stemp,&qo,stemp);
	fscanf(in,"%lf%s%lf%s%lf%s",&dLx,stemp,&dLy,stemp,&noise,stemp);
	fclose(in);
	
	/*size of y-dimension for real fftw variables*/
	N2=Ny/2+1,Ny2=N2*2; 
	
	theta=theta*DEG_TO_RAD; //convert to radians

	/*set up parallel transform -- create plan*/
	planf=rfftw2d_mpi_create_plan(MPI_COMM_WORLD,Nx, Ny,
	FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE); //add FFTW_MEASURE command
	plani=rfftw2d_mpi_create_plan(MPI_COMM_WORLD,Nx, Ny,
	FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE);
	//determine how memory is allocated on each processor
	rfftwnd_mpi_local_sizes(planf, &local_nlast, &local_last_start,
							&local_nlast2_after_trans,
							&local_last2_start_after_trans,
							&total_local_size);
	//display memory allocation results
	MPI_Barrier(MPI_COMM_WORLD);   
	printf("\nmyid: %d, local_last_start: %d, local_nlast2_after_trans: %d"
			,myid,local_last_start,
			local_last2_start_after_trans);
	//MPI_Barrier(MPI_COMM_WORLD);   //space evenly
	printf("\n myid: %d, local_nlast: %d",myid,local_nlast);
	//MPI_Barrier(MPI_COMM_WORLD);   //space evenly
	printf("\nmyid: %d, local_last_after: %d, total_local_size: %d"
			,myid,local_last2_start_after_trans,total_local_size);
  
	
	/*assign space to make arrays -- total_local_size excess for rn,c?*/
	rn = (double*)malloc ( sizeof ( double ) * total_local_size );
	c = (double*)malloc ( sizeof ( double ) * total_local_size );
	ql=(double*)malloc ( sizeof ( double ) * total_local_size );
	qlc=(double*)malloc ( sizeof ( double ) * total_local_size );
	qnc=(double*)malloc ( sizeof ( double ) * total_local_size );
	qn=(double*)malloc ( sizeof ( double ) * total_local_size );
	tdata = (double*)malloc ( sizeof ( double ) * total_local_size );
	tdatac = (double*)malloc ( sizeof ( double ) * total_local_size );
	work = (double*)malloc ( sizeof ( double ) * total_local_size );

	/* complex data pointing to same place as real data*/
	data = (fftw_complex*) tdata;
	datac = (fftw_complex*) tdatac;
	cc = (fftw_complex*) c;
	crn = (fftw_complex*) rn;

	/*initialize data*/
	initialize_fields(rn,c);

	/*set Fourier factors*/
	calculate_factors();	

	/*run main program*/
	for (ns=nstart+1;ns<(nnend+1);ns++){
	    
		/*calculate nonlinear parts in real space*/
		for (ix=0;ix<local_nlast;ix++)
			for (iy=0;iy<Ny2;iy++){
				tdatac[ix*Ny2+iy]=c[ix*Ny2+iy]*(rn[ix*Ny2+iy]*
				rn[ix*Ny2+iy]*BL2+u*c[ix*Ny2+iy]*c[ix*Ny2+iy])+
				ef*rn[ix*Ny2+iy];
				tdata[ix*Ny2+iy]=rn[ix*Ny2+iy]*(BL2*c[ix*Ny2+iy]
				*c[ix*Ny2+iy]+rn[ix*Ny2+iy]*(-t+v*rn[ix*Ny2+iy]))
				+ef*(c[ix*Ny2+iy]);
			}

		/*transform fields into Fourier space*/
		MPI_Barrier(MPI_COMM_WORLD);
		rfftwnd_mpi(planf, 1, rn, work, FFTW_NORMAL_ORDER);
		MPI_Barrier(MPI_COMM_WORLD);
		rfftwnd_mpi(planf, 1, c, work, FFTW_NORMAL_ORDER);
		MPI_Barrier(MPI_COMM_WORLD);
		/*transform nonlinear parts into Fourier space*/
		rfftwnd_mpi(planf, 1, tdata, work, FFTW_NORMAL_ORDER);
		MPI_Barrier(MPI_COMM_WORLD);
		rfftwnd_mpi(planf, 1, tdatac, work, FFTW_NORMAL_ORDER);

		/*fourier operations to increment density 
		linear part x linear multiplier
		nonlinear part x nonlinear multiplier */
		for (ix=0;ix<local_nlast;ix++)
			for (iy=0;iy<N2;iy++){
				crn[ix*N2+iy].re=crn[ix*N2+iy].re*ql[ix*Ny+iy]+
				(data[ix*N2+iy].re)*qn[ix*Ny+iy];
				crn[ix*N2+iy].im=crn[ix*N2+iy].im*ql[ix*Ny+iy]+
				(data[ix*N2+iy].im)*qn[ix*Ny+iy];
			}

		/*fourier operations to increment concentration */
		for (ix=0;ix<local_nlast;ix++)
		for (iy=0;iy<N2;iy++){
			cc[ix*N2+iy].re=cc[ix*N2+iy].re*qlc[ix*Ny+iy]+
			datac[ix*N2+iy].re*qnc[ix*Ny+iy];
			cc[ix*N2+iy].im=cc[ix*N2+iy].im*qlc[ix*Ny+iy]+
			datac[ix*N2+iy].im*qnc[ix*Ny+iy];
			}
		
		/*inverse Fourier transform results to update timestep*/
		MPI_Barrier(MPI_COMM_WORLD);
		rfftwnd_mpi(plani, 1, rn, work, FFTW_NORMAL_ORDER);
		MPI_Barrier(MPI_COMM_WORLD);
		rfftwnd_mpi(plani, 1, c, work, FFTW_NORMAL_ORDER);

		if (ns%nout==0){
			MPI_Barrier(MPI_COMM_WORLD);
			/*ouput data to file at each timestep*/
			output_to_file(rn, densf, ns);
			output_to_file(c, concf, ns);
		}
	} //end of main loop
	
	/*free memory allocated to arrays*/
	free(rn);	free(c);	free(qn);	free(qnc);
	free(ql);	free(qlc);	free(tdata);
	free(tdatac);	free(work);
	

	/*free memory of plans*/
	rfftwnd_mpi_destroy_plan(planf);
	rfftwnd_mpi_destroy_plan(plani);
	
	/*end program*/
	MPI_Finalize();
	
	return(0);
}

/*functions*/
int initialize_fields(double *rn, double *c){
	/*initialize variables*/
	
	/*variables first defined in the fields initialization*/
	double ax, ay1, ay2,Amp,rrr;
	/*various counter variables*/
	long int ix, iy, ic, ns,kk,ixh; 
	double xx,yy;

	/*file variables*/
	FILE *in, *out;
	
	/*string variables and miscellaneous*/
	char s[9],temp1[20],temp2[20];
	double temp;

	/*start from a previous run*/
	if(ntype==2){
		
		/*call init_random_seed(iwarm)*/
		srand(time(NULL)); /* initialize random number generator*/

		/*density name file*/
		kk=sprintf(s,"%ld",nstart);
		strcpy(temp1,s);
		strcpy(temp2,densf);
		
		//open density file
		in=fopen(strcat(temp2,strcat(temp1,".out")),"rt");        
	  
		/*read in density data*/
		for (ix=0;ix<local_last_start;ix++)
			for (iy=0;iy<Ny;iy++)
				fscanf(in,"%lf",&temp);
				
		for(ix=0;ix<local_nlast;ix++)
			for (iy=0;iy<Ny;iy++)
				fscanf(in,"%lf",&(rn[ix*Ny2+iy]));
		
		/*close file*/
		fclose(in);
		
		/*concentration name file*/
		kk=sprintf(s,"%ld",nstart);
		strcpy(temp1,s);
		strcpy(temp2,concf);

		//open concentration file
		in=fopen(strcat(temp2,strcat(temp1,".out")),"rt");        
		
		/*read in concentration data*/
		for (ix=0;ix<local_last_start;ix++)
			for (iy=0;iy<Ny;iy++)
				fscanf(in,"%lf",&temp);
				
		for(ix=0;ix<local_nlast;ix++)
			for (iy=0;iy<Ny;iy++)
				fscanf(in,"%lf",&(c[ix*Ny2+iy]));
		
		/*close file*/
		fclose(in);		  
	}

	/*initialize rn and c fields -- single solid seed in liquid*/
	if(ntype==5){
		srand(time(NULL)); /* initialize random number generator*/
		ns=nstart;	//initial run

		/*initialize variables used for ic's*/        
		Amp=0.5;
		ax=qo*dx;
		ay1=qo*dx/sqrt(3.0);
		ay2=2.0*ay1;
		theta=0;

		for(ix=0;ix<local_nlast;ix++){
			ixh=ix+local_last_start; //x position depends on node
			for (iy=0;iy<Ny;iy++){
				/*solid*/
				if((ixh>((Nx-dLx)*0.5))&&(ixh<((Nx+dLx)*0.5))&&
				(iy>((Ny-dLy)*0.5))&&(iy<((Ny+dLy)*0.5))){
					//transform coordinates to crystal orientation
					xx=ixh*cos(theta/2)-iy*sin(theta/2);
					yy=ixh*sin(theta/2)+iy*cos(theta/2);
					//oscillatory density field (1 mode approximation)
					rn[ix*Ny2+iy]=Amp*(0.5*cos(ay2*yy)-cos(ax*xx)*
					cos(ay1*yy))+rno;
					//noise in concentration field
					rrr=((double)rand()/RAND_MAX-0.5)*noise;
					//concentration field in solid
					c[ix*Ny2+iy]=co+rrr;
				}
				/*liquid*/
				else {
					//constant density field
					rn[ix*Ny2+iy]=rno;
					//noise in concentration field
					rrr=((double)rand()/RAND_MAX-0.5)*noise;
					//concentration field in liquid
					c[ix*Ny2+iy]=cl+rrr;
				}
           }
		}
  		MPI_Barrier(MPI_COMM_WORLD);
		
		/*output density data*/
		output_to_file(rn,densf,ns);

		/*output concentration data*/
		output_to_file(c,concf,ns);
		
		}
	return(0);
}

int output_to_file(double *x, char *name, long int ns){
/*output arrange to a file*/

	/*declare variables*/
	FILE *out;
	/*variables for local operations*/
	char s[9],temp1[20],temp2[20];
	int ix, iy,kk;

	/*get file name*/
	kk=sprintf(s,"%ld",ns);
	strcpy(temp1,s);
	strcpy(temp2,name);

	/*output data to file*/
	
	for (kk=0;kk<numprocs;kk++){
		if (myid==kk){
			/*open file for output*/
			if (myid==0)
				out=fopen(strcat(temp2,strcat(temp1,".out")),"wt");
			else
				out=fopen(strcat(temp2,strcat(temp1,".out")),"at");
				
			/*output data one processor at a time*/	
			for (ix=0;ix<local_nlast;ix++) for (iy=0;iy<Ny;iy++)
				fprintf(out,"%.10lf \n",x[ix*Ny2+iy]);
				
			/*close file*/
			fclose(out);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	
	return(0);
}

int calculate_factors(){
	int ix,ixh,iy;
	double qq; //wave number squared
	double wq,wqc;	//linear operator in explicit scheme
	
	/*determine all of the factors used in the Fourier operations*/ 
	/*variables initialized before fast fourier transform factors*/
	facx=2*PI/Nx;
	facy=2*PI/Ny;
	rli=8./(Ny*Nx);   
	enfac=1.0/(Ny*Nx);
	qfac=2.0/(dx*dx);
	qcon=1.0/(dx*dx);
	fac=rli/8.0;
	
	/*get factors for fourier transform*/
	for (ix=0;ix<local_nlast;ix++){
		ixh=ix+local_last_start;
			for (iy=0;iy<Ny;iy++){
				//laplacian operator in Fourier space del_k -> ~|k|^2
				qq=qcon*(cos(ixh*facx)+cos(iy*facy)+cos(ixh*facx)
				*cos(iy*facy)-3);

				//linear part of functional derivatives
				wq=BL+BX*qq*(2.0+qq);
				wqc=w-RK*qq;

				//operation on field in Fourier space for 
				//linear response
				ql[ix*Ny+iy]=exp(wq*qq*dt)*fac;
				qlc[ix*Ny+iy]=exp(wqc*qq*dt)*fac;
				
				//factor for time increment of fields due to 
				//nonlinear parts
				if(wqc==0)
					qnc[ix*Ny+iy]=qq*dt*fac;
				else
					qnc[ix*Ny+iy]=fac*(exp(qq*wqc*dt)-1.0)/wqc;

				if(wq==0)
					qn[ix*Ny+iy]=qq*dt*fac;
				else
					qn[ix*Ny+iy]=fac*(exp(qq*wq*dt)-1.0)/wq;
           }
    }
	return(0);
}

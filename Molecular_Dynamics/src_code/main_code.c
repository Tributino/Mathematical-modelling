#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 500

int main(){
	
	//declaring variables 
	srand(99);
	int i,j;
	float dt=0.01;
	int intsteps=0;
	float L=10.0;
	float ax[N],ay[N],az[N], vx[N],vy[N], vz[N], x[N],y[N], z[N], vhx[N], vhy[N], vhz[N];
	float F,t=0.0,tmax=100.0;
	float m=1.0,g=9.8,kT=1.0;
	float k=1.0;
	float aF=0.0;
	float Psim;
	float Pthe;
	float Pideal;
	float PI=3.14159265359;
	float r1,r2,r3,r4,gr1,gr2,gr3,rx,ry,rz,dist,Fx,Fy,Fz,box;
	float sigma=1.0; 
	float eps =1.0;
	float rc=pow(2.0,1.0/6.0)*sigma; 

	//open files
	FILE *outo,*outt,*outv, *outg;
	outo=fopen("3d.dat","w");
	outt=fopen("3d.xyz","w");
	outv=fopen("3d.prn","w");
	outg=fopen("3d.dat","w");

	//loop over box size L
	for (box=20; box<31; box++){
		L=box;
		t=0.0;

		// initialization of the particles
		for(i=0;i<N;i++){

	      		// picking uniform random position and verifying particles are not too close
	       		int placed = 0;
			while(!placed){
				x[i]=(float)rand()/RAND_MAX*L-L/2.0;
				y[i]=(float)rand()/RAND_MAX*L-L/2.0;
				z[i]=(float)rand()/RAND_MAX*L-L/2.0;
				int too_close =0;
				for(j=0;j<i;j++){
					rx=x[i]-x[j];
					ry=y[i]-y[j];
		            		rz=z[i]-z[j];
					dist=sqrt(pow(rx,2.0) + pow(ry,2.0) + pow(rz,2.0));
					if(dist<rc){
						too_close =1;
					}
				}
			if(!too_close){
		    		placed =1;
			}
			}

			// velocity initialization, following Boltzmann Distribution (1D)
			r1=(float)rand()/RAND_MAX;
			r2=(float)rand()/RAND_MAX;
			r3=(float)rand()/RAND_MAX;
			r4=(float)rand()/RAND_MAX;
			gr1=sqrt(kT/(m))*sqrt(-2.0*log(r1))*cos(2.0*PI*r2);
			gr2=sqrt(kT/(m))*sqrt(-2.0*log(r1))*sin(2.0*PI*r2);
			gr3=sqrt(kT/(m))*sqrt(-2.0*log(r3))*cos(2.0*PI*r4);
			vx[i]=gr1;
			vy[i]=gr2;
			vz[i]=gr3;
			fprintf(outv,"%f\n%f\n%f\n",vx[i],vy[i],vz[i]);
			ax[i]=0;
			ay[i]=0;
			az[i]=0;
		}

		aF=0;
		intsteps=0;

		//initialization of time step loop
		while(t<tmax){
	
			if(intsteps%10==0){
				fprintf(outt,"%i\n",N);
				fprintf(outt,"title\n");
			}
              
			//First two steps of Velocity Verlet algorithm (calculating half step velocity and position)
			for(i=0;i<N;i++){
				vhx[i]=vx[i]+0.5*ax[i]*dt;
				vhy[i]=vy[i]+0.5*ay[i]*dt;
				vhz[i]=vz[i]+0.5*az[i]*dt;

				x[i]=x[i]+vhx[i]*dt;
				y[i]=y[i]+vhy[i]*dt;
				z[i]=z[i]+vhz[i]*dt;
			}

			for(i=0;i<N;i++){

				ax[i] = 0;
				ay[i] = 0;
				az[i] = 0;

				//calculating repulsive force, and updating acceleration accord to WCA potential when particles are too close
				for(j=0;j<N;j++){
				
					rx=x[i]-x[j];
					ry=y[i]-y[j];
		            		rz=z[i]-z[j];
					dist=sqrt(pow(rx,2.0) + pow(ry,2.0) + pow(rz,2.0));
		            
					if (dist <= rc && i!=j){
						
						F=-4.0*eps*((-12.0*pow(sigma,12.0)/pow(dist,13.0))
						F+=(6.0*pow(sigma,6.0)/pow(dist,7.0)));	
						Fx=F*rx/dist;
						Fy=F*ry/dist;
						Fz=F*rz/dist;
						ax[i] += Fx/m;
						ay[i] += Fy/m;
						az[i] += Fz/m;

					}
				}				
			
				//final step of the velocity verlet algorithm: updating final/full velocity 
				vx[i]=vhx[i]+0.5*ax[i]*dt;
				vy[i]=vhy[i]+0.5*ay[i]*dt;
				vz[i]=vhz[i]+0.5*az[i]*dt;
				
				//boundary conditions - reflecting velocity whenever particles hit boundary
				//af is the force of all colisions - when particles hit box 
				if(x[i]<-L/2.0 || x[i]>L/2.0){
					vx[i]=-vx[i];
					aF+=2.0*m*fabs(vx[i])/dt;
				}
				if(y[i]<-L/2.0 || y[i]>L/2.0){
					vy[i]=-vy[i];
					aF+=2.0*m*fabs(vy[i])/dt;
				}
				if(z[i]<-L/2.0 || z[i]>L/2.0){
					vz[i]=-vz[i];
					aF+=2.0*m*fabs(vz[i])/dt;
				}
				if(intsteps%10==0)
					fprintf(outt,"a%i %f %f %f\n",i,x[i],y[i],z[i]);

			}

			t+=dt;
			intsteps++;

		}

	
		//Pressures calculation and printing;
		
		aF=aF/(float)intsteps;
		Psim=aF/(6.0*L*L);
		Pthe=N*kT/((L*L*L)-4.0*N*((4.0/3.0)*PI*pow((sigma/2.0),3.0)));
		Pideal=N*kT/(L*L*L);
		
		printf("P from sims: %f        P theoretical gas excluded volume: %f\n          P from ideal gas law: %f\n",Psim,Pthe,Pideal);		
		fprintf(outg,"%f %f\n",Psim, (L*L*L));
	}

}

// Read .trr (or .xtc) file and obtain rmsd
#include <limits> 		// DBL_MAX
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <cstring>
#include <fstream>
#include <ctime>
#include <malloc.h>
#include </usr/include/xdrfile/xdrfile_trr.h> // xdr include file 
#include </usr/include/xdrfile/xdrfile_xtc.h> // xdr include file 


    using namespace std;

int main(int argc, char **argv){
	
	/*variables*/
	// XTC variables
	XDRFILE *xd;		// the xtc file
	int natoms;	// number of total atoms
	int step;		// the step counter for the simulation
	float time;		// simulation time
	matrix box;		// box coordinates in a 3x3 matrix
	rvec *coor;		// atom coordinates in a 2-D matrix
	rvec *vel;		// atom coordinates in a 2-D matrix
	rvec *force;		// atom coordinates in a 2-D matrix
//	float prec;	
	float lambda;
	
	/*read xtc files*/
	read_trr_natoms(argv[1],&natoms);

	int firstframe = 10;				
	if(argc >= 3 ) firstframe = atoi(argv[2]);
	int lastframe = -1;		// ps
	if(argc >= 4 ) lastframe = atoi(argv[3]);

	int i,j,k;

	//other variables
	// output	
	char f[40];
	FILE *data;
	
	// user defined variables

	int nFrame = 0;
	double ***position;		// position [nFrame][natoms][dimension]
	double *rmsd;
	double dt = 0.0;
	int natomsPerMol = natoms/10;	// number of atoms in a molecule

	coor = (rvec *)malloc(natoms*sizeof(rvec));
	vel = (rvec *)malloc(natoms*sizeof(rvec));
	force = (rvec *)malloc(natoms*sizeof(rvec));

	//open trr(xtc) files, get the number of frames
	xd=xdrfile_open(argv[1],"r");
	k = 0;
	while( ! read_trr(xd, natoms, &step, &time, &lambda, box, coor, vel, force)){
		if(coor == 0){
			printf("Insufficient memory to load .trr file. \n");
			return 1;
		}
		if(dt == 0 && time > 0 ) dt = time;
	    if(time >= firstframe && (time <= lastframe || lastframe <= firstframe)  ) k++;
	}

	nFrame = k;
	printf("# Counted frame = %d...\n",nFrame);
	printf("# dt = %.2f...\n#...\n",dt);

	// allocate memory
	
	rmsd = (double *)malloc(sizeof(double)*nFrame/10);	
	for(i=0;i<nFrame/10;i++) rmsd[i] = 0.0;

	position = (double ***)malloc(sizeof(double **)*nFrame);
	for(i=0;i<nFrame;i++){
		position[i] = (double **)malloc(sizeof(double*)*natoms/natomsPerMol);
		for(j=0;j<natoms/natomsPerMol;j++) position[i][j] = (double *)malloc(sizeof(double)*3);
	}
	if(position == NULL){
		printf("# Error in memory allocation...\n");
		return 1;
	}
	
	// second loop to restore data
	xd=xdrfile_open(argv[1],"r");
	i = 0;
	while( ! read_trr(xd, natoms, &step, &time, &lambda, box, coor, vel, force)){
	    if(time >= firstframe && (time <= lastframe || lastframe <= firstframe)  ){
			for(j=0;j<natoms/natomsPerMol;j++){
				for(k=0;k<3;k++) position[i][j][k] = coor[j*natomsPerMol][k];
			}
			i++;
		}
	}

	printf("# Finish reading .trr (.xtc) file...(%d frames)...\n",i);

	// obtain rmsd(t)
	int count;
	for(int t=1;t<nFrame/10;t++){
		count = 0;
		for(i=0;i+t<nFrame;i++){
			count++;
			for(j=0;j<natoms/natomsPerMol;j++){
				for(k=0;k<3;k++){
					rmsd[t] += (position[i][j][k] - position[i+t][j][k])*(position[i][j][k] - position[i+t][j][k])/double(natoms/natomsPerMol);
				}
				if(rmsd[t] >= std::numeric_limits<double>::max()) printf("# error ( DBL_limit)\n");
			}
		}
		rmsd[t] =  rmsd[t]/double(count) ;
	}


	//output result
	sprintf(f,"rmsd.dat");
	data = fopen(f,"w");
	for(i = 1;i<nFrame/10;i++){
		fprintf(data,"%.3f\t%e\n",i*dt,rmsd[i]);
	}
	fclose(data);

	return 0;
}

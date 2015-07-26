// vim: set cin ts=4 sw=4 tw=80:
// g++ -Wall -O3 E-replicaA-partition-l.cc -o replicaA
// Run with command line arguments, e.g. ./testo betamin betamax delbeta
// Considering 2d Blume Capel model
//warming up system for first N_mc/10 loops
//averaging energy for the next N_mc updates
//incorporating SIMULATED ANNEALING

// axis2 = axis1
//axis1 × axis1 square lattice, cut into 2 rectangles of
//  sizes ell x axis1 & (axis1-ell) x axis1

#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <boost/multi_array.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/lexical_cast.hpp>
#include "Eigen/Core"
#include <vector>

// gen is a variable name
// Its data-type is boost::random::mt19937
boost::random::mt19937 gen;
// boost::random::mt19937 gen(std::time(0)); 
// time(0) changes seed every time you run
using boost::lexical_cast;
using boost::bad_lexical_cast;
using namespace std;

typedef
 boost::multi_array < int, 2 > array_2d;
// typedef keyword allows you to create an alias fo a data type

// Magnitude of parameters
double J = 1.0;
double D = 1.9;

unsigned int axis1 = 0;
unsigned int axis2 = axis1;
// above assigns length along each dimension of the 2d configuration

//No.of Monte Carlo updates we want
unsigned int N_mc = 1e7;
//No. of measurements per write to fil
unsigned int N_meas = 1e5;

//Function templates
int roll_coin(int a, int b);
double random_real(int a, int b);
double energy_tot(array_2d sitespin);
//double mag_tot(array_2d sitespin);
double nn_energy(array_2d sitespin, unsigned int row, unsigned int col);

long double calc_ratio(const array_2d& s1, const array_2d& s2, const double beta, const int ell);
void getBetas(vector<double>& betas, vector<int>& measure);
void swapSims(int sim1, int sim2, vector<array_2d*>& ss1, vector<array_2d*>& ss2, vector<double>& energy);
void clusterMove(array_2d& s1, array_2d& s2, const double beta, const int ell);

int main(int argc, char const * argv[])
{
	if (argc != 3)
	{
		cout << "Expecting two inputs: axis1, ell."
		     << endl << "Got " << argc - 1 << endl;
		return 1;
	}	

	unsigned int ell(0);

	try
	{
		axis1 = lexical_cast<unsigned int>(argv[1]);
		ell = lexical_cast<unsigned int>(argv[2]);
	}
	catch (const bad_lexical_cast & x)
	{
		cout << "Cannot convert input for axis1 to unsigned int" << endl;
		return 2;
	}
	
	axis2 = axis1;

	string axis_str = lexical_cast<string>(axis1);
	string ell_str = lexical_cast<string>(ell);
	//ofstream fout(string("EmA" + axis_str + "p" + ell_str + ".dat").c_str());	
	//ofstream fratio(string("Ratio" + axis_str + "p" + ell_str + ".dat").c_str());	
    // Opens a file for output
    //
    
    vector<double> betas;
    vector<int> meas;
    vector<ofstream*> fouts;
    getBetas(betas, meas);
    int numB = betas.size();


	//define replica 1 spin configuration array
	array_2d tsitespin1(boost::extents[axis1][axis2]);
	//define replica 2 spin configuration array
	array_2d tsitespin2(boost::extents[axis1][axis2]);

	//For subsystem A,both replicas have same spin configuration
	for (unsigned int i = 0; i < axis1; ++i)
    {   for (unsigned int j = 0; j < axis2; ++j)
        {
        tsitespin1[i][j] = roll_coin(-1, 1);
        tsitespin2[i][j] = tsitespin1[i][j];
        }
    }

    vector<array_2d*> ss1;
    vector<array_2d*> ss2;
    vector<double> energy;
    vector<long double> ratio_sum;

    // For simplicity, push back the same random copy to all the temperatures
    for(int b=0;b<numB;b++){
        ss1.push_back(new array_2d(tsitespin1));
        ss2.push_back(new array_2d(tsitespin2));
        energy.push_back(2.0*energy_tot(tsitespin1));
        fouts.push_back(new ofstream(string("ratio_" + lexical_cast<string>(betas[b])).c_str()));
        ratio_sum.push_back(0);
    }
    //std::cout << ss1[0] << std::endl;
    //std::cout << (*ss1[0])[0][0] << std::endl;
    //return 1;

	//calculate avg energy for replica spin config at temp 1/beta
	//logic: for a[n1][n2], a[n1] is n1 copies of 1d array of length n2

    //Main simulation loop
    unsigned int sys_size = axis1 * axis2;
    unsigned int row, col, label;
    double r(0), acc_ratio(0) ;
    double en_sum(0);
    int spin(0),newspin(0),choice[2]={0,0},choice_ind;
    unsigned int eff_sys_size = 2*sys_size - ell * axis2 ;


    for (unsigned int i = 1; i <=1e5+N_mc; ++i){
        for(int b=0;b<numB;b++){
			for (unsigned int j = 1; j <=eff_sys_size; ++j){
                //Choose a random spin site for the entire 2 replica system
                double energy_diff(0);
                label = roll_coin(1,2*sys_size);

                //if the random spin site is located in layer 1
                if (label <= sys_size)
                {
                    if (label % axis2 == 0)
                    {
                        row = (label / axis2) - 1;
                        col = axis2 -1 ;
                    }
                    else
                    {
                        col = label % axis2 - 1;
                        row = (label-col-1)/axis2;
                    }

                    spin = (*ss1[b])[row][col];
                    if (spin==0) 
                    {	choice[0]=-1;
                        choice[1]=1;
                    }
                    if (spin==-1) 
                    {	choice[0]=0;
                        choice[1]=1;
                    }
                    if (spin==1) 
                    {	choice[0]=-1;
                        choice[1]=0;
                    }

                    choice_ind = roll_coin(0,1);
                    newspin = choice[choice_ind];
     
                    energy_diff =-nn_energy(*ss1[b],row,col);
                    energy_diff -=D*spin*spin;
                    (*ss1[b])[row][col]=newspin;
                    energy_diff +=nn_energy(*ss1[b],row,col);
                    energy_diff +=D*newspin*newspin;


                    if (row < ell)
                      { energy_diff -=nn_energy(*ss2[b],row,col);
                        energy_diff -=D*spin*spin;
                        (*ss2[b])[row][col]=newspin;
                        energy_diff +=nn_energy(*ss2[b],row,col);
                        energy_diff +=D*newspin*newspin;
                      }

                    //Generate a random no. r such that 0 < r < 1
                    r = random_real(0, 1);
                    acc_ratio = exp(-1.0 * energy_diff *betas[b]);

                    //Spin flipped if r <= acceptance ratio
                    if (r <= acc_ratio)
                    {

                            energy[b] += energy_diff;
                    }
                    else 
                    {   (*ss1[b])[row][col] =spin;
                        if (row < ell)
                            (*ss2[b])[row][col]=spin;
                    }
                }

                //if the random spin site is located in layer 2
                if (label > sys_size)
                {
                    label -= sys_size;

                    if (label % axis2 == 0)
                    {
                        row = (label / axis2) - 1;
                        col = axis2 -1 ;
                    }
                    else
                    {
                        col = label % axis2 - 1;
                        row = (label-col-1)/axis2;
                    }

                    spin = (*ss2[b])[row][col];
                     if (spin==0) 
                        {	choice[0]=-1;
                            choice[1]=1;
                        }
                     if (spin==-1) 
                        {	choice[0]=0;
                            choice[1]=1;
                        }
                     if (spin==1) 
                        {	choice[0]=-1;
                            choice[1]=0;
                        }

                    choice_ind = roll_coin(0,1);
                    newspin = choice[choice_ind];
     
                    energy_diff =-nn_energy(*ss2[b],row,col);
                    energy_diff -=D*spin*spin;
                    (*ss2[b])[row][col]=newspin;
                    energy_diff +=nn_energy(*ss2[b],row,col);
                    energy_diff +=D*newspin*newspin;


                    if (row < ell)
                      { energy_diff -=nn_energy(*ss1[b],row,col);
                        energy_diff -=D*spin*spin;
                        (*ss1[b])[row][col]=newspin;
                        energy_diff +=nn_energy(*ss1[b],row,col);
                        energy_diff +=D*newspin*newspin;
                      }


                    //Generate a random no. r such that 0 < r < 1
                    r = random_real(0, 1);
                    acc_ratio = exp(-1.0 * energy_diff *betas[b]);

                    //Spin flipped if r <= acceptance ratio
                    if (r <= acc_ratio)
                    {

                            energy[b] += energy_diff;
                    }
                    else 
                    {   (*ss2[b])[row][col] =spin;
                        if (row < ell)
                            (*ss1[b])[row][col]=spin;
                    }
                }
            }
            // We'll do axis 1 cluster moves, just as an attempt
            /*
			for (unsigned int j = 1; j <=axis1; ++j){
                std::cout << betas[b] << ", " << j << " ,Energy = " << energy[b] << std::endl;
                clusterMove(*ss1[b], *ss2[b], betas[b], ell);
                // Need to fully recalculate energy after these moves
                energy[b] = energy_tot(*ss1[b]) + energy_tot(*ss2[b]);
            }
            */
            //if (i> 1e5) en_sum += energy;
            if (meas[b] == 1){
                if (i > N_meas){
                    ratio_sum[b] += calc_ratio(*ss1[b], *ss2[b], betas[b], ell);
                    if ((i % N_meas) == 0){
                        (*fouts[b]) << ratio_sum[b] / N_meas << endl;
                        ratio_sum[b] = 0;
                    }
                }
            }
        }
        //fout << beta << '\t' << en_sum / N_mc << endl;
        //fratio << beta << '\t' << ratio_sum / N_mc << endl;
        //Parallel tempering loop
        for(int b=0;b<(numB*numB);b++){
            int t = roll_coin(0,numB-2);
            double dB = betas[t+1] - betas[t];
            double dE = energy[t+1] - energy[t];
            r = random_real(0, 1);
            acc_ratio = exp(dB*dE);
            //std::cout << "r     " << r << std::endl;
            //std::cout << "acc   " << acc_ratio << std::endl;
            if(r < acc_ratio){
                //std::cout << "e1 " << energy[t] << std::endl;
                //std::cout << "e2 " << energy[t+1] << std::endl;
                swapSims(t, t+1, ss1, ss2, energy);
                //std::cout << "e1 " << energy[t] << std::endl;
                //std::cout << "e2 " << energy[t+1] << std::endl;
            }
        }
    }


	return 0;
}



//function to generate random integer
// between 2 integers a & b, including a & b
int roll_coin(int a, int b)
{
	boost::random::uniform_int_distribution <> dist(a, b);

	return dist(gen);

}

//function to generate random real no.
// between 2 integers a & b, including a & excluding b

double random_real(int a, int b)
{
	boost::random::uniform_real_distribution <> dist(a, b);
	// uniform_real_distribution: continuous uniform distribution 
	//on some range [min, max) of real number
	return dist(gen);

}

//function to calculate total energy
//for a given spin configuration
//with open boundary conditions

double energy_tot(array_2d sitespin)
{
	double energy = 0;

	for (unsigned int i = 0; i < axis1 - 1; ++i)
	{
		for (unsigned int j = 0; j < axis2 - 1; ++j)
		{
			energy +=-J*sitespin[i][j]*sitespin[i+1][j];

			energy +=-J * sitespin[i][j] * sitespin[i][j + 1];
			energy += D*sitespin[i][j]*sitespin[i][j];

		}
	}

	//open boundary conditions -- end spins disconnected
	

	return energy;
}


//Calculating interaction energy change for spin 1
//at random site->(row,col) with its nearest neighbours
double nn_energy(array_2d sitespin, unsigned int row, unsigned int col)
{
	double nn_en = 0;

	if (row > 0 && row < axis1 - 1)
	{
		nn_en -=J*sitespin[row][col]*sitespin[row-1][col];
		nn_en -=J*sitespin[row][col]*sitespin[row+1][col];
	}

	if (col > 0 && col < axis2 - 1)
	{
		nn_en -=J*sitespin[row][col]*sitespin[row][col-1];
		nn_en -=J*sitespin[row][col]*sitespin[row][col+1];
	}

	if (row == 0)
	{
		// open b.c. -- don't need this: nn_en -=J*sitespin[0][col]*sitespin[axis1-1][col];
		nn_en -=J*sitespin[0][col]*sitespin[1][col];

	}

	if (row == axis1 - 1)
	{
		nn_en -=J*sitespin[axis1-1][col]*sitespin[axis1-2][col];
		// open b.c. -- don't need this:nn_en -=J*sitespin[axis1-1][col]*sitespin[0][col];

	}

	if (col == 0)
	{
		// open b.c. -- don't need this:nn_en -=J*sitespin[row][0]*sitespin[row][axis2-1];
		nn_en -=J*sitespin[row][0]*sitespin[row][1];

	}

	if (col == axis2 - 1)
	{
		nn_en -=J*sitespin[row][axis2-1]*sitespin[row][axis2-2];
		// open b.c. -- don't need this:nn_en -=J*sitespin[row][axis2-1]*sitespin[row][0];

	}
	return nn_en;
}


// Method for calcuiating the ratio of partition functions using the transfer matrix method
// ell contains the current number of spins in region A
long double calc_ratio(const array_2d& s1, const array_2d& s2, const double beta, const int ell){
    Eigen::Matrix<double, 3,3> tmat; // Temporary matrix for filling and multiplying to results
    Eigen::Matrix<double, 3,3> t_top; // Spins in top layer
    Eigen::Matrix<double, 3,3> t_bot; // Spins in bottom layer
    Eigen::Matrix<double, 3,3> t_con; // Spins with both layers connected

    t_top.setIdentity(3,3);
    t_bot.setIdentity(3,3);
    t_con.setIdentity(3,3);

    long double ttop = 0.0;
    long double tbot = 0.0;
    long double tcon = 0.0;

    long double coeff = 0.0;

    // Number of adjacent -1 (m), zero (0), and +1 (p) spins to the 1st and 2nd spin of the transfer matrix
    int nm_top1 = 0;
    int np_top1 = 0;
    int nm_bot1 = 0;
    int np_bot1 = 0;

    int nm_top2 = 0;
    int np_top2 = 0;
    int nm_bot2 = 0;
    int np_bot2 = 0;

    // We will do the calculation along a fixed row, for all the spins in the column
    // The notation for accessing spins is s1[row][column]
    // Neighbours are +-1 in each direction, plus periodic boundary conditions
    // axis1 is global and defines the length of the system in the row direction
    // axis2 defines the length in the column direction
    
    // ell will be from 0 to L-1 (it could be L, but since we are always adding spins that isn't necessary)
    // We are doing the transfer matrix over the spins in row=ell
    int row_up = (ell+1)%axis1;
    int row_down = (ell-1+axis1)%axis1;

    for(int col=0; col<axis2; col++){
        if(col==0){
            nm_top1 = 0;
            np_top1 = 0;
            nm_bot1 = 0;
            np_bot1 = 0;

            if(s1[row_up][col]==-1) nm_bot1++;
            else if(s1[row_up][col]==1) np_bot1++;

            if(s1[row_down][col]==-1) nm_bot1++;
            else if(s1[row_down][col]==1) np_bot1++;

            if(s2[row_up][col]==-1) nm_top1++;
            else if(s2[row_up][col]==1) np_top1++;

            if(s2[row_down][col]==-1) nm_top1++;
            else if(s2[row_down][col]==1) np_top1++;
        }
        else{
            nm_top1 = nm_top2;
            np_top1 = np_top2;
            nm_bot1 = nm_bot2;
            np_bot1 = np_bot2;
        }
        // First set the occupation of the spins in both layers
        nm_top2 = 0;
        np_top2 = 0;
        nm_bot2 = 0;
        np_bot2 = 0;

        int col2 = (col+1)%axis2;

        if(s1[row_up][col2]==-1) nm_bot2++;
        else if(s1[row_up][col2]==1) np_bot2++;

        if(s1[row_down][col2]==-1) nm_bot2++;
        else if(s1[row_down][col2]==1) np_bot2++;

        if(s2[row_up][col2]==-1) nm_top2++;
        else if(s2[row_up][col2]==1) np_top2++;

        if(s2[row_down][col2]==-1) nm_top2++;
        else if(s2[row_down][col2]==1) np_top2++;

        // Now we know the number of adjacent spins (in the neighboring rows, ignoring the columns) for the spin in position col and col+1
        // We will split half of the interaction term to each spin, since for each bond in the column that we are doing the transfer matrix over
        // will have each physical spin occuring twice
        //
        // The two energy terms we have are 
        // -J \sum_{<i,j>} s_i s_j
        // +D s_i^2
        // Table of energies:
        // s1   s2      E
        // 1    1       2D - J
        // 1    0       D
        // 1    -1      2D + J
        // 0    1       D
        // 0    0       0
        // 0    -1      D
        // -1   1       2D + J
        // -1   0       D
        // -1   -1      2D - J
        // 
        // Keep in mind that we split the D contributions up into the two
        // transfer matrices, while the J term of the interaction is only
        // in one matrix. Contributions from external spins are also split up
        // over two matrices.
        //
        // tmat is a (3,3) matrix, where the rows correspond to spin col = -1,0,1
        // and the columns correspond to spin col2 = -1,0,1
        // We calculate the energy using the above quantities, then exponentiate the terms to get the weights of each configuration
        // By taking the product of matrices and finally a trace, we effectively integrate over the spin chain with the rest of the system frozen

        // First the calculation for the bottom layer in isolation
        // -1, -1 state
        tmat(0,0) = exp(-0.5*beta*(D*2 - J*(2 + nm_bot1 - np_bot1 + nm_bot2 - np_bot2)));
        // -1, 0 state
        tmat(0,1) = exp(-0.5*beta*(D*1 - J*(0 + nm_bot1 - np_bot1)));
        // -1, 1 state
        tmat(0,2) = exp(-0.5*beta*(D*2 - J*(-2 + nm_bot1 - np_bot1 + np_bot2 - nm_bot2)));
        // 0, -1 state
        tmat(1,0) = exp(-0.5*beta*(D*1 - J*(0 + nm_bot2 - np_bot2)));
        // 0, 0 state
        tmat(1,1) = exp(-0.5*beta*(D*0 - J*(0)));
        // 0, 1 state
        tmat(1,2) = exp(-0.5*beta*(D*1 - J*(0 + np_bot2 - nm_bot2)));
        // 1, -1 state
        tmat(2,0) = exp(-0.5*beta*(D*2 - J*(-2 + np_bot1 - nm_bot1 + nm_bot2 - np_bot2)));
        // 1, 0 state
        tmat(2,1) = exp(-0.5*beta*(D*1 - J*(0 + np_bot1 - nm_bot1)));
        // 1, 1 state
        tmat(2,2) = exp(-0.5*beta*(D*2 - J*(2 + np_bot1 - nm_bot1 + np_bot2 - nm_bot2)));

        t_bot *= tmat;
        coeff = t_bot.maxCoeff();
        tbot += log(coeff);
        t_bot /= coeff;

        // Now the top layer in isolation
        // -1, -1 state
        tmat(0,0) = exp(-0.5*beta*(D*2 - J*(2 + nm_top1 - np_top1 + nm_top2 - np_top2)));
        // -1, 0 state
        tmat(0,1) = exp(-0.5*beta*(D*1 - J*(0 + nm_top1 - np_top1)));
        // -1, 1 state
        tmat(0,2) = exp(-0.5*beta*(D*2 - J*(-2 + nm_top1 - np_top1 + np_top2 - nm_top2)));
        // 0, -1 state
        tmat(1,0) = exp(-0.5*beta*(D*1 - J*(0 + nm_top2 - np_top2)));
        // 0, 0 state
        tmat(1,1) = exp(-0.5*beta*(D*0 - J*(0)));
        // 0, 1 state
        tmat(1,2) = exp(-0.5*beta*(D*1 - J*(0 + np_top2 - nm_top2)));
        // 1, -1 state
        tmat(2,0) = exp(-0.5*beta*(D*2 - J*(-2 + np_top1 - nm_top1 + nm_top2 - np_top2)));
        // 1, 0 state
        tmat(2,1) = exp(-0.5*beta*(D*1 - J*(0 + np_top1 - nm_top1)));
        // 1, 1 state
        tmat(2,2) = exp(-0.5*beta*(D*2 - J*(2 + np_top1 - nm_top1 + np_top2 - nm_top2)));

        t_top *= tmat;
        coeff = t_top.maxCoeff();
        ttop += log(coeff);
        t_top /= coeff;

        // Now the calculation of both layers connected
        // -1, -1 state
        tmat(0,0) = exp(-0.5*beta*(D*4 - J*(4 + nm_bot1 - np_bot1 + nm_bot2 - np_bot2 + nm_top1 - np_top1 + nm_top2 - np_top2)));
        // -1, 0 state
        tmat(0,1) = exp(-0.5*beta*(D*2 - J*(0 + nm_bot1 - np_bot1 + nm_top1 - np_top1)));
        // -1, 1 state
        tmat(0,2) = exp(-0.5*beta*(D*4 - J*(-4 + nm_bot1 - np_bot1 + np_bot2 - nm_bot2 + nm_top1 - np_top1 + np_top2 - nm_top2)));
        // 0, -1 state
        tmat(1,0) = exp(-0.5*beta*(D*2 - J*(0 + nm_bot2 - np_bot2 + nm_top2 - np_top2)));
        // 0, 0 state
        tmat(1,1) = exp(-0.5*beta*(D*0 - J*(0)));
        // 0, 1 state
        tmat(1,2) = exp(-0.5*beta*(D*2 - J*(0 + np_bot2 - nm_bot2 + np_top2 - nm_top2)));
        // 1, -1 state
        tmat(2,0) = exp(-0.5*beta*(D*4 - J*(-4 + np_bot1 - nm_bot1 + nm_bot2 - np_bot2 + np_top1 - nm_top1 + nm_top2 - np_top2)));
        // 1, 0 state
        tmat(2,1) = exp(-0.5*beta*(D*2 - J*(0 + np_bot1 - nm_bot1 + np_top1 - nm_top1)));
        // 1, 1 state
        tmat(2,2) = exp(-0.5*beta*(D*4 - J*(4 + np_bot1 - nm_bot1 + np_bot2 - nm_bot2 + np_top1 - nm_top1 + np_top2 - nm_top2)));

        t_con *= tmat;
        coeff = t_con.maxCoeff();
        tcon += log(coeff);
        t_con /= coeff;
    }
    ttop += log(t_top.trace());
    tbot += log(t_bot.trace());
    tcon += log(t_con.trace());

    return exp(tcon - ttop - tbot);
    //return (t_con.trace() / (t_top.trace() * t_bot.trace()));
}

// Reads betas.dat and reads in the set of betas as well as which we should
// measure
void getBetas(vector<double>& betas, vector<int>& measure){
    std::string filename = "betas.dat";
    int m;
    double b;
    std::fstream inFile(filename.c_str());
    betas.resize(0);
    measure.resize(0);
    while(inFile >> b >> m){
        betas.push_back(b);
        measure.push_back(m);
    }
}

// Swap simulations sim1 and sim2 after a successfull parallel tempering move
void swapSims(int sim1, int sim2, vector<array_2d*>& ss1, vector<array_2d*>& ss2, vector<double>& energy){
    array_2d* temp1 = ss1[sim1];
    array_2d* temp2 = ss2[sim1];
    ss1[sim1] = ss1[sim2];
    ss2[sim1] = ss2[sim2];
    ss1[sim2] = temp1;
    ss2[sim2] = temp2;
    double te = energy[sim1];
    energy[sim1] = energy[sim2];
    energy[sim2] = te;
}

void addToCluster(int z, int sx, int sy, int l, double Padd, array_2d& f1, array_2d& f2, array_2d& s1, array_2d& s2, int ell){
    // First check if the spin matches
    // Then check if the spin is in the list already
    //std::cout << "Trying to add: " << l << ", " << sx << ", " << sy << std::endl;
    double r = 0.0;
    if (l==0){
        if(s1[sx][sy] != z) return;
        if(f1[sx][sy] != 1) return;
        r = random_real(0, 1);
        if (r<Padd){
            f1[sx][sy] = -1;
            if(sx<ell){
                f2[sx][sy] = -1;
            }
        }
    }
    else{
        if(s2[sx][sy] != z) return;
        if(f2[sx][sy] != 1) return;
        r = random_real(0, 1);
        if (r<Padd){
            f2[sx][sy] = -1;
            if(sx<ell){
                f1[sx][sy] = -1;
            }
        }
    }
    if (r<Padd){
        // Adding neighbors
        if ((sx-1) < 0){
            addToCluster(z, sx-1+axis1, sy, l, Padd, f1, f2, s1, s2, ell);
        }
        else{
            addToCluster(z, sx-1, sy, l, Padd, f1, f2, s1, s2, ell);
        }
        if ((sx+1) >= axis1){
            addToCluster(z, sx+1-axis1, sy, l, Padd, f1, f2, s1, s2, ell);
        }
        else{
            addToCluster(z, sx+1, sy, l, Padd, f1, f2, s1, s2, ell);
        }
        if ((sy-1) < 0){
            addToCluster(z, sx, sy-1+axis2, l, Padd, f1, f2, s1, s2, ell);
        }
        else{
            addToCluster(z, sx, sy-1, l, Padd, f1, f2, s1, s2, ell);
        }
        if ((sy+1) >= axis2){
            addToCluster(z, sx, sy+1-axis2, l, Padd, f1, f2, s1, s2, ell);
        }
        else{
            addToCluster(z, sx, sy+1, l, Padd, f1, f2, s1, s2, ell);
        }
        if(sx<ell){
            // Adding neighbors in the other layer
            l = (l+1)%2;
            if ((sx-1) < 0){
                addToCluster(z, sx-1+axis1, sy, l, Padd, f1, f2, s1, s2, ell);
            }
            else{
                addToCluster(z, sx-1, sy, l, Padd, f1, f2, s1, s2, ell);
            }
            if ((sx+1) >= axis1){
                addToCluster(z, sx+1-axis1, sy, l, Padd, f1, f2, s1, s2, ell);
            }
            else{
                addToCluster(z, sx+1, sy, l, Padd, f1, f2, s1, s2, ell);
            }
            if ((sy-1) < 0){
                addToCluster(z, sx, sy-1+axis2, l, Padd, f1, f2, s1, s2, ell);
            }
            else{
                addToCluster(z, sx, sy-1, l, Padd, f1, f2, s1, s2, ell);
            }
            if ((sy+1) >= axis2){
                addToCluster(z, sx, sy+1-axis2, l, Padd, f1, f2, s1, s2, ell);
            }
            else{
                addToCluster(z, sx, sy+1, l, Padd, f1, f2, s1, s2, ell);
            }
        }
    }
}

// Cluster move that builds a Wolff cluster flipping S=1 <--> S=-1
// Does not interact with S=0 states (can't solve EQs this way)
void clusterMove(array_2d& s1, array_2d& s2, const double beta, const int ell){
    // Choose a random spin
    // If it is not S=1 or S=-1, finish
    // If it is, build the cluster, adding neighbors with probability P_add if
    // they match the spin configuration
    // Once the cluster has stopped building, flip all the spins in the cluster
    double Padd = 1 - exp(-2*beta*J);
    int start_spin = roll_coin(0,axis1*axis2*2-1);
    int start_layer = start_spin/(axis1*axis2);
    start_spin = start_spin % (axis1*axis2);
    int sx, sy;
    sx = start_spin % axis1;
    sy = start_spin / axis1;
    int spin_conf = 0;
    if (start_layer == 0){
        if(s1[sx][sy] == 0) return;
        spin_conf = s1[sx][sy];
    }
    else{
        if(s2[sx][sy] == 0) return;
        spin_conf = s2[sx][sy];
    }
	array_2d f1(boost::extents[axis1][axis2]);
	array_2d f2(boost::extents[axis1][axis2]);
	for (unsigned int i = 0; i < axis1; ++i)
    {   for (unsigned int j = 0; j < axis2; ++j)
        {
            f1[i][j] = 1;
            f2[i][j] = 1;
        }
    }
    addToCluster(spin_conf, sx, sy, start_layer, Padd, f1, f2, s1, s2, ell);
    // When this is done, f1 and f2 contain the spins to flip
	for (unsigned int i = 0; i < axis1; ++i)
    {   for (unsigned int j = 0; j < axis2; ++j)
        {
            s1[i][j] = s1[i][j]*f1[i][j];
            s2[i][j] = s2[i][j]*f2[i][j];
        }
    }
}

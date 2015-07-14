// vim: set cin ts=4 sw=4 tw=80:
// g++ -Wall -O3 Ereplica-vs-beta.cc -o testo
// Run with command line arguments, e.g. ./testo betamin betamax delbeta
// Considering 2d Blume Capel model
//warming up system for first N_mc/10 loops
//averaging energy for the next N_mc updates
//incorporating SIMULATED ANNEALING
//Previous result: kT/J=0.695, D/J=1.965 at critical point

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

// Magnitude of J
double J = 1.0;
double D = 1.965;
unsigned int axis1 = 0;
unsigned int axis2 = axis1;
// above assigns length along each dimension of the 2d configuration

//No.of Monte Carlo updates we want
unsigned int N_mc = 1e6;

// Size of the region, must be set by command line
// Corresponds to the number of STRIPS in the system
// So this value runs from 0 to axis1, inclusive
unsigned int regionSize = 0;

//Function templates
int roll_coin(int a, int b);
double random_real(int a, int b);
double energy_tot(array_2d sitespin);
//double mag_tot(array_2d sitespin);
double nn_energy(array_2d sitespin, unsigned int row, unsigned int col);
long double calc_ratio(const array_2d& s1, const array_2d& s2, const double beta);

bool inA(const int row, const int col=0){
    if (row<regionSize){
            return true;
    }
    return false;
}

int main(int argc, char const * argv[])
{
	if (argc != 6)
	{
		cout << "Expecting five inputs: bmin, bmax, dbeta, axis1, regionSize"
		     << endl << "Got " << argc - 1 << endl;
		return 1;
	}	

    double beta(0);
    double bmin(0);
    double bmax(0);
    double dbeta(0);

	try
	{
		bmin = lexical_cast<double>(argv[1]);
		bmax = lexical_cast<double>(argv[2]);
		dbeta = lexical_cast<double>(argv[3]);
	}
	catch (const bad_lexical_cast & x)
	{
		cout << "Cannot convert input to double" << endl;
		return 2;
	}

	try
	{
		axis1 = lexical_cast<unsigned int>(argv[4]);
	}
	catch (const bad_lexical_cast & x)
	{
		cout << "Cannot convert input for axis1 to unsigned int" << endl;
		return 2;
	}

	try
	{
		regionSize = lexical_cast<unsigned int>(argv[5]);
	}
	catch (const bad_lexical_cast & x)
	{
		cout << "Cannot convert input for regionSize to unsigned int" << endl;
		return 2;
	}
	
	axis2 = axis1;

	//string axis_str = lexical_cast<string>(axis1);
	//ofstream fout(string("Em" + axis_str + ".dat").c_str());	
    // Opens a file for output
	

	//define replica 1 spin configuration array
	array_2d sitespin1(boost::extents[axis1][axis2]);
	//define replica 2 spin configuration array
	array_2d sitespin2(boost::extents[axis1][axis2]);

	//For subsystem A,both replicas have same spin configuration
	for (unsigned int i = 0; i < axis1; ++i)
		{   for (unsigned int j = 0; j < axis2; ++j)
			{
			sitespin1[i][j] = roll_coin(-1, 1);
			sitespin2[i][j] = sitespin1[i][j];
			}
		}


	double energy =2.0*energy_tot(sitespin1);

	//calculate avg energy for replica spin config at temp 1/beta
	//logic: for a[n1][n2], a[n1] is n1 copies of 1d array of length n2


	for (beta =bmin;beta<bmax;beta += dbeta)
	{
        string beta_str = lexical_cast<string>(beta);
        ofstream fratio(string("ratio_" + beta_str).c_str());

		unsigned int sys_size = axis1 * axis2;
		unsigned int row, col, label;
		double r(0), acc_ratio(0) ;
	
		double en_sum(0);
		long double ratio_sum(0);
        int ratio_counter(0);
		int spin(0),newspin(0),choice[2]={0,0},choice_ind;


		for (unsigned int i = 1; i <=1e5+N_mc; ++i)
		{
			for (unsigned int j = 1; j <=2*sys_size; ++j)
			{
				
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

        		spin = sitespin1[row][col];
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
 
				energy_diff =-nn_energy(sitespin1,row,col);
				energy_diff -=D*spin*spin;
				sitespin1[row][col]=newspin;
				energy_diff +=nn_energy(sitespin1,row,col);
				energy_diff +=D*newspin*newspin;


				if (inA(row))
				  { energy_diff -=nn_energy(sitespin2,row,col);
					energy_diff -=D*spin*spin;
					sitespin2[row][col]=newspin;
					energy_diff +=nn_energy(sitespin2,row,col);
					energy_diff +=D*newspin*newspin;
			      }

				//Generate a random no. r such that 0 < r < 1
				r = random_real(0, 1);
				acc_ratio = exp(-1.0 * energy_diff *beta);

				//Spin flipped if r <= acceptance ratio
				if (r <= acc_ratio)
				{

						energy += energy_diff;
				}
				else 
				{   sitespin1[row][col] =spin;
					if (inA(row))
						sitespin2[row][col]=spin;
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

        		spin = sitespin2[row][col];
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
 
				energy_diff =-nn_energy(sitespin2,row,col);
				energy_diff -=D*spin*spin;
				sitespin2[row][col]=newspin;
				energy_diff +=nn_energy(sitespin2,row,col);
				energy_diff +=D*newspin*newspin;


				if (inA(row))
				  { energy_diff -=nn_energy(sitespin1,row,col);
					energy_diff -=D*spin*spin;
					sitespin1[row][col]=newspin;
					energy_diff +=nn_energy(sitespin1,row,col);
					energy_diff +=D*newspin*newspin;
			      }


				//Generate a random no. r such that 0 < r < 1
				r = random_real(0, 1);
				acc_ratio = exp(-1.0 * energy_diff *beta);

				//Spin flipped if r <= acceptance ratio
				if (r <= acc_ratio)
				{

						energy += energy_diff;
				}
				else 
				{   sitespin2[row][col] =spin;
					if (inA(row))
						sitespin1[row][col]=spin;
				}
			}
		}

		if (i> 1e5) en_sum += energy;
		if (i> 1e5){
            ratio_sum += calc_ratio(sitespin1, sitespin2, beta);
            ratio_counter++;
        }
        if (ratio_counter >= 1e5){
            fratio << ratio_sum / ratio_counter << endl;
            ratio_counter = 0;
            ratio_sum = 0;
        }
	}

	//fout << beta << '\t' << en_sum / N_mc << endl;
	//fratio << ratio_sum / N_mc << endl;
    fratio.close();
	}

	//fout.close();
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
//with periodic boundary conditions

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

	//periodic boundary conditions
	for (unsigned int j = 0; j < axis2; ++j)
		{	
			energy +=-J * sitespin[axis1 - 1][j] * sitespin[0][j];
			energy +=D*sitespin[axis1-1][j]*sitespin[axis1-1][j];
		}

	for (unsigned int i = 0; i < axis1; ++i)
		{	
			energy +=-J * sitespin[i][axis2-1] * sitespin[i][0];
			energy +=D*sitespin[i][axis2-1]*sitespin[i][axis2-1];
		}

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
		nn_en -=J*sitespin[0][col]*sitespin[axis1-1][col];
		nn_en -=J*sitespin[0][col]*sitespin[1][col];

	}

	if (row == axis1 - 1)
	{
		nn_en -=J*sitespin[axis1-1][col]*sitespin[axis1-2][col];
		nn_en -=J*sitespin[axis1-1][col]*sitespin[0][col];

	}

	if (col == 0)
	{
		nn_en -=J*sitespin[row][0]*sitespin[row][axis2-1];
		nn_en -=J*sitespin[row][0]*sitespin[row][1];

	}

	if (col == axis2 - 1)
	{
		nn_en -=J*sitespin[row][axis2-1]*sitespin[row][axis2-2];
		nn_en -=J*sitespin[row][axis2-1]*sitespin[row][0];

	}
	return nn_en;
}

// Method for calcuiating the ratio of partition functions using the transfer matrix method
// regionSize contains the current number of spins in region A
long double calc_ratio(const array_2d& s1, const array_2d& s2, const double beta){
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
    
    // regionSize will be from 0 to L-1 (it could be L, but since we are always adding spins that isn't necessary)
    // We are doing the transfer matrix over the spins in row=regionSize
    int row_up = (regionSize+1)%axis1;
    int row_down = (regionSize-1+axis1)%axis1;

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
        tbot += log(t_bot.maxCoeff());
        t_bot /= t_bot.maxCoeff();

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
        ttop += log(t_top.maxCoeff());
        t_top /= t_top.maxCoeff();

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
        tcon += log(t_con.maxCoeff());
        t_con /= t_con.maxCoeff();
    }
    ttop += log(t_top.trace());
    tbot += log(t_bot.trace());
    tcon += log(t_con.trace());

    return exp(tcon - ttop - tbot);
}

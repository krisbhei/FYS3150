#include "mpi.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include "libising.h"
#include <stdlib.h>
#include <cmath>
#include <random>
#include <string>


using namespace std;
ofstream ofile;

void writeExpectedValuesTwoSpin(int dim, double temp, double expectations[5],long trials)
{
    //Start: Initialize computed values
    double normalizing = 1./((double)trials);
    double E = expectations[0]*normalizing;
    double E2 = expectations[1]*normalizing;

    double M = expectations[2]*normalizing;
    double M2 = expectations[3]*normalizing;
    double absM = expectations[4]*normalizing;

    normalizing = 1./dim/dim;

    double E_variance = (E2 - E*E)*normalizing;
    double M_variance = (M2 - absM*absM)*normalizing;

    double E_comp = E*normalizing;
    double heatCapacity_comp =  E_variance/temp/temp;
    double avgAbsM_comp = absM*normalizing;
    double X_comp = M_variance/temp;
    //End: Initialize computed values

    //Start: Analytical values
    double beta = 1.0/temp;

    double Z        = 4*(cosh(8*beta) + 3);
    double dividend = 1./(cosh(8*beta) + 3);

    E        = -8*sinh(8*beta)*dividend;
    E2       = 64*cosh(8*beta)*dividend;

    absM     = 2*(exp(8*beta) + 2)*dividend;
    M2       = 8*(exp(8*beta) + 1)*dividend;

    E_variance        = beta*beta*(64*(1+3*cosh(8*beta)))*dividend*dividend*normalizing;
    M_variance        = beta*(M2 - absM*absM)*normalizing;

    double E_analytical = E*normalizing;
    double heatCapacity_analytical =  E_variance/temp/temp;
    double avgAbsM_analytical = absM*normalizing;
    double X_analytical = M_variance/temp;
    //End: Analytical values

    ofile << " & Computed values & Analytical values \\\\ \\hline " << endl;
    ofile << "Value of the energy & " << E_comp << " & "<< E_analytical << "\\\\"<<endl;
    ofile << "Heat capacity & " << heatCapacity_comp << " & " <<heatCapacity_analytical <<"\\\\"<< endl;
    ofile << "Mean value of the magnetization & " << avgAbsM_comp<< " & " << avgAbsM_analytical<<"\\\\"<< endl;
    ofile << "Susceptibility & " << X_comp << " & " << X_analytical << endl;
}

void twoSpinTest()
{
    const int L = 2;
    const double T = 1.0;
    double trials[] = {1E12};//{1E2,1E3,1E4,1E5};//,1E6,1E7,1E8};
    for (double trial : trials)
    {
        int ** spins = init_matr(L);
        trial = (long) trial;
        double expectations[5];
        for (int i = 0 ; i < 5 ; i++)expectations[i] = 0;

        double energy = 0, magnetization = 0;
        initialize(L,spins,energy, magnetization);

        metropolis(spins,L,T,expectations,energy,magnetization,0,trial,-1);

        string filename = string("twoSpinLatex");
        stringstream ss;
        ss << setprecision(8) << trial;
        filename += string("_trials=")+ss.str();
        filename += string(".txt");

        ofile.open(filename);
        writeExpectedValuesTwoSpin(L,T,expectations,trial);
        ofile.close();
    }
    cout << "End twoSpin test" << endl;
    return;
}

void metropolisProbability(int ** spins,int dim, int trials, double T, double w[17], double energy, double magnetization,double * searchEnergies,double * hist,int i_max,string filename)
{
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_real_distribution<double> distr(0.0,1.0);

    double norm = 1./dim/dim;

    double normalizedEnergyNonRandom=0;
    double normalizedEnergyRandom=0;
    double normalizedMagnetizationNonRandom=0;
    double normalizedMagnetizationRandom=0;

    double E_Prev = 0;
    double absM_Prev = 0;
    double E2 = 0;

    double E = energy;
    double absM = magnetization;

    double tol = 1E-5;
    int numCyclesBeforeLikely = 1;

    //Start: Go toward the likely state
    while(((fabs((E - E_Prev))/numCyclesBeforeLikely)*norm > tol) || (((fabs((absM - absM_Prev)))/numCyclesBeforeLikely)*norm > tol))
    {
        metropolisOneCycle(dim,spins,energy,magnetization,w);
        E_Prev = E;
        absM_Prev = absM;

        E += energy;
        E2 += energy*energy;
        absM += fabs(magnetization);
        ++numCyclesBeforeLikely;
    }
    //End: Go toward the likely state

    cout << numCyclesBeforeLikely << endl;
    tol = 1E-5;

    double E_min = -2*dim*dim;

    //Start: counting the number of searchEnergy
    for(int cycles = numCyclesBeforeLikely+1 ; cycles <= trials ; cycles ++)
    {
        metropolisOneCycle(dim,spins,energy,magnetization,w);
        E += energy;
        double histIndex = (energy - E_min)/4.;
        hist[(int)histIndex] += 1;

        E2 += energy*energy;
        absM += fabs(magnetization);
    }
    //End: counting the number of searchEnergy

    //Start: write to file the computed values
    ofile.open(filename);
    ofile << T << endl;
    for(int i = 0 ; i < i_max ; i ++)ofile << setw(10) << (E_min + 4*i)*norm;
    ofile << endl;
    for(int i = 0 ; i < i_max ; i ++)
    {
        ofile << hist[i] << endl;
        hist[i] = 0;
    }

    E2 /= trials;
    E /= trials;

    double energySTD = sqrt(E2- E*E);
    cout << energySTD << endl;
    ofile << energySTD << endl;

    ofile.close();
    //End: Write to file
    return;
}

void probableEnergy()
{
    const int L = 20;
    int trials = 1E7;

    double e_min = -2*L*L;
    double e_max = 2*L*L;
    int i_max = (800.-e_min)/4.+1;
    double * searchEnergies = new double[i_max];

    for(int i = 0 ; i < i_max ; i ++)
    {
        searchEnergies[i] = e_min + 4*i;
    }
    double * histNonRandom = new double[i_max];
    double * histRandom = new double[i_max];
    for(int i = 0 ; i < i_max ; i ++)
    {
        histNonRandom[i] = 0;
        histRandom[i] = 0;
    }

    for (double T : {1.,2.4})
    {
                //cout << T << endl;

                double w[17];
                for(int i = 0; i < 17 ; i++) w[i] = 0;
                for(int i = -8; i < 9 ; i+=4) w[i+8] = exp(-((double)i)/T);

                int ** spinsNonRandom = init_matr(L);
                double energyNonRandom = 0, magnetizationNonRandom = 0;
                initialize(L,spinsNonRandom,energyNonRandom, magnetizationNonRandom);

                double energyRandom = 0, magnetizationRandom = 0;
                int ** spinsRandom = init_matr(L);
                initializeRandom(L,spinsRandom,energyRandom,magnetizationRandom);

                string filename = string("searchEnergies");
                filename += string("_temp=")+to_string(T);
                metropolisProbability(spinsNonRandom,L,trials,T,w,energyNonRandom,magnetizationNonRandom,searchEnergies,histNonRandom,i_max,filename+string("_ordered.dat"));
                metropolisProbability(spinsRandom,L,trials,T,w,energyRandom,magnetizationRandom,searchEnergies,histRandom,i_max,filename+string("_random.dat"));

    }

    return;
}

void metropolisLikelyState(int dim, int trials, double T)
{
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_real_distribution<double> distr(0.0,1.0);

    double w[17];
    for(int i = 0; i < 17 ; i++) w[i] = 0;
    for(int i = -8; i < 9 ; i+=4) w[i+8] = exp(-((double)i)/T);

    int ** spinsNonRandom = init_matr(dim);
    int ** spinsRandom = init_matr(dim);

    double energyNonRandom = 0, magnetizationNonRandom = 0;
    double energyRandom = 0, magnetizationRandom = 0;

    initialize(dim,spinsNonRandom,energyNonRandom, magnetizationNonRandom);
    initializeRandom(dim,spinsRandom,energyRandom,magnetizationRandom);

    int acceptedNonRandom = 0;
    int acceptedRandom = 0;

    //Initializing expectation values to calculate variance and plot
    double E_Random = 0;
    double E2_Random = 0;
    double absM_Random = 0;

    double E2_NonRandom = 0;
    double E_NonRandom = 0;
    double absM_NonRandom = 0;
    //End initializing

    double norm = 1./dim/dim;
    for(int cycle = 1 ; cycle <= trials ; cycle ++ )
    {
        for(int i = 0 ; i < dim*dim ; i++)
        {
            //Start: Spend a MC-cycle at a non random configuration
            int r_x = (int) (distr(gen)*(double)dim);
            int r_y = (int) (distr(gen)*(double)dim);

            int deltaEnergy = 2*spinsNonRandom[r_y][r_x]*
                    (spinsNonRandom[r_y][periodic(r_x,dim,-1)] +
                    spinsNonRandom[periodic(r_y,dim,-1)][r_x] +
                    spinsNonRandom[r_y][periodic(r_x,dim,1)] +
                    spinsNonRandom[periodic(r_y,dim,1)][r_x]);

            if( distr(gen) <= w[deltaEnergy + 8])
            {
                spinsNonRandom[r_y][r_x] *= -1;
                magnetizationNonRandom += (double) 2*spinsNonRandom[r_y][r_x];
                energyNonRandom += (double) deltaEnergy;

                ++acceptedNonRandom;
            }
            //End: Spend a MC-cycle at a non random configuration

            //Start: Spend a MC-cycle at a random configuration
            r_x = (int) (distr(gen)*(double)dim);
            r_y = (int) (distr(gen)*(double)dim);

            deltaEnergy = 2*spinsRandom[r_y][r_x]*
                    (spinsRandom[r_y][periodic(r_x,dim,-1)] +
                    spinsRandom[periodic(r_y,dim,-1)][r_x] +
                    spinsRandom[r_y][periodic(r_x,dim,1)] +
                    spinsRandom[periodic(r_y,dim,1)][r_x]);

            if( distr(gen) <= w[deltaEnergy + 8])
            {
                spinsRandom[r_y][r_x] *= -1;
                magnetizationRandom += (double) 2*spinsRandom[r_y][r_x];
                energyRandom += (double) deltaEnergy;
                ++acceptedRandom;
            }
            //End: Spend a MC-cycle at a random configuration
        }
        E_NonRandom += energyNonRandom;
        E2_NonRandom += energyNonRandom*energyNonRandom;
        absM_NonRandom += fabs(magnetizationNonRandom);

        E_Random += energyRandom;
        E2_Random += energyRandom*energyRandom;
        absM_Random += fabs(magnetizationRandom);

        if (cycle%100 == 0 )
        {
            ofile << (E_NonRandom/((double)cycle))*norm << setw(20) << (absM_NonRandom/((double)cycle))*norm<< setw(20) << cycle <<setw(20)<< (acceptedNonRandom/((double)cycle))*norm*100;
            ofile << endl;
            ofile << (E_Random/((double)cycle))*norm << setw(20) << (absM_Random/((double)cycle))*norm << setw(20) << cycle <<setw(20)<< (acceptedRandom/((double)cycle))*norm*100;
            ofile << endl;
        }
    }

    return;
}

void mostLikelyState()
{
    const int L = 20;
    int trials = 1E7;

    for (double T = 2.4; T <= 2.4 ; T +=1.4)
    {

                string filename = string("lmostLikelyState");
                stringstream ss;
                ss << setprecision(8) << trials;
                filename += string("_trials=")+ss.str();
                ss.str(string());
                ss << setprecision(8) << T;
                filename+=string("_temp=")+ss.str();
                filename += string(".dat");

                ofile.open(filename);
                ofile << T << setw(10) << trials << endl;
                metropolisLikelyState(L,trials,T);
                ofile.close();

    }
    return;
}

void writeExpectedValuesPhase(double T,int numSpins,int trialsPrProc,int num_processors,double totalExpectations[5])
{
    double normalizing = 1./trialsPrProc/num_processors;
    double inverse_temp = 1./T;
    double norm_numSpins = 1./numSpins/numSpins;

    double avgE = totalExpectations[0]*normalizing;
    double avgE2 = totalExpectations[1]*normalizing;
    double heatCapacity = (avgE2 - avgE*avgE)*norm_numSpins*inverse_temp*inverse_temp;

    double avgM2 = totalExpectations[3]*normalizing;
    double avgAbsM = totalExpectations[4]*normalizing;
    double susceptibility = (avgM2 - avgAbsM*avgAbsM)*norm_numSpins*inverse_temp;

    ofile << avgE*norm_numSpins << setw(10) << heatCapacity << setw(10);
    ofile << avgAbsM*norm_numSpins << setw(10) << susceptibility << endl;

}
void phaseTransitions()
{
    MPI_Init (NULL, NULL);

    int num_processors,this_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &num_processors);
    MPI_Comm_rank (MPI_COMM_WORLD, &this_rank);

    double start_T = 2.;
    double end_T = 2.4;
    double step_T = 0.015;

    int trials = 1E6;

    int trialsPrProc = (int)((double)trials/num_processors);
    int this_cycleStart = this_rank*trialsPrProc + 1;
    int this_cycleEnd = (this_rank+1)*trialsPrProc;
    if((this_rank == num_processors - 1) && (this_cycleEnd < trials)) this_cycleEnd = trials;

    if(this_rank == 0)
    {

        string filename = string("phaseTransitions");

        stringstream ss;
        ss << setprecision(4) << start_T;
        filename += string("_Tstart=") + ss.str();
        ss.str(string());
        ss << setprecision(4) << end_T;
        filename += string("_Tend=") + ss.str();
        ss.str(string());
        ss << setprecision(4) << step_T;
        filename += string("_Tstep=") + ss.str();

        filename += string(".dat");
        ofile.open(filename);
        for(double T = start_T ; T <= end_T ; T += step_T) ofile <<setw(10)<< T;
        ofile << endl;
    }

    for(int numSpins : {40,60,100,140})
    {

        for(double T = start_T ; T <= end_T ; T += step_T)
        {

            int ** spins = init_matr(numSpins);

            double expectations[5],totalExpectations[5];
            for(int i = 0 ; i < 5 ; i ++)
            {
                totalExpectations[i] = 0;
                expectations[i] = 0;
            }

            double energy = 0,magnetization = 0;
            initialize(numSpins,spins,energy,magnetization);

            metropolis(spins,numSpins,T,expectations,energy,magnetization,this_cycleStart,this_cycleEnd,this_rank);

            MPI_Reduce(&expectations,&totalExpectations,5,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
            if(this_rank == 0)
            {
                writeExpectedValuesPhase(T,numSpins,trialsPrProc,num_processors,totalExpectations);

            }
        }

    }
    if(this_rank == 0) ofile.close();

    MPI_Finalize();
    return;
}
void extractCriticalTemperature()
{
    MPI_Init (NULL, NULL);

    int num_processors,this_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &num_processors);
    MPI_Comm_rank (MPI_COMM_WORLD, &this_rank);


    int trials = 1E6;

    int trialsPrProc = (int)((double)trials/num_processors);
    int this_cycleStart = this_rank*trialsPrProc + 1;
    int this_cycleEnd = (this_rank+1)*trialsPrProc;
    if((this_rank == num_processors - 1) && (this_cycleEnd < trials)) this_cycleEnd = trials;
    double temp_start = 2.;
    double temp_end = 2.4;
    double steps[]= {0.025,0.001,0.025};



    for(int numSpins : {40,60,100,140})
    {
        double prevHeatC = 0;
        double critTemp = 0;
        if(this_rank == 0)
        {
            cout << numSpins << endl;
            string filename = string("criticalTemp");

            stringstream ss;
            ss << setprecision(4) << temp_start;
            filename += string("_Tstart=") + ss.str();
            ss.str(string());
            ss << setprecision(4) << temp_end;
            filename += string("_Tend=") + ss.str()+to_string(numSpins);
            ss.str(string());

            filename += string(".dat");
            ofile.open(filename);
            for(int i = 0 ; i < 3 ; i++)
            {
               for(double T = temp_start; T <= temp_end ; T += steps[i]) ofile <<setw(10)<< T;
            }

            ofile << endl;
            ofile << numSpins << endl;
        }

        double start_Ts[] = {2.1,2.2,2.4};
        double end_Ts[] = {2.2,2.4,2.5};
        for( int i = 0 ; i < 3 ; i++)
        {
            for(double T = start_Ts[i] ; T < end_Ts[i] ; T += steps[i])
            {
                int ** spins = init_matr(numSpins);

                double expectations[5],totalExpectations[5];
                for(int i = 0 ; i < 5 ; i ++)
                {
                    totalExpectations[i] = 0;
                    expectations[i] = 0;
                }
                double energy = 0, magnetization = 0;
                initialize(numSpins,spins,energy, magnetization);


               // void metropolis(int** spins,int dim, double T,double expectations[5],double energy, double magnetization,int cycleStart,int cycleEnd,int rank)
                metropolis(spins,numSpins,T,expectations,energy,magnetization,this_cycleStart,this_cycleEnd,this_rank);
                MPI_Reduce(&expectations,&totalExpectations,5,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
                if(this_rank == 0)
                {
                    writeExpectedValuesPhase(T,numSpins,trialsPrProc,num_processors,totalExpectations);
                    double normalizing = 1./trialsPrProc/num_processors;
                    double heatCapacity = ( totalExpectations[1]*normalizing -  totalExpectations[0]* totalExpectations[0]*normalizing*normalizing)/T/T;
                   // cout << T << endl;
                    if (heatCapacity > prevHeatC)
                    {
                        critTemp = T;
                        prevHeatC = heatCapacity;
                    }
                }
            }

        }

        if(this_rank == 0)
        {
            //cout << critTemp << endl;
            ofile << critTemp << endl;
            ofile.close();
            prevHeatC = 0;
            critTemp = 0;
        }

    }

    MPI_Finalize();
    return;
}

int main(int argc, char *argv[])
{
    twoSpinTest();
    //mostLikelyState();
    //probableEnergy();
    //phaseTransitions();
    //extractCriticalTemperature();
    return 0;
}

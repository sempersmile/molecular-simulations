#include <iostream>
#include <fstream>
#include <vector>
#include <time.h>
#include <random>

using namespace std;

/***************************************
* Global Variables
***************************************/
int N = 13;
int STEPS = 100;

/***************************************
* Function Declarations
***************************************/
vector<vector<double>> read_initial_positions(string inFile_name);
vector<vector<double>> initialize_random_velocities(int N);
vector<vector<double>> rescale_velocities(vector<vector<double>> vel, double T);
void write_final_positions(string outFile_name, vector<vector<double>> pos, int N);

double calculate_Ekin(); //TODO
double calculate_T(double E_kin);

/***************************************
* Main Program
***************************************/
int main()
{
    // Declare variables
    double Etot = 0; // Total energy
    double Epot = 0; // Potential energy
    double Ekin = 0; // Kinetic energy
    double T; // Temperature
    vector<vector<double>> pos(N, vector<double>(3)); // positions
    vector<vector<double>> vel(N, vector<double>(3)); // velocities
    string inFile_name = "initial_positions_13atoms.txt";
    string outFile_name = "final_positions_13atoms.xyz";
    //ifstream inFile="initial_positions_13atoms.txt"; // for reading in initial positions
    //ofstream outFile="final_positions_13atoms.xyz"; // for writing final positions

    cout << "Temperature/K: ";
    cin >> T;

    // ------------------------------------
    // Initialization
    // ------------------------------------

    // Read initial positions from file
    pos = read_initial_positions(inFile_name);
    // Initialize velocities randomly
    vel = initialize_random_velocities(N);
    vel = rescale_velocities(vel, T);

    // ------------------------------------
    // MD cycle
    // ------------------------------------
    for (int i=0; i<STEPS; i++)
    {
        
    }

    // Write final positions to file
    write_final_positions(outFile_name, pos, N);

    return 0;
}

/***************************************
* Function definitions
***************************************/

vector<vector<double>> read_initial_positions(string inFile_name)
{
    vector<vector<double>> pos(N, vector<double>(3));
    ifstream inFile;
    inFile.open(inFile_name);
    if (!inFile) 
    {
        cerr << "Unable to open file containing initial positions. Stopping program";
        exit(1);
    }
    int i = 0; // Iterator for atoms
    int j = 0; // Iterator for coordinates x,y,z
    double value;
    while(inFile >> value)
    {
        // If all coordinates of atom i have been saved --> Go to next atom
        if (j == 3)
        {
            i += 1;
            j -= 3;
        }
        if (value!=0)
        {
            pos[i][j] = value;
        }
        j += 1;
    }
    inFile.close();

    return pos;
}

vector<vector<double>> initialize_random_velocities(int N)
{
    random_device rd{};
    mt19937 generator{rd()};

    vector<vector<double>> vel(N, vector<double>(3));

    normal_distribution<double> distribution{0.0,5.0};

    for (int i=0; i<N; i++)
    {
        for (int j=0; j<3; j++)
        {
            vel[i][j]=distribution(generator);
        }
    }
    return vel;
}


vector<vector<double>> rescale_velocities(vector<vector<double>> vel, double T)
{
    
}


void write_final_positions(string outFile_name, vector<vector<double>> pos, int N)
{
    ofstream outFile;
    outFile.open(outFile_name);
    if (!outFile) 
    {
        cerr << "Unable to open file containing final positions. Stopping program";
        exit(1);
    }
    outFile << N << "\n\n";
    for (int i=0; i<N; i++)
    {
        outFile << "C"; // Atom type (choose 'C' here)
        for (int j=0; j<3; j++)
        {
            outFile << "\t" << pos[i][j];
        }
        outFile << "\n";
    }
    outFile.close();
}
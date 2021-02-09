#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <random>
#include <tuple>

using namespace std;

/***************************************
* Global Variables
***************************************/
int N = 13;
int NSTEPS = 1000;
double dt = 1;

/***************************************
* Function Declarations
***************************************/
vector<vector<double>> read_initial_positions(string inFile_name);
vector<vector<double>> initialize_random_velocities();
vector<vector<double>> rescale_velocities(vector<vector<double>> vel, double T);
tuple<double,vector<vector<double>>> calculate_forces(vector<vector<double>> pos);
void write_final_positions(string outFile_name, vector<vector<double>> pos);
double calculate_Ekin(vector<vector<double>> vel);
double calculate_T(double Ekin);


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
    vector<double> Etots(NSTEPS);
    vector<double> Epots(NSTEPS);
    vector<double> Ekins(NSTEPS);
    vector<double> Ts(NSTEPS);
    vector<vector<double>> pos(N, vector<double>(3)); // positions
    vector<vector<double>> vel(N, vector<double>(3)); // velocities
    vector<vector<double>> fcs(N, vector<double>(3)); // forces
    string inFile_name = "initial_positions_13atoms.txt";
    string outFile_name = "final_positions_13atoms.xyz";

    // Ask temperature from user
    cout << "Temperature (in reduced units): ";
    cin >> T;

    // ------------------------------------
    // Initialization
    // ------------------------------------

    // Read initial positions from file
    pos = read_initial_positions(inFile_name);
    // Initialize velocities randomly
    vel = initialize_random_velocities();
    // vel = rescale_velocities(vel, T); // Don't rescale velocities in canonical ensemble


    // ------------------------------------
    // MD cycle
    // ------------------------------------
    for (int i=0; i<NSTEPS; i++)
    {
        // Calculate kinetic energy and temperature
        Ekin = calculate_Ekin(vel);
        T = calculate_T(Ekin)

        // Calculate forces and potential energy
        tie(Epot, fcs) = calculate_forces(pos);

        // Bookkeeping of the relevant variables
        Etot = Ekin + Etot;
        Etots[i] = Etot;
        Ekins[i] = Ekin;
        Epots[i] = Epot;
        Ts[i] = T;


    }

    // Write final positions to file
    write_final_positions(outFile_name, pos);
    // Write final velocities to file
    write_final_positions("final_velocities_13atoms.txt", vel);

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

vector<vector<double>> initialize_random_velocities()
{
    random_device rd{};
    mt19937 generator{rd()};

    vector<vector<double>> vel(N, vector<double>(3));

    normal_distribution<double> distribution{0.0,1.0};

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
    //TODO Do we need to normalize velocity centre of mass to 0?
    double sum_vel2 = 0;
    for (int i=0; i<N; i++)
    {
        for (int j=0; j<3; j++)
        {
            sum_vel2 = sum_vel2 + vel[i][j]*vel[i][j];
        }
    }
    double scaling_factor = sqrt(3*N*T/sum_vel2);
    for (int i=0; i<N; i++)
    {
        for (int j=0; j<3; j++)
        {
            vel[i][j] = vel[i][j]*scaling_factor;
        }
    }
    return vel;
}


tuple<double,vector<vector<double>>> calculate_forces(vector<vector<double>> pos)
{
    vector<vector<double>> fcs(N, vector<double>(3));
    // Initialize forces to zero
    for (int i=0; i<N; i++)
    {
        for (int j=0; j<3; j++)
        {
            fcs[i][j]=0;
        }
    }
    // Loop over all atom pairs
    double E_pot = 0;
    for (int i=0; i<N; i++)
    {
        for (int k=i+1; k<N; k++)
        {
            double r = sqrt((pos[i][0] - pos[k][0])*(pos[i][0] - pos[k][0])
                    + (pos[i][1] - pos[k][1])*(pos[i][1] - pos[k][1])
                    + (pos[i][2] - pos[k][2])*(pos[i][2] - pos[k][2]));
            E_pot = E_pot + 4*(1/(pow(r,12)) - 1/(pow(r,6)));
            for (int j=0; j<3; j++)
            {
                double r_j = pos[i][j] - pos[k][j];
                double force_j = 48*r_j/(r*r) * ((1/pow(r,12)) - 0.5*1/(pow(r,6)));
                fcs[i][j] = fcs[i][j] + force_j;
                fcs[k][j] = fcs[k][j] - force_j;
            }
        }
    }

    return tuple<double,vector<vector<double>>>{E_pot, fcs};
}


void write_final_positions(string outFile_name, vector<vector<double>> pos)
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

double calculate_Ekin(vector<vector<double>> vel)
{
    double sumv2 = 0;
    for (int i=0; i<N; i++)
    {
        for (int j=0; j<3; j++)
        {
            sumv2 = sumv2 + vel[i][j]*vel[i][j];
        }
    }
    return 0.5*sumv2;
}

double calculate_T(double Ekin)
{
    return 2/3 * 1/N * Ekin;
}

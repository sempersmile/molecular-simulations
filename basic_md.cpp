#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

/***************************************
* Global variables
***************************************/
int N = 13;

/***************************************
* Function declarations
***************************************/
double calculate_Ekin();
double calculate_T(double E_kin);


int main()
{
    // Declare variables
    double Etot = 0; // Total energy
    double Epot = 0; // Potential energy
    double Ekin = 0; // Kinetic energy
    double T; // Temperature
    vector<vector<double>> pos(N, vector<double>(3)); // positions
    vector<vector<double>> vel(N, vector<double>(3)); // velocities
    ifstream inFile; // for reading in initial positions
    ofstream outFile; // for writing final positions

    // Read from file containing initial positions
    inFile.open("initial_positions_13atoms.txt");
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
            cout << i << j << endl;
            cout << pos[i][j] << endl;
        }
        j += 1;
    }
    inFile.close();

    // Write to file containing final positions
    outFile.open("final_positions_13atoms.xyz");
    if (!outFile) 
    {
        cerr << "Unable to open file containing final positions. Stopping program";
        exit(1);
    }
    int i = 0; // Iterator for atoms
    int j = 0; // Iterator for coordinates x,y,z
    double value;
    outFile << N << "\n\n";
    for (i=0; i<N; i++)
    {
        outFile << "C\t" // Atom type (choose 'C' here)
        for (j=0; j<3; j++)
        {
            // If all coordinates of atom i have been saved --> Go to next atom
            if (j == 4)
            {
                i += 1;
                j -= 3;
            }
            if (value!=0)
            {
                outFile << pos[i][j] << "\t";
            }
            j += 1;
        }
        outFile << "\n";
    }
    outFile.close();
    

    return 0;
}

/***************************************
* Function definitions
***************************************/
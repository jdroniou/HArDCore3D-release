#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <bits/stdc++.h>

int main(int argc, char * argv[])
{
     std::ifstream inFile;

     inFile.open("outputs/data_rates.dat");
     if (!inFile)
     {
          std::cerr << "Unable to open file\n";
          exit(1);
     }

    //  std::string header_line;
    //  std::getline(inFile, header_line);

     size_t num_objects = std::atoi(argv[1]);
     std::string item;
     std::vector<std::string> all_items;
     std::vector<double> each_items[num_objects];

     while (inFile >> item) 
     {    
          all_items.push_back(item);
     }

     inFile.close();

     size_t num_items = all_items.size();
     for (size_t i = num_objects; i < num_items; i++) {
          double temp = std::stod(all_items[i]);
          for (size_t j = 0; j < num_objects; j++) {
               if(i % num_objects == j) {
                    each_items[j].push_back(temp);
                    break;
               }
          } 
     }
    //  exit(1);
     num_items /= num_objects;
     num_items--;

     std::vector<double> meshsizes = each_items[0];
     std::vector<double> energyerrors = each_items[1];
    //  std::vector<double> cell_l2errors = each_items[2];
     std::vector<double> pressureerrors = each_items[2];
    //  std::vector<double> alt_h1errors = each_items[4];
    //  std::vector<double> jumperrors = each_items[5];
     std::vector<double> l2errors = each_items[3];

     std::string l2rate;
    //  std::string cell_l2rate;
     std::string pressurerate;
    //  std::string alt_h1rate;
    //  std::string jumprate;
     std::string energyrate;

     const int precision = 4;
     const int width = 12;

     std::cout.precision(precision);

     std::cout << "\n ---- Compute convergence rates ----\n";
     std::cout << "Pressure Error:\n";
     std::cout << std::setw(width) << "Mesh size";
     std::cout << std::setw(width) << "Error";
     std::cout << std::setw(width) << "Rate";
     std::cout << "\n";
     for (int i = 0; i < num_items; i++)
     {
          if (i >= 1)
          {
               l2rate = std::to_string((log(l2errors[i] / l2errors[i - 1])) / (log(meshsizes[i] / meshsizes[i - 1])));
          }
          std::cout << std::setw(width) << meshsizes[i];
          std::cout << std::setw(width) << l2errors[i];
          std::cout << std::setw(width) << l2rate;
          std::cout << "\n";
     }
    //  std::cout << "----------------------------------------\n";
    //  std::cout << "Cell L2 error:\n";
    //  std::cout << std::setw(width) << "Mesh size";
    //  std::cout << std::setw(width) << "Error";
    //  std::cout << std::setw(width) << "Rate";
    //  std::cout << "\n";
    //  for (int i = 0; i < num_items; i++)
    //  {
    //       if (i >= 1)
    //       {
    //            cell_l2rate = std::to_string((log(cell_l2errors[i] / cell_l2errors[i - 1])) / (log(meshsizes[i] / meshsizes[i - 1])));
    //       }
    //       std::cout << std::setw(width) << meshsizes[i];
    //       std::cout << std::setw(width) << cell_l2errors[i];
    //       std::cout << std::setw(width) << cell_l2rate;
    //       std::cout << "\n";
    //  }
     std::cout << "----------------------------------------\n";
     std::cout << "Magnetic Energy Error:\n";
     std::cout << std::setw(width) << "Mesh size";
     std::cout << std::setw(width) << "Error";
     std::cout << std::setw(width) << "Rate";
     std::cout << "\n";
     for (int i = 0; i < num_items; i++)
     {
          if (i >= 1)
          {
               pressurerate = std::to_string((log(pressureerrors[i] / pressureerrors[i - 1])) / (log(meshsizes[i] / meshsizes[i - 1])));
          }
          std::cout << std::setw(width) << meshsizes[i];
          std::cout << std::setw(width) << pressureerrors[i];
          std::cout << std::setw(width) << pressurerate;
          std::cout << "\n";
     }
    //  std::cout << "----------------------------------------\n";
    //  std::cout << "Alt H1 error:\n";
    //  std::cout << std::setw(width) << "Mesh size";
    //  std::cout << std::setw(width) << "Error";
    //  std::cout << std::setw(width) << "Rate";
    //  std::cout << "\n";
    //  for (int i = 0; i < num_items; i++)
    //  {
    //       if (i >= 1)
    //       {
    //            alt_h1rate = std::to_string((log(alt_h1errors[i] / alt_h1errors[i - 1])) / (log(meshsizes[i] / meshsizes[i - 1])));
    //       }
    //       std::cout << std::setw(width) << meshsizes[i];
    //       std::cout << std::setw(width) << alt_h1errors[i];
    //       std::cout << std::setw(width) << alt_h1rate;
    //       std::cout << "\n";
    //  }
    //  std::cout << "----------------------------------------\n";
    //  std::cout << "Jump error:\n";
    //  std::cout << std::setw(width) << "Mesh size";
    //  std::cout << std::setw(width) << "Error";
    //  std::cout << std::setw(width) << "Rate";
    //  std::cout << "\n";
    //  for (int i = 0; i < num_items; i++)
    //  {
    //       if (i >= 1)
    //       {
    //            jumprate = std::to_string((log(jumperrors[i] / jumperrors[i - 1])) / (log(meshsizes[i] / meshsizes[i - 1])));
    //       }
    //       std::cout << std::setw(width) << meshsizes[i];
    //       std::cout << std::setw(width) << jumperrors[i];
    //       std::cout << std::setw(width) << jumprate;
    //       std::cout << "\n";
    //  }
     std::cout << "----------------------------------------\n";
     std::cout << "Velocity Energy Error:\n";
     std::cout << std::setw(width) << "Mesh size";
     std::cout << std::setw(width) << "Error";
     std::cout << std::setw(width) << "Rate";
     std::cout << "\n";
     for (int i = 0; i < num_items; i++)
     {
          if (i >= 1)
          {
               energyrate = std::to_string((log(energyerrors[i] / energyerrors[i - 1])) / (log(meshsizes[i] / meshsizes[i - 1])));
          }
          std::cout << std::setw(width) << meshsizes[i];
          std::cout << std::setw(width) << energyerrors[i];
          std::cout << std::setw(width) << energyrate;
          std::cout << "\n";
     }

     return 0;
}

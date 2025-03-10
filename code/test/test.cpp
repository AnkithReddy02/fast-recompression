#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;

int main() {
    // Name of the file to be opened
    std::string filename = "einstein.en.txt";

    // Open the file in input mode
    std::ifstream file(filename);

    // Check if the file was opened successfully
    if (!file.is_open()) {
        std::cerr << "Error: Could not open the file " << filename << std::endl;
        return 1;
    }

    // Create a stringstream to hold the file contents
    std::stringstream buffer;

    // Read the file into the stringstream
    buffer << file.rdbuf();

    // Store the contents in a string
    std::string fileContents = buffer.str();

    // Close the file
    file.close();

    // Output the file contents
    // std::cout << "File contents:\n" << fileContents << std::endl;

    cout << fileContents[34] << ' ' << fileContents[35] << endl;

    return 0;
}

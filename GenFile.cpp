#include<bits/stdc++.h>
#include <unistd.h>
using namespace std;

class MutateChromosome {

private:
    string fasta_file;
    string chromosome;
    double mutation_rate;
    uint64_t targetFileSize;
    uint64_t copies;
    string OUTPUT_FILE = "mutated_genome.1GiB.0m_r";

    const vector<char> nucleotides_upper = {'A', 'T', 'C', 'G'};
    const vector<char> nucleotides_lower = {'a', 't', 'c', 'g'};
    mt19937 gen;
    
public:
    MutateChromosome(const string& fasta_file, double mutation_rate = 0.0001, uint64_t targetFileSize = (uint64_t)3 * 1024 * 1024 * 1024) : fasta_file(fasta_file), mutation_rate(mutation_rate), targetFileSize(targetFileSize), gen(random_device{}()) {
        OUTPUT_FILE = fasta_file + "." + to_string(targetFileSize/(1024 * 1024 * 1024)) + "GiB." + to_string(mutation_rate)+"mr";
	readChromosomeFromFASTA();
        cout << "Chromosome Size: " << chromosome.length() << endl;
        sleep(2);
        copies = targetFileSize / chromosome.length();
        // copies = 2;
        generateLargeMutatedFile();
    }

    string mutateSequence() {
        string mutated_seq = chromosome;

        random_device rd;
        (rd());
        uniform_real_distribution<> dis(0, 1);
        uniform_int_distribution<> nucleotide_dist_upper(0, 3);
        uniform_int_distribution<> nucleotide_dist_lower(0, 3);

        for(uint64_t i = 0; i < mutated_seq.size(); ++i) {
            if(dis(gen) < mutation_rate) {
		// cout << "Mutated" << endl;
                if(isupper(mutated_seq[i])) { 
                    char new_nucleotide = nucleotides_upper[nucleotide_dist_upper(gen)];

                    while(new_nucleotide == mutated_seq[i]) {
                        new_nucleotide = nucleotides_upper[nucleotide_dist_upper(gen)];
                    }
                    
                    mutated_seq[i] = new_nucleotide;
                } 
                else if(islower(mutated_seq[i])) {
                    char new_nucleotide = nucleotides_lower[nucleotide_dist_lower(gen)];

                    while(new_nucleotide == mutated_seq[i]) {
                        new_nucleotide = nucleotides_lower[nucleotide_dist_lower(gen)];
                    }
                    
                    mutated_seq[i] = new_nucleotide;
                }
            }
        }

        return mutated_seq;
    }

    void generateLargeMutatedFile() {
        ofstream file(OUTPUT_FILE);
        if(!file.is_open()) {
            cerr << "Error: Could not open output file." << endl;
            return;
        }

        for(uint64_t i = 0; i < copies; ++i) {
            cout << "COPY: " << i + 1 << endl;
            string mutated_seq = mutateSequence();
            file << mutated_seq;
        }
        file.close();
        cout << "Mutated genome file generated successfully." << endl;
    }

    void readChromosomeFromFASTA() {
        ifstream file(fasta_file);
        if(!file.is_open()) {
            cerr << "Error: Could not open FASTA file." << endl;
            exit(1);
        }

        string line;
        bool end = false;

        while(getline(file, line)) {
            if(line[0] != '>') {
                chromosome += line;
            }
            else {
                if(end) {
                    break;
                }

                end = true;
            }
        }

        file.close();
        return;
    }
};

int main(int argc, char** argv) {
    if(argc < 1) {
	cerr << "./GenFile input_file";
        exit(1);
    }
    string input_file(argv[1]);
    MutateChromosome mc(input_file);

}


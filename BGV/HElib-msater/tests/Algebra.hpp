#include <helib/NumbTh.h>
#include <helib/sample.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <set>
#include <unordered_map>
#include <cmath>
#include <iomanip>
#include <sys/stat.h> // for stat()
#include <filesystem>
#include <cassert>
#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace helib;
using namespace std;
using namespace NTL;

long getCoord(long i, long d, const Vec<long> prodvec){
    return (i % prodvec[d]) / prodvec[d+1];
}

long L(long a, long b){
    return (a % b) - b;
}

long basePrime(long n){
    Vec<Pair<long, long>> factors;
    factorize(factors, n);
    return factors[0].a;
}

void Matrix_secret(Mat<Vec<long>>& result, const Vec<long> mvec, const Vec<long> mhatvec, const Vec<long> m_prodvec, const Vec<long> phim_prodvec){
    
    long m = m_prodvec[0];
    long phim = phim_prodvec[0];
    long k = mvec.length();
    long pow = 1L << k;

    cout << "Generate thin rotation matrix of secret...";

    for (long i=0; i < phim; i++){
        for (long j = 0; j < phim; j++){   
            Vec<long> acc;
            for (long c=0; c < pow; c++){

                long ind = 0;
                long tau_sum = 0;

                for (long d = 0; d < k; d++){
                    long ind_d; 
                    long tau_d = bit(c, d);

                    long i_d = getCoord(i, d, phim_prodvec);
                    long j_d = getCoord(j, d, phim_prodvec);

                    if (tau_d == 0) ind_d = ((i_d - j_d) % mvec[d] + mvec[d]) % mvec[d];
                    else{
                        long p_d = basePrime(mvec[d]);
                        ind_d = ((L(i_d, mvec[d]/p_d) - j_d) % mvec[d] + mvec[d]) % mvec[d];
                    }

                    ind += ind_d * mhatvec[d];
                    tau_sum += tau_d;

                }

                long index = ind % m;
                
                if ((tau_sum % 2) == 0) append(acc, index+1);
                else append(acc, - index -1);

            }
            result[i][j] = acc;
        }
    }

    cout << " done (rows=" << phim << ", cols=" << phim << ")" << endl;
}

void Matrix_secret_thin(Mat<Vec<long>>& result, const Vec<long> mvec, const Vec<long> mhatvec, const Vec<long> m_prodvec, const Vec<long> phim_prodvec){
    
    long m = m_prodvec[0];
    long phim = phim_prodvec[0];
    long k = mvec.length();
    long pow = 1L << k;

    long nslots = result.NumRows();
    long degree = phim/nslots;

    cout << "Generate rotation matrix of secret index...";
    for (long i=0; i < nslots; i++){

        for (long j = 0; j < phim; j++){   
            Vec<long> acc;
            for (long c=0; c < pow; c++){

                long ind = 0;
                long tau_sum = 0;

                for (long d = 0; d < k; d++){
                    long ind_d; 
                    long tau_d = bit(c, d);

                    long i_d = getCoord(degree *i, d, phim_prodvec);
                    long j_d = getCoord(j, d, phim_prodvec);

                    if (tau_d == 0) ind_d = ((i_d - j_d) % mvec[d] + mvec[d]) % mvec[d];
                    else{
                        long p_d = basePrime(mvec[d]);
                        ind_d = ((L(i_d, mvec[d]/p_d) - j_d) % mvec[d] + mvec[d]) % mvec[d];
                    }

                    ind += ind_d * mhatvec[d];
                    tau_sum += tau_d;

                }

                long index = ind % m;
                
                if ((tau_sum % 2) == 0) append(acc, index+1);
                else append(acc, - index -1);

            }
            result[i][j] = acc;
        }

    }

    cout << " done (slots=" << nslots << ", phim=" << phim << ")" << endl;
}

void Matrix_ctxt(Mat<Vec<long>>& result, Mat<long>& matrix_irwin_hall, const Mat<Vec<long>>& secret_perm_matrix)
{
    long phim = result.NumRows();
    long m = result.NumCols();
    
    cout << "Rearrange with regard to ctxt coefficient index...";
    for (long i=0; i<phim; ++i){
        for (long ind = 0; ind < m; ++ind){
            Vec<long> acc;
            //find ind+1 or -ind-1 among secret_perm_matrix[i][j] (0 <=j < phim)
            for (long j=0; j < phim; ++j){
                //insert j+1 or -j-1 
                for (long k=0; k < secret_perm_matrix[i][j].length(); ++k){
                    if (secret_perm_matrix[i][j][k] == ind+1) append(acc, j+1);
                    else if (secret_perm_matrix[i][j][k] == -ind-1) append(acc, -j-1);
                }
            }
            result[i][ind] = acc;
            matrix_irwin_hall[i][ind] = acc.length();
        }
    }

    cout << " done (rows=" << phim << ", cols=" << m << ")" << endl;
}

void Matrix_ctxt_thin(Mat<Vec<long>>& result, Mat<long>& matrix_irwin_hall, const Mat<Vec<long>>& secret_perm_matrix_thin)
{
    long nslots = result.NumRows();
    long m = result.NumCols();
    long phim = phi_N(m);
    
    cout << "Arrange with regard to ctxt coefficient index...";
    for (long i=0; i < nslots; ++i){
        for (long ind = 0; ind < m; ++ind){
            Vec<long> acc;
            //find ind+1 or -ind-1 among secret_perm_matrix[i][j] (0 <=j < phim)
            for (long j=0; j < phim; ++j){
                //insert j+1 or -j-1 
                for (long k=0; k < secret_perm_matrix_thin[i][j].length(); ++k){
                    if (secret_perm_matrix_thin[i][j][k] == ind+1) append(acc, j+1);
                    else if (secret_perm_matrix_thin[i][j][k] == -ind-1) append(acc, -j-1);
                }
            }
            result[i][ind] = acc;
            matrix_irwin_hall[i][ind] = acc.length();
        }
    }

    cout << " done (rows=" << nslots << ", cols=" << m << ")" << endl;
}

template<typename T>
void getCoeff(Vec<T>& result, long i, const Vec<T>& ctxt_pwfl_coeff, const Mat<Vec<long>>& ctxt_perm_matrix)
{
    long m = result.length();
    long phim = phi_N(m);
    
    #pragma omp parallel for
	for (long j = 0; j < m; ++j){
		T acc = T(0);
		
		for (long k=0; k < ctxt_perm_matrix[i][j].length(); ++k){
			long index = ctxt_perm_matrix[i][j][k];

			if (index > 0) acc = acc + ctxt_pwfl_coeff[index-1];
			else if (index <0) acc = acc - ctxt_pwfl_coeff[-index-1];
			
		}

		result[j] = acc;
	}
}

template<typename T>
void getCoeff_thin(Vec<T>& result, long i, const Vec<T>& ctxt_pwfl_coeff, const Mat<Vec<long>>& ctxt_perm_matrix_thin)
{
    long m = result.length();
    long phim = phi_N(m);

    // Each j reads a distinct row of the permutation matrix and writes a
    // distinct element of result, so the iterations are independent.
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
	for (long j = 0; j < m; ++j){
		T acc = T(0);
		
		for (long k=0; k < ctxt_perm_matrix_thin[i][j].length(); ++k){
			long index = ctxt_perm_matrix_thin[i][j][k];

			if (index > 0) acc = acc + ctxt_pwfl_coeff[index-1];
			else if (index <0) acc = acc - ctxt_pwfl_coeff[-index-1];
			
		}

		result[j] = acc;
	}
}

template<typename T>
void printVec(const Vec<T>& list) {
    cout << "[ ";
    for (long i = 0; i < list.length(); i++) {
        cout <<  list[i] << " ";
    }
    cout << "]" << endl;
}

template<typename T>
void printMat(const Mat<T>& matrix) {
    cout << "[\n";
    for (long i = 0; i < matrix.NumRows(); i++) {
        printVec(matrix[i]);
    }
    cout << "]" << endl;
}

bool fileExists(const std::string& filename) {
    struct stat buffer;
    return (stat(filename.c_str(), &buffer) == 0);
}

template<typename T>
void loadData(const std::string& load_dir, std::vector<Vec<T>>& data)
{
    ifstream infile(load_dir);
    string line;

    while (getline(infile, line)) {
        stringstream ss(line);
        Vec<T> row;
        row.SetLength(0);

        T val;
        while (ss >> val) {
            row.append(val);
        }

        data.push_back(row); // store per row
    }

    infile.close();
}

template<typename T>
void saveVec(const Vec<T>& vec, const string& folder_path) 
{
    // store continuedly (append)
    std::ofstream outFile(folder_path, std::ios::app);  
    outFile << setprecision(10) << fixed;
    if (!outFile.is_open()) {
        cerr << "Failed to open file: " << folder_path << endl;
        return;
    }

    for (long i = 0; i < vec.length(); ++i) {
        outFile << vec[i];
        if (i != vec.length() - 1)
            outFile << " ";
    }
    outFile << "\n";
    outFile.close();
}

template<typename T>
void saveMat(const Mat<T>& mat, const string& folder_path) 
{
    ofstream outFile(folder_path);
    outFile << setprecision(10) << fixed;
    if (!outFile.is_open()) {
        cerr << "Failed to open file: " << folder_path << endl;
        return;
    }

    for (long i = 0; i < mat.NumRows(); ++i) {
        for (long j = 0; j < mat.NumCols(); ++j) {
            outFile << mat[i][j];
            if (j != mat.NumCols() - 1)
                outFile << " ";
        }
        outFile << "\n";
    }

    outFile.close();
}

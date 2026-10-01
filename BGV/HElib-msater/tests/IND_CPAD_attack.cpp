#include <NTL/BasicThreadPool.h>
#include <cassert>
#include <helib/helib.h>
#include <helib/matmul.h>
#include <helib/debugging.h>
#include <helib/ArgMap.h>
#include <helib/sample.h>
#include <filesystem>
#include <sys/stat.h>
#include <sys/types.h>
#include <fstream>
#include <iostream> 
#include <chrono>
#include <ctime>
#include <iomanip>
#include <sstream>
#include "Algebra.hpp"
#if defined(_OPENMP)
#include <omp.h>
#endif

NTL_CLIENT
using namespace helib;
using namespace std;
using namespace NTL;

static bool noPrint = false;
static bool dry = false; // a dry-run flag
static bool debug = 0;   // a debug flag
static bool waitForInput = false;
static bool progress = true;
static long progressStep = 10;
static bool stopAfterFirstFailure = false;
static long stopAfterFailCount = 0;
// static int scale = 10;
static int scale = 8;

class TeeBuf : public std::streambuf
{
public:
	TeeBuf(std::streambuf* lhs, std::streambuf* rhs) : lhs(lhs), rhs(rhs) {}

protected:
	int overflow(int c) override
	{
		if (c == traits_type::eof())
			return traits_type::not_eof(c);

		const int left = lhs ? lhs->sputc(static_cast<char>(c)) : c;
		const int right = rhs ? rhs->sputc(static_cast<char>(c)) : c;
		if (left == traits_type::eof() || right == traits_type::eof())
			return traits_type::eof();
		return c;
	}

	int sync() override
	{
		const int left = lhs ? lhs->pubsync() : 0;
		const int right = rhs ? rhs->pubsync() : 0;
		return (left == 0 && right == 0) ? 0 : -1;
	}

private:
	std::streambuf* lhs;
	std::streambuf* rhs;
};

static long mValues[][14] = {
//   { p, phi(m), m, d, m1,  m2, m3,   g1,    g2,    g3,ord1,ord2,ord3, c_m}
//   {  2,    48,   105, 12,  3,  35,  0,    71,    76,    0,  2,  2,   0, 200}, // prime-power-factorization (x)
//   {  2,   600,  1023, 10, 11,  93,  0,   838,   584,    0, 10,  6,   0, 100}, // m=(3)*11*{31} m/phim(m)=1.7    C=24  D=2 E=1, prime-power-factorization (x)
//   {  2,  1200,  1705, 20, 11, 155,  0,   156,   936,    0, 10,  6,   0, 100}, // m=(5)*11*{31} m/phim(m)=1.42   C=34  D=2 E=2, prime-power-factorization (x)
//   {  2,  1728,  4095, 12,  7,  5, 117,  2341,  3277, 3641,  6,  4,   6, 100}, // m=(3^2)*5*7*{13} m/phim(m)=2.36 C=26 D=3 E=2, prime-power-factorization (x)
//   {  2,  2304,  4641, 24,  7,  3, 221,  3979,  3095, 3760,  6,  2,  -8, 300}, // m=3*7*(13)*{17} :-( m/phim(m)=2.01 C=45 D=4 E=3, prime-power-factorization (x)
//   {  2,  4096,  4369, 16, 17, 257,  0,   258,  4115,    0, 16, -16,   0, 100}, // m=17*(257) :-( m/phim(m)=1.06 C=61 D=3 E=4, prime-power-factorization (o)
//   {  2, 12800, 17425, 40, 41, 425,  0,  5951,  8078,    0, 40, -8,   0, 100}, // m=(5^2)*{17}*41 m/phim(m)=1.36 C=93  D=3 E=3, prime-power-factorization (x)
//   {  2, 15004, 15709, 22, 23, 683,  0,  4099, 13663,    0, 22, 31,   0, 100}, // m=23*(683) m/phim(m)=1.04      C=73  D=2 E=1, prime-power-factorization (o)
// Parameter A*
     {  2, 16384, 21845, 16, 17,   5,257,  8996, 17477, 21591, 16, 4, -16, 200}, // m=5*17*(257) :-( m/phim(m)=1.33 C=65 D=4 E=4, prime-power-factorization (o)
//   {  2, 18000, 18631, 25, 31, 601,  0, 15627,  1334,    0, 30, 24,   0, 100}, // m=31*(601) m/phim(m)=1.03      C=77  D=2 E=0, prime-power-factorization (o)
//   {  2, 18816, 24295, 28, 43, 565,  0, 16386, 16427,    0, 42, 16,   0, 100}, // m=(5)*43*{113} m/phim(m)=1.29  C=84  D=2 E=2, prime-power-factorization (x)
//   {  2, 21168, 27305, 28, 43, 635,  0, 10796, 26059,    0, 42, 18,   0, 100}, // m=(5)*43*{127} m/phim(m)=1.28  C=86  D=2 E=2, prime-power-factorization (x)
// Parameter B*. Both sets are live; the N argument picks between them, since the
// search below takes the first row whose phi(m) is at least N. N=0 or 16384 gives
// A*, N=16385 up to 23040 gives B*.
     {  2, 23040, 28679, 24, 17,  7, 241, 15184,  4098,28204, 16,  6, -10, 200}, // m=7*17*(241) m/phim(m)=1.24    C=63  D=4 E=3, prime-power-factorization (o)
//   {  2, 24000, 31775, 20, 41, 775,  0,  6976, 24806,    0, 40, 30,   0, 100}, // m=(5^2)*{31}*41 m/phim(m)=1.32 C=88  D=2 E=2, prime-power-factorization (x)
//   {  2, 26400, 27311, 55, 31, 881,  0, 21145,  1830,    0, 30, 16,   0, 100}, // m=31*(881) m/phim(m)=1.03      C=99  D=2 E=0, prime-power-factorization (o)
//   {  2, 27000, 32767, 15, 31,  7, 151, 11628, 28087,25824, 30,  6, -10, 200}, // prime-power-factorization (o) 
//   {  2, 31104, 35113, 36, 37, 949,  0, 16134,  8548,    0, 36, 24,   0, 200}, // m=(13)*37*{73} m/phim(m)=1.12  C=94  D=2 E=2, prime-power-factorization (x)
//   {  2, 31104, 35113, 36, 73, 13, 37, 32709, 10805, 0, 72, 12, 0, 200}, // Custom, m=(13)*37*{73} m/phim(m)=1.12  C=94  D=2 E=2, prime-power-factorization (o)
//   {  2, 34848, 45655, 44, 23, 1985, 0, 33746, 27831,    0, 22, 36,   0, 100}, // m=(5)*23*{397} m/phim(m)=1.31  C=100 D=2 E=2, prime-power-factorization (x)
//   {  2, 42336, 42799, 21, 127, 337, 0, 25276, 40133,    0,126, 16,   0, 200}, // m=127*(337) m/phim(m)=1.01     C=161 D=2 E=0, prime-power-factorization (o)
//   {  2, 45360, 46063, 45, 73, 631,  0, 35337, 20222,    0, 72, 14,   0, 100}, // m=73*(631) m/phim(m)=1.01      C=129 D=2 E=0, prime-power-factorization (o)
//   {  2, 46080, 53261, 24, 17, 13, 241, 43863, 28680,15913, 16, 12, -10, 100}, // m=13*17*(241) m/phim(m)=1.15   C=69  D=4 E=3, prime-power-factorization (o)
//   {  2, 49500, 49981, 30, 151, 331, 0,  6952, 28540,    0,150, 11,   0, 100}, // m=151*(331) m/phim(m)=1        C=189 D=2 E=1, prime-power-factorization (o)
//   {  2, 54000, 55831, 25, 31, 1801, 0, 19812, 50593,    0, 30, 72,   0, 100}, // m=31*(1801) m/phim(m)=1.03     C=125 D=2 E=0, prime-power-factorization (o)
//   {  2, 60016, 60787, 22, 89, 683,  0,  2050, 58741,    0, 88, 31,   0, 200}, // m=89*(683) m/phim(m)=1.01      C=139 D=2 E=1, prime-power-factorization (o)
//   {  7,    36,    57,  3,  3,  19,  0,    20,    40,    0,  2, -6,   0, 100}, // m=3*(19) :-( m/phim(m)=1.58 C=14 D=3 E=0, prime-power-factorization (o)
//   { 17,    48,   105, 12,  3,  35,  0,    71,    76,    0,  2,  2,   0, 100}, // m=3*(5)*{7} m/phim(m)=2.18 C=14 D=2 E=2, prime-power-factorization (x)
//   { 17,   576,  1365, 12,  7,   3, 65,   976,   911,  463,  6,  2,   4, 100}, // m=3*(5)*7*{13} m/phim(m)=2.36  C=22  D=3, prime-power-factorization (x)
//   { 17, 18000, 21917, 30, 101, 217, 0,  5860,  5455,    0, 100, 6,   0, 100}, // m=(7)*{31}*101 m/phim(m)=1.21  C=134 D=2, prime-power-factorization (x)
//   { 17, 30000, 34441, 30, 101, 341, 0,  2729, 31715,    0, 100, 10,  0, 100}, // m=(11)*{31}*101 m/phim(m)=1.14 C=138 D=2, prime-power-factorization (x)
//   { 17, 40000, 45551, 40, 101, 451, 0, 19394,  7677,    0, 100, 10,  0, 200}, // m=(11)*{41}*101 m/phim(m)=1.13 C=148 D=2, prime-power-factorization (x)
//	 { 17, 40000, 45551, 40, 101, 11, 41, 19394,  37270,    0, 100, 10,  0, 200}, // CUSTOM m=(11)*{41}*101 m/phim(m)=1.13 C=148 D=2, prime-power-factorization (o)
//   { 17, 46656, 52429, 36, 109, 481, 0, 46658,  5778,    0, 108, 12,  0, 100}, // m=(13)*{37}*109 m/phim(m)=1.12 C=154 D=2, prime-power-factorization (x)
//   { 17, 54208, 59363, 44, 23, 2581, 0, 25811,  5199,    0, 22, 56,   0, 100}, // m=23*(29)*{89} m/phim(m)=1.09  C=120 D=2, prime-power-factorization (x)
//   { 17, 70000, 78881, 10, 101, 781, 0, 67167, 58581,    0, 100, 70,  0, 100}, // m=(11)*{71}*101 m/phim(m)=1.12 C=178 D=2, prime-power-factorization (x)
//   {127,   576,  1365, 12,  7,   3, 65,   976,   911,  463,  6,  2,   4, 100}, // m=3*(5)*7*{13} m/phim(m)=2.36   C=22  D=3, prime-power-factorization (x)
//   {127,  1200,  1925, 20,  11, 175, 0,  1751,   199,    0, 10, 6,    0, 100}, //  m=(5^2)*{7}*11 m/phim(m)=1.6   C=34 D=2, prime-power-factorization (x)
//   {127,  2160,  2821, 30,  13, 217, 0,   652,   222,    0, 12, 6,    0, 100}, // m=(7)*13*{31} m/phim(m)=1.3     C=46 D=2, prime-power-factorization (x)
//   {127, 18816, 24295, 28, 43, 565,  0, 16386, 16427,    0, 42, 16,   0, 100}, // m=(5)*43*{113} m/phim(m)=1.29   C=84  D=2, prime-power-factorization (x)
//   {127, 26112, 30277, 24, 17, 1781, 0, 14249, 10694,    0, 16, 68,   0, 100}, // m=(13)*17*{137} m/phim(m)=1.15  C=106 D=2, prime-power-factorization (x)
//	 {127, 31752, 32551, 14, 43,  757, 0,  7571, 28768,    0, 42, 54, 0, 100}, // m=43*(757) :-( m/phim(m)=1.02   C=161 D=3, prime-power-factorization (o), next candidate
//   {127, 46656, 51319, 36, 37, 1387, 0, 48546, 24976,    0, 36, -36,  0, 200}, //m=(19)*37*{73}:-( m/phim(m)=1.09 C=141 D=3, prime-power-factorization (x)
//	 {127, 46656, 51319, 36, 37, 19, 73, 48546, 35114, 13358, 36, 18,  2, 200}, // Custom, m=(19)*37*{73}:-( m/phim(m)=1.09 C=141 D=3, prime-power-factorization (o)
//   {127, 49392, 61103, 28, 43, 1421, 0,  1422, 14234,    0, 42, 42,   0, 200}, // m=(7^2)*{29}*43 m/phim(m)=1.23  C=110 D=2, prime-power-factorization (x)
//   {127, 54400, 61787, 40, 41, 1507, 0, 30141, 46782,    0, 40, 34,   0, 100}, // m=(11)*41*{137} m/phim(m)=1.13  C=112 D=2, prime-power-factorization (x)
//   {127, 72000, 77531, 30, 61, 1271, 0,  7627, 34344,    0, 60, 40,   0, 100}  // m=(31)*{41}*61 m/phim(m)=1.07   C=128 D=2, prime-power-factorization (x)
};

#define num_mValues (sizeof(mValues)/(14*sizeof(long)))

void TestIt(long idx, long p, long r, long L, long c, long skHwt, long skHwt_sp, int iter, int build_cache=0)
{
	helib::setQuietMode(noPrint);
	const auto attackStart = std::chrono::steady_clock::now();
	// Where the time goes, so that any attempt at threading knows what to aim at.
	double t_recrypt = 0.0, t_getcoeff = 0.0;
	auto elapsedString = [&attackStart]() {
		const auto now = std::chrono::steady_clock::now();
		const auto secs = std::chrono::duration_cast<std::chrono::seconds>(now - attackStart).count();
		const long h = static_cast<long>(secs / 3600);
		const long m = static_cast<long>((secs % 3600) / 60);
		const long s = static_cast<long>(secs % 60);
		return std::to_string(h) + "h " + std::to_string(m) + "m " + std::to_string(s) + "s";
	};

	Vec<long> mvec;
	vector<long> gens;
	vector<long> ords;

	long phim = mValues[idx][1];
	long m = mValues[idx][2];

	assert(GCD(p, m) == 1);

	append(mvec, mValues[idx][4]);
	if (mValues[idx][5]>1) append(mvec, mValues[idx][5]);
	if (mValues[idx][6]>1) append(mvec, mValues[idx][6]);
	gens.push_back(mValues[idx][7]);
	if (mValues[idx][8]>1) gens.push_back(mValues[idx][8]);
	if (mValues[idx][9]>1) gens.push_back(mValues[idx][9]);
	ords.push_back(mValues[idx][10]);
	if (abs(mValues[idx][11])>1) ords.push_back(mValues[idx][11]);
	if (abs(mValues[idx][12])>1) ords.push_back(mValues[idx][12]);

    long k = mvec.length();

    Vec<long> phimvec;
    phimvec.SetLength(k);
    for (int i=0; i<k; ++i){
        phimvec[i] = phi_N(mvec[i]);
    } 
    
    Vec<long> mhatvec, m_prodvec, phim_prodvec;
    mhatvec.SetLength(k);
    m_prodvec.SetLength(k+1);
    phim_prodvec.SetLength(k+1);
    
    for (long w = 0; w < k; w++) {
        long prod = 1;
        long phiprod = 1;
        for (long i=w; i<k; i++){
            prod *= mvec[i];
            phiprod *= phimvec[i];
        }
        m_prodvec[w] = prod;
        phim_prodvec[w] = phiprod;

        mhatvec[w] = m/mvec[w];
    }

    m_prodvec[k] = 1;
    phim_prodvec[k] = 1;

	if (!noPrint) {
	cout << "*** TestIt";
	if (isDryRun()) cout << " (dry run)";
	cout << ": p=" << p
	 << ", r=" << r
	 << ", L=" << L
	 << ", t=" << skHwt
	 << ", c=" << c
	 << ", m=" << m
	 << " (=" << mvec << "), gens="<<gens<<", ords="<<ords
	 << endl;
	cout << "Computing key-independent tables..." << std::flush;
	}
	
	setTimersOn();
	setDryRun(false); // Need to get a "real context" to test bootstrapping

	double t = -GetTime();

	Context context = helib::ContextBuilder<helib::BGV>()
				   .m(m)
				   .p(p)
				   .r(r)
				   .gens(gens)
				   .ords(ords)
				   .scale(scale ? scale : 10 /*10 is default*/)
				   .buildModChain(false)
				   .build();

	context.buildModChain(L,
						c,
						/*willBeBootstrappable=*/true,
						/*t=*/skHwt);

	if (!noPrint) {
	  std::cout << "security=" << context.securityLevel() << std::endl;
	  std::cout << "# small primes = " << context.getSmallPrimes().card()
				<< "\n";
	  std::cout << "# ctxt primes = " << context.getCtxtPrimes().card() << "\n";
	  std::cout << "# bits in ctxt primes = "
				<< long(context.logOfProduct(context.getCtxtPrimes()) /
							log(2.0) +
						0.5)
				<< "\n";
	  std::cout << "# special primes = " << context.getSpecialPrimes().card()
				<< "\n";
	  std::cout << "# bits in special primes = "
				<< long(context.logOfProduct(context.getSpecialPrimes()) /
							log(2.0) +
						0.5)
				<< "\n";
	  std::cout << "# prime factors of m = " << context.getNFactors() << std::endl;
	  std::cout << "scale=" << context.getScale() << std::endl;
	}

	context.enableBootStrapping(mvec, build_cache, /*alsoThick=*/false);
	// save time...disable some fat boot precomputation
	t += GetTime();

	//if (skHwt>0) context.rcData.skHwt = skHwt;
	if (!noPrint) {
        std::cout << " done in " << t << " seconds\n";
        std::cout << "  e=" << context.getRcData().e
                << ", e'=" << context.getRcData().ePrime
                << ", t=" << context.getRcData().skHwt
                << ", t'=" << skHwt_sp << "\n  ";
                
        context.printout();
	}

	setDryRun(dry); // Now we can set the dry-run flag if desired

	long p2r = context.getAlMod().getPPowR();

	string folder = "m" + std::to_string(m) + "_p" + std::to_string(p) + "_e" + std::to_string(context.getRcData().e) 
				+ "_ep" + std::to_string(context.getRcData().ePrime) + "_r" + std::to_string(r);

	if (!noPrint)
		cerr << "folder name = " <<  folder << std::endl;

	const std::filesystem::path output_root = std::filesystem::path(".") / "Collected_data";
	string path_dir = (output_root / folder).string();

	if (!std::filesystem::exists(path_dir)) {
		std::filesystem::create_directories(path_dir);
	}

	const auto now = std::chrono::system_clock::now();
	const std::time_t now_tt = std::chrono::system_clock::to_time_t(now);
	std::tm now_tm{};
#if defined(_WIN32)
	localtime_s(&now_tm, &now_tt);
#else
	localtime_r(&now_tt, &now_tm);
#endif

	std::stringstream ss;
	ss << std::put_time(&now_tm, "%Y%m%d_%H%M%S");
	string time_str = ss.str();

	string path_secret = path_dir + "/secret_keys_" + time_str + ".txt";
    std::ofstream secretKeyFile(path_secret, std::ios::out | std::ios::trunc);
    if (!secretKeyFile) {
      std::cerr << "Error opening secret key output file!" << std::endl;
      return;
    }
    secretKeyFile << "# Secret key log\n";
    secretKeyFile << "# folder=" << folder << "\n";
    secretKeyFile << "# p=" << p << ", r=" << r << ", L=" << L
                  << ", c=" << c << ", t=" << skHwt
                  << ", skHwt_sp=" << skHwt_sp << "\n";

	t = -GetTime();
	if (!noPrint) cout << "Generating keys, " << std::flush;
	SecKey secretKey(context);
	TeeBuf teeBuf(noPrint ? nullptr : std::cout.rdbuf(), secretKeyFile.rdbuf());
	std::streambuf* coutBufBackup = std::cout.rdbuf(&teeBuf);
	secretKey.GenSecKey(); 
	addSome1DMatrices(secretKey); // compute key-switching matrices that we need
	addFrbMatrices(secretKey);
	if (!noPrint) cout << "computing key-dependent tables..." << std::flush;

	// Generate a sparse bootstrapping key for BTS
	if (!noPrint)
		std::cout << "\nBTS sparse secretkey generation start..." << std::endl;
	secretKey.genRecryptData(skHwt_sp);
	std::cout.rdbuf(coutBufBackup);
	secretKeyFile << "# Secret key generation done\n";
    secretKeyFile.flush();

	t += GetTime();
	if (!noPrint) cout << " done in "<<t<<" seconds\n";

	PubKey publicKey = secretKey;

	long d = context.getOrdP();
	long nslots = context.getPhiM() / d;

	// GG defines the plaintext space Z_p[X]/GG(X)
	NTL::ZZX GG;
	GG = context.getAlMod().getFactorsOverZZ()[0];
	std::shared_ptr<helib::EncryptedArray> ea(std::make_shared<helib::EncryptedArray>(context, GG));
	helib::setupDebugGlobals(&secretKey, ea);

	NTL::zz_p::init(p2r);
	NTL::Vec<NTL::zz_p> val0(NTL::INIT_SIZE, nslots);
	for (auto& x : val0)
		x = to_zz_p(0);

	std::ostringstream coeffName;
	coeffName << "coeff_" << std::put_time(&now_tm, "%Y%m%d_%H%M%S") << ".txt";

    string path_coeff = path_dir + "/" + coeffName.str();
	std::ofstream coeffFile(path_coeff, std::ios::app);
	
    if (!fileExists(path_coeff)) {
        std::ofstream outFile(path_coeff);
        outFile.close();
    }

	if (!coeffFile) {
	  std::cerr << "Error opening one of the files!" << std::endl;
	  return;
	}
	
	
	if (waitForInput) {
		std::cout << "To continue, press Enter...\n" << std::flush;
		std::cin.get();
	}

	if (progress) {
		std::cerr << "[phase] Building secret permutation matrix... elapsed=" << elapsedString() << std::endl;
	}
	
    Mat<Vec<long>> secret_perm_matrix_thin;
    secret_perm_matrix_thin.SetDims(nslots, phim);
	
    //denote s_ind => (ind +1) & -s_ind => (-ind-1) for 0 <= ind < m
    Matrix_secret_thin(secret_perm_matrix_thin, mvec, mhatvec, m_prodvec, phim_prodvec);

	if (progress) {
		std::cerr << "[phase] Building ciphertext permutation matrix... elapsed=" << elapsedString() << std::endl;
	}

    Mat<Vec<long>> ctxt_perm_matrix_thin;
    Mat<long> matrix_irwin_hall_thin;
    ctxt_perm_matrix_thin.SetDims(nslots, m);
    matrix_irwin_hall_thin.SetDims(nslots, m);

    Matrix_ctxt_thin(ctxt_perm_matrix_thin, matrix_irwin_hall_thin, secret_perm_matrix_thin);

    string path_irwin_hall = path_dir + "/irwin_hall_" + to_string(m) + ".txt";
	if (!noPrint) cerr << "Save matrix! ...";
    if (!fileExists(path_irwin_hall)) {
        saveMat(matrix_irwin_hall_thin, path_irwin_hall);
    }
	if (!noPrint) cerr << "done!" << endl;

	if (progress) {
		std::cerr << "[phase] Entering attack iterations... elapsed=" << elapsedString() << std::endl;
	}
    
    xdouble q = power(to_xdouble(p), context.getRcData().e)+1;
	long lastProgressBucket = -1;
	long lastFailProgressBucket = -1;
	long totalFailingIterations = 0;
	long totalFailSlots = 0;
	long totalCollectedRows = 0;

	for (long w=0; w< iter; ++w){
		if (progress && stopAfterFailCount <= 0 && iter > 0 && progressStep > 0) {
			const long pct = ((w + 1) * 100) / iter;
			const long bucket = pct / progressStep;
			if (bucket != lastProgressBucket || (w + 1) == iter) {
				lastProgressBucket = bucket;
				std::cerr << "[progress] " << pct << "% (" << (w + 1) << "/" << iter
						  << "), elapsed=" << elapsedString() << std::endl;
			}
		}

		if (!noPrint) {
			cerr << endl;
			cerr << "==================== " << w << " =====================" << endl;
		}

		std::vector<NTL::ZZX> val1;
		val1.resize(nslots);
		for (long i = 0; i < nslots; i++) {
			val1[i] = NTL::conv<NTL::ZZX>(NTL::conv<NTL::ZZ>(rep(val0[i])));
		}

		helib::Ctxt c1(publicKey);
		ea->encrypt(c1, publicKey, val1);

		c1.multiplyBy(c1);
		// cout << "Capacity = " << c1.bitCapacity() << endl;

		helib::Ctxt c2(c1);
		
		const DoubleCRT& recryptKey = secretKey.getRecryptKey();

		vector<Vec<ZZ>> pwrflParts; 
		NTL::Vec<NTL::xdouble> scaled_ptxt;
		vector<Vec<xdouble>> round_error;
		vector<Vec<ZZ>> v_pwfl;

		const auto t_rc0 = std::chrono::steady_clock::now();
		publicKey.thinReCryptAtk(c2, recryptKey, pwrflParts, scaled_ptxt, round_error, v_pwfl);
		t_recrypt += std::chrono::duration<double>(std::chrono::steady_clock::now() - t_rc0).count();

		std::vector<NTL::ZZX> val2;
		ea->decrypt(c2, secretKey, val2);
		
		// cout << "\nval1 = \n";
		// for (int j=0; j<val1.size(); ++j){
		//   cout << val1[j] << ", ";
		// }
		// cout << "\nval2 = \n";
		// for (int j=0; j<val2.size(); ++j){
		//   cout << val2[j] << ", ";
		// }
		// cout << endl;

		if (val1 == val2) {
			if (!noPrint)
				cout << "GOOD" << endl;
		}
		else{
			totalFailingIterations++;
			if (progress && stopAfterFailCount > 0) {
				if (progressStep > 0) {
					const long failPct = (totalFailingIterations * 100) / stopAfterFailCount;
					const long failBucket = failPct / progressStep;
					if (failBucket != lastFailProgressBucket || totalFailingIterations == stopAfterFailCount) {
						lastFailProgressBucket = failBucket;
						std::cerr << "[progress] fail: " << totalFailingIterations
								  << " / " << stopAfterFailCount
								  << " (iter " << (w + 1) << "/" << iter
								  << "), elapsed=" << elapsedString() << std::endl;
					}
				} else {
					std::cerr << "[progress] fail: " << totalFailingIterations
							  << " / " << stopAfterFailCount
							  << " (iter " << (w + 1) << "/" << iter
							  << "), elapsed=" << elapsedString() << std::endl;
				}
			}
			if (!noPrint)
				cout << "BAD! BTS failure indeed happens at slots:" << endl;

			// Indices where val1[j] != val2[j]
			Vec<long> fail_indices;
			for (long j = 0; j < val1.size(); ++j) {
			  if (val1[j] != val2[j]) {
				fail_indices.append(j);
			  }
			}
			totalFailSlots += fail_indices.length();
			
			long vecLength = v_pwfl[1].length();

			Vec<xdouble> pwfl_coeff0, pwfl_coeff1;
			pwfl_coeff0.SetLength(vecLength);
			pwfl_coeff1.SetLength(vecLength);

            // Version 1
			// for (long j=0; j < vecLength; ++j){
			// 	pwfl_coeff0[j] = to_xdouble(v_pwfl[0][j]);
			// 	pwfl_coeff1[j] = to_xdouble(v_pwfl[1][j]);
			// }

            // Version 2
			for (long j=0; j < vecLength; ++j){
				pwfl_coeff0[j] = to_xdouble(pwrflParts[0][j])/q + round_error[0][j]*(1 - 1/q) + to_xdouble(v_pwfl[0][j]);
				pwfl_coeff1[j] = to_xdouble(pwrflParts[1][j])/q + round_error[1][j]*(1 - 1/q) + to_xdouble(v_pwfl[1][j]);
			}

			for (int j=0; j < fail_indices.length(); ++j){

				long ind = fail_indices[j];
				if (!noPrint)
					cerr << ind << ", ";

				Vec<xdouble> result;
				append(result, xdouble(ind));
				append(result, pwfl_coeff0[d*ind]);
				
				Vec<xdouble> rot_coeff;
				rot_coeff.SetLength(m);
				const auto t_gc0 = std::chrono::steady_clock::now();
				getCoeff_thin(rot_coeff, ind, pwfl_coeff1, ctxt_perm_matrix_thin);
				t_getcoeff += std::chrono::duration<double>(std::chrono::steady_clock::now() - t_gc0).count();
            	result.append(rot_coeff);
				saveVec(Vec<xdouble>(result), path_coeff);
				totalCollectedRows++;
			}

			if (!noPrint && fail_indices.length() > 0)
				cerr << endl;

			if (stopAfterFirstFailure) {
				std::cout << "[stop] Stopping after first observed failure at iteration "
						  << w << std::endl;
				break;
			}

			if (stopAfterFailCount > 0 && totalFailingIterations >= stopAfterFailCount) {
				std::cout << "[stop] Reached target failing iterations: "
						  << totalFailingIterations << " at iteration " << w << std::endl;
				break;
			}

		}

	}

	if (!noPrint || progress) {
		std::cout << "[summary] failing iterations: " << totalFailingIterations << "/" << iter << std::endl;
		std::cout << "[summary] total fail slots collected: " << totalFailSlots << std::endl;
		std::cout << "[summary] total coeff rows written: " << totalCollectedRows << std::endl;
		std::cout << "[summary] time in thinReCryptAtk: " << t_recrypt
		          << " s, in getCoeff_thin: " << t_getcoeff << " s" << std::endl;
	}

	if (progress) {
		std::cerr << "[done] Attack run finished. total elapsed=" << elapsedString() << std::endl;
	}
}

// extern long fhe_disable_intFactor;
// extern long fhe_disable_fat_boot;
// extern long fhe_force_chen_han;

int main(int argc, char *argv[]){
  ArgMap amap;

  // Param 1 from HElib
  long p=2;
  long r=1;
  long c=3;
  long L=400;
  long N=0;
  long t=32;
  long nthreads=1;
  long skHwt_sp=32;
  int iter = 1000;

  long seed=0;
  long useCache=1;

  amap.arg("p", p, "plaintext base");

  amap.arg("r", r,  "exponent");
  amap.note("p^r is the plaintext-space modulus");

  amap.arg("c", c, "number of columns in the key-switching matrices");
  amap.arg("L", L, "# of levels in the modulus chain");
  amap.arg("N", N, "lower-bound on phi(m)");
  amap.arg("t", t, "Hamming weight of recryption secret key", "heuristic");
  amap.arg("dry", dry, "dry=1 for a dry-run");
  amap.arg("nthreads", nthreads, "number of threads");
	amap.arg("iter", iter, "number of attack iterations");
	amap.arg("skHwt_sp", skHwt_sp, "Hamming weight of sparse bootstrapping key");
  amap.arg("seed", seed, "random number seed");
  amap.arg("noPrint", noPrint, "suppress printouts");
  amap.arg("useCache", useCache, "0: zzX cache, 1: DCRT cache");
	amap.arg("waitForInput", waitForInput, "pause for Enter before matrix generation");
	amap.arg("progress", progress, "print phase and percent progress");
	amap.arg("progressStep", progressStep, "progress print step in percent");
	amap.arg("stopAfterFirstFailure", stopAfterFirstFailure, "stop immediately after the first failed iteration");
	amap.arg("stopAfterFailCount", stopAfterFailCount, "stop after this many failing iterations have been observed (0 disables)");

  amap.arg("force_bsgs", fhe_test_force_bsgs);
  amap.arg("force_hoist", fhe_test_force_hoist);

  //  amap.arg("disable_intFactor", fhe_disable_intFactor);
  amap.arg("chen_han", fhe_force_chen_han);

  amap.arg("debug", debug, "generate debugging output");
  amap.arg("scale", scale, "scale parameter");

  amap.parse(argc, argv);

  if (seed)
    SetSeed(ZZ(seed));

	SetNumThreads(nthreads);

  for (long i=0; i<(long)num_mValues; i++)
    if (mValues[i][0]==p && mValues[i][1]>=N) {
      TestIt(i,p,r,L,c,t,skHwt_sp, iter, useCache);
      break;
    }

  return 0;

}

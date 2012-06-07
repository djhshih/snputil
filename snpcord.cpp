#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdexcept>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>

#include <boost/bind.hpp>

std::string progname = "snpcord";



namespace chp {

	using namespace std;
	namespace qi = boost::spirit::qi;
	namespace ascii = boost::spirit::ascii;
	namespace ph = boost::phoenix;

	/** SNP call reader. Parses SNP call files. */
	class Reader {

	public:

		/** Constructor.
		 *  \param capacity initial capacity for data vectors
		 */
		Reader(int capacity=0) {
			if (capacity > 0) {
				markers.reserve(capacity);
			}
		}

		~Reader() {}

		void read(const string& filename) {
			ifstream file;
			file.open(filename.c_str(), ios::in);
			if (!file.is_open()) throw runtime_error("Failed to open input file.");
			read(file);
			file.close();
		}

		void read(ifstream& file) {
			string line;
			const char delim = '\t';

			// read header row
			while (true) {
				getline(file, line);
				if (file.eof()) break;

				// ignore comment line
				if (line[0] != '#') {
					// treat first non-comment line as header
					parse_header(line);
					break;
				}
			}

			init_calls();

			// read in data
			while (true) {
				getline(file, line);
				if (file.eof()) break;

				// ignore comment line
				if (line[0] == '#') continue;

				// first element is marker name
				size_t k = find_nth_of(line, delim);
				string s(line.begin(), line.begin() + k);
				markers.push_back(s);

				// reset sample index
				sample_idx = 0;	
				// parse data in remaining part of line
				parse_calls(line.begin() + k + 1, line.end());
			}

		}

		void print() {
			write(std::cout);
		}

		void write(const string& filename) {
			ofstream file;
			file.open(filename.c_str(), ios::out);
			if (!file.is_open()) throw runtime_error("Faield to open output file.");
			write(file);
			file.close();
		}

		/*
		void write(ostream& file) {
			for (size_t i = 0; i < samples.size(); ++i) {
				// report number of valid calls
				size_t n_valid_calls = 0;
				for (size_t j = 0; j < calls[i].size(); ++j) {
					if (calls[i][j] != -1) ++n_valid_calls;
				}
				file << samples[i] << ": " << n_valid_calls << " calls" << endl;
				// print calls
				for (size_t j = 0; j < calls[i].size(); ++j) {
					file << calls[i][j] << ' ';
				}
				file << endl;
			}
		}
		*/

		/**
		 * Output SNP genotype concordance rates
		 */
		void write(ostream& file) {
			calc_concordance();
			size_t n = samples.size();

			// print header row
			for (size_t i = 0; i < n; ++i) {
				file << samples[i] << delim;
			}
			file << endl;

			for (size_t i = 0; i < n; ++i) {
				// print header column
				file << samples[i] << delim;
				// print data
				for (size_t j = 0; j < n; ++j) {
					file << concordances[i][j] << delim;
				}
				file << endl;
			}
		}

	private:

		bool parse_header(const string& header) {
			size_t k, kk;
			// skip the first element
			k = find_nth_of(header, delim) + 1;
			while (k != -1 && k < header.size()) {
				kk = find_nth_of(header, delim, 1, k);
				if (kk == -1) {
					kk = header.size();
				}
				string s(header.begin() + k, header.begin() + kk);
				samples.push_back(s);
				k = kk + 1;
			}
		}

		template <typename iterator>
		bool parse_calls(iterator first, iterator last) {
			using qi::short_;
			
			bool r = qi::phrase_parse(first, last,
				(
					short_[boost::bind(&Reader::push_back_call, this, _1)] >>
					*( short_[boost::bind(&Reader::push_back_call, this, _1)] )
				),
				ascii::space
			);

			if (first != last) // partial match
				return false;
			return r;
		}

		void push_back_call(short const& x) {
			calls[sample_idx++].push_back(x);
		}

		// find the nth instance of character c in string str
		size_t find_nth_of(const string& str, char c, size_t n = 1, size_t pos = 0) {
			size_t k = pos-1;
			for (size_t i = 0; i < n; ++i) {
				k = str.find_first_of(c, k+1);
			}
			return k;
		}

		// initialize calls vectors
		void init_calls() {
			size_t n = samples.size();
			calls.resize(n);
			for (size_t i = 0; i < n; ++i) {
				// reserve space
				calls[i].reserve(markers.capacity());
			}
		}

		void calc_concordance() {
			size_t n = samples.size();
			// allocate space
			concordances.resize(n);
			for (size_t i = 0; i < n; ++i) {
				concordances[i].resize(n);
			}

			// calculate concordance
			for (size_t i = 0; i < n; ++i) {
				// concordance with itself is 1
				concordances[i][i] = 1;
				for (size_t j = i+1; j < n; ++j) {
					float r = concordance(calls[i], calls[j]);
					concordances[i][j] = concordances[j][i] = r;
				}
			}
		}

		typedef vector<short> genotype;

		// Return the concordance rate.
		float concordance(genotype x, genotype y) {

			if (x.size() != y.size()) {
				throw invalid_argument("Genotype vector sizes must match");
			}

			size_t n = x.size();
			size_t nvalids = 0, nmatches = 0;
			for (size_t i = 0; i < n; ++i) {
				// ignore SNP if either sample has missing call	
				if (x[i] == -1 || y[i] == -1) continue;
				
				if (x[i] == y[i]) {
					// SNP genotype match
					++nmatches;
				}
				++nvalids;
			}

			return float(nmatches) / nvalids;
		}

		vector<string> markers;
		vector<string> samples;
		vector<genotype> calls;
		vector< vector<float> > concordances;

		size_t sample_idx;

		static const char delim = '\t';

	};

}



int main(int argc, char *argv[]) {

	using namespace std;

	po::variables_map vm;
	po::options_description opts;
	po::positional_options_description popts;

	string input_fname, output_fname;
	size_t capacity;

	// setup command line options
	
	opts.add_options()
		("help,h", "print help message")
		("input,i", po::value<string>(&input_fname), "input SNP calls file")
		("output,o", po::value<string>(&output_fname), "output file")
		("capacity,c", po::value<size_t>(&capacity)->default_value(1e6), "expected length of the data arrays")
	;
	popts.add("input", 1);

	po::store(
		po::command_line_parser(argc, argv)
			.options(opts)
			.positional(popts)
			.run(),
		vm
	);
	po::notify(vm);

	// process command linen options
	
	if (vm.count("help")) {
		cout << "usage: " << progname << " [options] <input file>" << endl;
		cout << opts << endl;
		return 0;
	}

	if (!vm.count("input")) {
		throw invalid_argument("No input file specified.");
	}

	// read and write data
	chp::Reader reader(capacity);
	reader.read(input_fname);

	if (output_fname.size() > 0) {
		reader.write(output_fname);
	} else {
		reader.print();
	}

	return 0;
}


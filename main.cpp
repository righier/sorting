#include "utils.h"
#include "args.h"
#include "sort.h"
#include "timer.h"

using namespace std;

int main(int argc, char **argv) {

	Args args(argc, argv);

	int N = args.getInt("-n", 100000000);
	string S = args.getString("--algo", "insertion");
	int nTests = args.getInt("--tests", 10);
	bool printArray = args.getBool("--print");
	bool verbose = args.getBool("--verbose");
	int threshold = args.getInt("--threshold", 50);
	int buckets = args.getInt("--buckets", 256);

	cout << "N: " << N << endl;
	cout << "num tests: " << nTests << endl;
	cout << "algo: " << S << endl;
	cout << "threshold: " << threshold << endl;
	cout << "buckets: " << buckets << endl;
	cout << "----------------------------" << endl << endl;

	Random rng(31);

	Sort *sorter = nullptr;

	if (S == "gcc") {
		sorter = new GccSort();
	} else if (S == "insertion") {
		sorter = new InsertionSort();
	} else if (S == "quick") {
		sorter = new QuickSort();
	} else if (S == "hybrid") {
		sorter = new HybridSort(threshold);
	} else if (S == "bucket") {
		sorter = new BucketSort(buckets, threshold);
	} else if (S == "adabucket") {
		sorter = new AdaptiveBucketSort(threshold);
	} else if (S == "bucket2") {
		sorter = new BucketSort2(buckets, threshold);
	} else {
		cerr << "UNKNOWN SORT ALGORITHM" << endl;
		return 1;
	}

	double totalTime = 0;
	for (int i = 1; i <= nTests; i++) {

		vector<u64> v = randomVector(rng, N);

		if (verbose){
			cout << "generated: " << i << endl;
		}

		if (printArray) {
			cout << v << endl;
		}

		Timer<> timer;
		timer.start();

		u64 output = sorter->sort(v);

		totalTime += timer.delta();

		cout << output << endl;

		if (printArray) {
			cout << v << endl;
		}

		if (verbose) {
			cout << "done: " << i << endl << endl;
		}

	}

	cout << "average time: " << (totalTime / nTests) * 1000 << endl;

	delete sorter;

	return 0;
}